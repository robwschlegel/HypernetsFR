# code/0_functions.R
# Code shared across the project. Sourced by every other script; never run standalone.


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(DBI)
library(RSQLite)
library(ncdf4)
library(tidync)
library(FNN) # Needed for fastest nearest neighbor searching
library(fuzzyjoin) # For joining data based on nearest neighbor searching
library(geosphere) # For determining distance between points
library(ggtext) # For rich text labels
library(ggimage) # For adding .jpg files to figures
library(patchwork) # For complex paneling of figures
library(future)
library(furrr)
library(ggh4x) # For independent per-facet axis limits (facetted_pos_scales()) alongside a true
                # 1:1 x/y aspect ratio per panel -- see plot_global_nm())
library(sf) # For the clean_water_sf point-in-polygon pixel filter -- see db_export_matchups_site()


# Setup -------------------------------------------------------------------

# Disable scientific notation
# NB: Necessary for correct time stamp conversion
options(scipen = 9999)

plan(multicore, workers = max(1L, parallel::detectCores() - 2L))


# Utilities ---------------------------------------------------------------

# Define the wavelength (nm) band colour palette
colour_nm_func <- function(sensor_Y){
  if(sensor_Y == "PACE"){
    labels_nm <- c("351-400", "401-450", "451-500", "501-550", "551-600", "601-650", "651-700", "701-750", "751-800", "801-900", "901-1050")
    colour_nm <- c("darkviolet", "violet", "blue", "darkgreen", "yellow", "orange", "red", "firebrick", "sienna", "black", "#777777")
  } else if(sensor_Y == "AQUA"){
    labels_nm <- c("412","443","469","488","531","547","555","645","667","678")
    colour_nm <- c("darkviolet","blue","cyan","green","yellowgreen","yellow","orange","red","firebrick","sienna")
  } else if(sensor_Y == "SNPP"){
    labels_nm <- c("410","443","486","551","671")
    colour_nm <- c("darkviolet","blue","cyan","yellowgreen","red")
  } else if(sensor_Y %in% c("JPSS1", "JPSS2")){
    labels_nm <- c("411","445","489","556","667")
    colour_nm <- c("darkviolet","blue","cyan","yellowgreen","red")
  } else if(sensor_Y %in% c("S3A", "S3B", "S3", "S3_all")){
    labels_nm <- c("400", "412", "442", "490", "510", "560", "620", "665", "673", "681", "709", "754", "779", "865", "885", "1020")
    colour_nm <- c("darkviolet","blueviolet","blue4","blue2",
                   "#00BFFF","#00FF7F","#ADFF2F","#FFFF00",
                   "#FFD700","#FFA500","#FF8C00","#FF4500",
                   "#FF0000","#8B0000","#4b0000","#777777")
  } else {
    stop(paste0("Incorrect value for 'sensor_Y' : ",sensor_Y))
  }
  names(colour_nm) <- labels_nm
  return(colour_nm)
}


# Function that assembles file directory based on desired variable and sensors
file_path_build <- function(site_name, sat_name){
  file_path <- paste0("~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/", site_name, "/RHOW_HYPERNETS_vs_", sat_name)
}

# Inverse of file_path_build(): recover a matchup file's site_name from its full path
# (".../FR/<site_name>/RHOW_HYPERNETS_vs_<sat_name>/<file>.csv"). Needed because raw matchup
# filenames are not unique *across* sites -- e.g. the same satellite pass can produce an
# identically-named export at both MAFR and THFR when their HYPERNETS scan schedules happen to
# align on the same 15-min slots (confirmed: 43 such collisions for OLCI alone). Any code that
# joins/matches on file_name across a multi-site file list (see code/2_outliers.R) must include
# site_name in the key, or it silently mixes rows from different sites together.
path_site_name <- function(path){
  basename(dirname(dirname(path)))
}

# Load a single matchup file and create mean values from all replicates
# file_name <- "/home/calanus/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/MAFR/RHOW_HYPERNETS_vs_JPSS1/JPSS1_20240531T120600_vs_HYPERNETS_20240531T114500_RHOW_C.csv"
# file_name <- "/home/calanus/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/MAFR/RHOW_HYPERNETS_vs_SNPP/SNPP_20240611T131200_vs_HYPERNETS_20240611T124500_RHOW_C.csv"
# file_name <- file_path
load_matchup_mean <- function(file_name){
  
  # message(paste0("Started loading : ", file_name))
  
  # Load the csv file
  suppressMessages(
    df_match <- read_delim(file_name, delim = ";", col_types = "ccccnnic")
  )
  colnames(df_match)[1] <- "sensor"
  
  # Get means per file
  # NB: Satellite matchups have a different structure than in situ matchups
  # NB: For the moment, JPSS1 files have 'rhow weighted' not 'weighted', so this is accounted for below
  # This needs to be fixed in Hypernets_matchups
  if(any(df_match$data_type %in% c("weighted", "rhow weighted", "computed weighted"))){
    df_mean <- df_match |>
      mutate(sensor = gsub(" 1$| 2$| 3$| 4$| 5$| 6$| 7$| 8$| 9$", "", sensor)) |>
      filter(sensor != "Hyp_nosc") |> 
      filter(data_type %in% c("weighted", "rhow weighted", "computed weighted"))
  } else {
    df_mean <- df_match |> 
      filter(grepl(" 1", sensor)) |> 
      mutate(sensor = gsub(" 1$| 2$| 3$| 4$| 5$| 6$| 7$| 8$| 9$", "", sensor)) |>
      filter(sensor != "Hyp_nosc")
  }
  
  # Remove unneeded columns
  df_mean <- df_mean |> 
      dplyr::select(-radiometer_id, -data_type, -type, -pixel_pos, -variability_centered)

  # Double check that only two rows of data have been selected
  if(nrow(df_mean) > 2){
    stop(paste0("More than 2 rows returned for : ", file_name))
  }
  if(nrow(df_mean) < 2){
    cat(paste0("Less than 2 rows returned for : ", file_name,",\n defaulting to unweighted in situ data for now..."))
    df_mean <- df_match |>
      mutate(sensor = gsub(" 1$| 2$| 3$| 4$| 5$| 6$| 7$| 8$| 9$", "", sensor)) |>
      filter(sensor != "Hyp_nosc") |> 
      filter(data_type %in% c("computed", "computed weighted")) |> 
      dplyr::select(-radiometer_id, -data_type, -type, -pixel_pos, -variability_centered)
  }

  # Exit
  # message(paste0("Finished loading : ", file_name))
  return(df_mean)
}

# Load a single matchup file directly into long format
# file_name <- file.path(folder_path, file_uniq_list$file_name)[1]
# file_name <- "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/MAFR/RHOW_HYPERNETS_vs_S3A/S3A_20240531T101658_vs_HYPERNETS_20240531T100000_RHOW.csv"
# file_name <- "/home/calanus/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/MAFR/RHOW_HYPERNETS_vs_JPSS1/JPSS1_20240531T120600_vs_HYPERNETS_20240531T114500_RHOW_C.csv"
# file_name <- file_list_clean[[1]]
load_matchup_long <- function(file_name){
  
  df_mean <- load_matchup_mean(file_name)
  
  # Pivot longer
  df_long <- df_mean |> 
    pivot_longer(cols = matches("1|2|3|4|5|6|7|8|9"), names_to = "wavelength", values_to = "value") |> 
    dplyr::select(-day, -time, -longitude, -latitude) |>
    # filter(value <= 1) |>  # Remove erroneously high values
    pivot_wider(names_from = sensor, values_from = value) |>
    na.omit() |> 
    mutate(wavelength = as.numeric(wavelength),
           file_name = basename(file_name), .before = "wavelength")
  
  # Exit
  return(df_long)
}

# Load the variance data for a matchup file
load_matchup_var <- function(file_name){

  # Load the csv file
  suppressMessages(
    df_match <- read_delim(file_name, delim = ";", col_types = "ccccnnic")
  )
  colnames(df_match)[1] <- "sensor"
  
  # Calculate uncertainties per wavelength
  df_var <- df_match |> 
    mutate(sensor = gsub(" 1$| 2$| 3$| 4$| 5$| 6$| 7$| 8$| 9$", "", sensor),
           var_name = "rhow") |>
    filter(!(sensor %in% c("Hyp_nosc"))) |> 
    # distinct() |> 
    dplyr::select(-day, -time, -latitude, -longitude, -radiometer_id, -type) |> 
    mutate(data_type = case_when(data_type == "rhow" ~ "var_value", TRUE ~ data_type)) |> 
    pivot_longer(cols = matches("1|2|3|4|5|6|7|8|9"), names_to = "wavelength", values_to = "value") |> 
    pivot_wider(names_from = data_type, values_from = value) |>
    # na.omit() |> 
    mutate(max_sd_diff = abs(var_value - std_max),
           min_sd_diff = abs(var_value - std_min),
           # Max and min should be the same, but this addresses any rounding issues
           sd = (max_sd_diff + min_sd_diff)/2,
           cv = sd/var_value,
           wavelength = as.numeric(wavelength))
}

# Load all files in a given folder
load_matchups_folder <- function(site_name, sat_name, long = FALSE){
  
  # Create file path
  folder_path <- file_path_build(site_name, sat_name)
  
  # List all files in directory
  file_list <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)
  
  # Remove stats output files
  file_list_clean <- file_list[!grepl("all|global", file_list)]
  
  # Load data
  if(long){
    match_base <- furrr::future_map_dfr(file_list_clean, load_matchup_long, .options = furrr_options(seed = TRUE)) |> 
       mutate(site_name = site_name, .before = wavelength)
  } else {
    match_base <- furrr::future_map_dfr(file_list_clean, load_matchup_mean, .options = furrr_options(seed = TRUE)) |> 
       mutate(site_name = site_name, .before = sensor)
  }

  # Exit
  return(match_base)
}

# Convenience function to get lon/lat coords from HYPERNETS .nc files
load_HYPERNETS_coords <- function(file_name){
  ncdump::NetCDF(file_name)$attribute$global[c("site_longitude", "site_latitude")]
}

# Process HYPERNETS files to get mean and sd per wavelength per sequence
# file_name <- L1C_HYPERNETS_files[2]
proc_HYPERNETS_L1C <- function(file_name, stat_calc = TRUE){

  print(file_name)

  # Get file info for use later
  suppressMessages(nc_info <- ncdump::NetCDF(file_name))
  
  # Load observation values
  df_meta <- tidync(file_name) |> 
    activate("D1") |> 
    hyper_tibble()

  # Loda measurement variables
  df_base <- tidync(file_name) |> 
    hyper_tibble() |> 
    mutate(Rrs = reflectance / pi, .before = "downwelling_radiance") |> 
    left_join(df_meta, by = "scan") |> 
    dplyr::select(wavelength, scan, rhof_wind:rhof_vza, Rrs:reflectance_nosc)
  
  # Calculate stats if desired
  if(stat_calc){

    df_res <- df_base |> 
      pivot_longer(rhof_wind:reflectance_nosc) |> 
      mutate(wavelength = round(as.numeric(wavelength))) |> 
      summarise(n = n(),
                mean = mean(value, na.rm = TRUE),
                sd = sd(value, na.rm = TRUE), .by = c("wavelength", "name")) |> #, "scan")) |> 
      # filter(name != "reflectance_nosc") |> 
      mutate(cv = sd / mean,
            lon = nc_info$attribute$global$site_longitude,
            lat = nc_info$attribute$global$site_latitude,
            date = as.POSIXct(sub("SEQ", "", nc_info$attribute$global$source_file),
              format = "%Y%m%dT%H%M%S", tz = "UTC"),
            system = "HYPERNETS") |> 
      mutate(name = case_when(name == "reflectance" ~ "Rhow",
                              name == "reflectance_nosc" ~ "Rhow_nosc",
                              name == "water_leaving_radiance" ~ "Lw",
                              name == "downwelling_radiance" ~ "Ld",
                              name == "upwelling_radiance" ~ "Lu",
                              name == "irradiance" ~ "Ed",
                              TRUE ~ name)) |> 
      dplyr::select(system, lon, lat, date, name, wavelength, n, mean, sd, cv)

  } else {

    df_res <- df_base |> 
      pivot_longer(rhof_wind:reflectance_nosc) |> 
      mutate(wavelength = round(as.numeric(wavelength)),
             lon = nc_info$attribute$global$site_longitude,
             lat = nc_info$attribute$global$site_latitude,
             date = as.POSIXct(sub("SEQ", "", nc_info$attribute$global$source_file),
               format = "%Y%m%dT%H%M%S", tz = "UTC"),
            system = "HYPERNETS") |> 
      mutate(name = case_when(name == "reflectance" ~ "Rhow",
                              name == "reflectance_nosc" ~ "Rhow_nosc",
                              name == "water_leaving_radiance" ~ "Lw",
                              name == "downwelling_radiance" ~ "Ld",
                              name == "upwelling_radiance" ~ "Lu",
                              name == "irradiance" ~ "Ed",
                              TRUE ~ name)) |> 
      dplyr::select(system, lon, lat, date, name, wavelength, value)
  }
  
  # Exit
  return(df_res)
}

# Check the amount of variance in satellite files and return a message if there is an issue
# file_name <- "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/MAFR/RHOW_HYPERNETS_vs_S3A/S3A_20240531T101658_vs_HYPERNETS_20240531T100000_RHOW.csv"
sat_var_check <- function(file_name, cv_limit = 30){
  
  # Load the csv file
  suppressMessages(
    df_match <- read_delim(file_name, delim = ";", col_types = "ccccnnic")
  )
  colnames(df_match)[1] <- "sensor"
  
  # PACE files don't have weighted mean values
  if(!("weighted" %in% df_match$data_type)){
    df_match <- df_match |> 
      mutate(data_type = case_when(data_type == "rhow" ~ "weighted", TRUE ~ data_type))
  }
  
  # Old legacy code
  # Check for variance in W_nm columns
  # df_check <- df_match |> 
  #   mutate(sensor = gsub(" 1$| 2$| 3$| 4$| 5$| 6$| 7$| 8$| 9$", "", sensor)) |>
  #   filter(!(sensor %in% c("Hyp_nosc", "Hyp"))) |>
  #   dplyr::select(-day, -time, -latitude, -longitude, -radiometer_id, -type, -pixel_pos) |> 
  #   pivot_longer(cols = matches("1|2|3|4|5|6|7|8|9"), names_to = "wavelength", values_to = "value") |> 
  #   pivot_wider(names_from = data_type, values_from = value) |>
  #   na.omit() |> 
  #   mutate(max_sd_diff = abs(weighted - std_max),
  #          min_sd_diff = abs(weighted - std_min),
  #          # Max and min should be the same, but this addresses any rounding issues
  #          sd = (max_sd_diff + min_sd_diff),
  #          cv = sd/weighted,
  #          wavelength = as.numeric(wavelength)) |> 
  #   filter(wavelength <= 600, wavelength >= 400)
  
  # Old legacy code
  # If variables are too different, issue a warning and omit file from being loaded
  # if(nrow(df_check) > 0){
  #   if(mean(df_check$variability_centered, na.rm = TRUE) >= cv_limit){
  #     # warning(paste0("Weighted mean has too much variance in file : ", file_name))
  #     # return(basename(file_name))
  #     file_check <- basename(file_name)
  #   } else {
  #     file_check <- NULL
  #   }
  # } else {
  #   file_check <- NULL
  # }

  # Get the existing vcariance column
  df_check <- df_match |> 
    dplyr::select(variability_centered) |> 
    na.omit() |> 
    distinct()

  return(data.frame(file_name = file_name, cv = abs(df_check$variability_centered), cv_limit = cv_limit))
}

# Determine which sites currently have data on disk for a given satellite platform
# NB: MAFR (Gironde Estuary, highly turbid) and THFR (lagoon, clear-water comparison site
# to MAFR, see Doxaran et al. 2024's analogous Berre lagoon vs Gironde Estuary contrast)
# are both present on disk as of 2026-07-14. Additional sites are picked up automatically
# the moment their data folder exists on disk. 
# THFR_NE (added 2026-08-28) and THFR_poly (added 2026-08-28) are independent derived sites
# alongside MAFR/THFR -- NE-quadrant and clean-water-polygon-filtered THFR re-analyses
# respectively (see the TEMPORARY section in this file and meta/pixel_explore.R), each
# intentionally kept as its own site for direct comparison rather than replacing THFR anywhere
# in the pipeline.
available_sites <- function(sat_name){
  candidate_sites <- c("MAFR", "THFR", "THFR_NE", "THFR_poly")
  site_present <- vapply(candidate_sites, function(s) dir.exists(file_path_build(s, sat_name)), logical(1))
  sites_found <- candidate_sites[site_present]
  if(length(sites_found) == 0) stop(paste0("No site data found on disk for sensor: ", sat_name))
  return(sites_found)
}

# Site-specific matchup time-window limit (minutes), added 2026-07-10.
# NB: unlike Doxaran et al. 2024 (who additionally varied the *spatial* matchup criterion per site --
# a 3x3-pixel box at Berre vs. nearest-pixel-only at Gironde), this pipeline keeps ONE spatial rule
# (nearest pixel + a fixed dist_limit = 10 km ceiling, see process_sensor()) for every site, and
# varies only the TIME window: MAFR's fast tidal turbidity dynamics need a tighter window than THFR's
# comparatively stable lagoon water, mirroring Doxaran et al. 2024's Gironde (+/-15 min) vs Berre
# (+/-30 min) choice. See code/3_sensitivity.R for the empirical check behind these two numbers, and
# manuscript/roadmap.md for the open methodological question this resolves.
# Values checked against real MAFR and THFR data via code/3_sensitivity.R (2026-07-14).
site_diff_time_limit <- function(site_name){
  dplyr::case_when(
    site_name == "MAFR" ~ 15,
    site_name == "THFR" ~ 30,
    TRUE ~ 30 # fallback for any future/unrecognised site
  )
}

# Average multiple same-day HYPERNETS-scan-vs-satellite-overpass matchups into one representative
# daily value per wavelength, for use at the global-stats stage. Added 2026-07-10 per project
# decision: rather than simply tightening/loosening the diff_time filter per site, matchups that pass
# the site-specific site_diff_time_limit() threshold should be averaged per day before
# global per-wavelength statistics are computed, rather than treated as independent data points --
# this also directly addresses the "not all matchups are independent of one another" caveat raised in
# the Tara "in review" paper's Conclusion.
#
# Validated 2026-07-14 against real MAFR data for all four sensor families (OLCI, MODIS, VIIRS,
# OCI/PACE). Date-extraction logic reuses the exact convention used in global_scatterplot() -- i.e.
# split file_name on "_", take the 2nd element, split that on "T", take the 1st element -- confirmed
# to work for S3A/S3B/JPSS1/JPSS2/SNPP/AQUA/PACE naming patterns at both MAFR and THFR.
#
# df: expects the long-format data.frame produced by load_matchup_long()/load_matchups_folder(long =
#     TRUE), i.e. one row per file_name x wavelength, already filtered to wavelength %in% W_nm and to
#     files passing site_diff_time_limit() (that filtering happens upstream, in global_stats()).

daily_average_matchups <- function(df, site_name){
  df |>
    mutate(match_date = sapply(str_split(file_name, "_"), "[[", 2),
           match_date = sapply(str_split(match_date, "T"), "[[", 1),
           match_date = as.Date(match_date, format = "%Y%m%d")) |>
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
              n_scans_averaged = dplyr::n(),
              .by = c("match_date", "wavelength")) |>
    mutate(site_name = site_name, .before = "match_date")
}

# Create a grid of sensor to ply over
sensor_grid <- function(sensor_Z){

   # Satellite names per sensor
    if(sensor_Z == "MODIS"){
      sensor_Y <- c("AQUA")
    } else if(sensor_Z == "OCI"){
      sensor_Y <- c("PACE")
    } else if(sensor_Z == "VIIRS"){
      sensor_Y <- c("SNPP", "JPSS1", "JPSS2")
    } else if(sensor_Z == "OLCI"){
      sensor_Y <- c("S3A", "S3B")
    } else {
      stop("Incorrect name given for sensor_Z")
    }

  # Print sensors for ease of use
  message("Sensor name : ", paste0(sensor_Z, collapse = ", ")); message("Sat name(s) : ",paste0(sensor_Y, collapse = ", "))

  # Create grid for mdply()
  # NB: site list is determined per sensor_Y via available_sites() so THFR is included
  # automatically once its data exist on disk (see available_sites() above); MAFR-only today.
  site_list <- unique(unlist(lapply(sensor_Y, available_sites)))
  ply_grid <- expand_grid(site_name = site_list, sensor_Y = sensor_Y) |> distinct()
}

# Output desired wavelengths based on sensor_Y
# TODO: Increase wm range into IR
W_nm_out <- function(sensor_Y){
  if(sensor_Y == "PACE"){
    W_nm <- 350:1150
  } else if(sensor_Y == "AQUA"){
    W_nm <- c(412, 443, 469, 488, 531, 547, 555, 645, 667, 678)
  } else if(sensor_Y == "SNPP"){
    W_nm = c(410, 443, 486, 551, 671)
  } else if(sensor_Y == "JPSS1" | sensor_Y == "JPSS2"){
    W_nm <- c(411, 445, 489, 556, 667)
  } else if(sensor_Y %in% c("S3A", "S3B", "S3", "S3_all")){
    W_nm <- c(400, 412, 442, 490, 510, 560, 620, 665, 673, 681, 709, 754, 779, 865, 885, 1020)
  } else {
    stop(paste0("Incorrect value for 'sensor_Y' : ",sensor_Y))
  }
}

# Get possible MODIS files
# NB: At the moment this is optimized to work with just one day of data
MODIS_dl <- function(prod_id, dl_date, bbox, usrname, psswrd, dl_files = TRUE, dl_dir = "data/MODIS"){
  
  # If download is FALSE, just print possible files
  if(!dl_files){
    message("Data files : ")
    luna::getNASA(prod_id, dl_date, dl_date, aoi = bbox, download = FALSE)
    message("Mask files : ")
    luna::getNASA("MOD44W", dl_date, dl_date, aoi = bbox, download = FALSE)
  } else {
    message("Data files : ")
    luna::getNASA(prod_id, dl_start, dl_start, aoi = bbox, download = TRUE, overwrite = FALSE,
                  path = dl_dir, username = earth_up$usrname, password = earth_up$psswrd)
    message("Mask files : ")
    luna::getNASA("MOD44W", dl_start, dl_start, aoi = bbox, download = TRUE, overwrite = FALSE,
                  path = dl_dir, username = usrname, password = psswrd)
  }
}

# Process MODIS data in a batch
# file_names <- list.files(path = "data/MODIS", pattern = "MOD", full.names = TRUE)
# file_names <- list.files(path = "data/MODIS", pattern = "MYD", full.names = TRUE)
MODIS_proc <- function(file_names, bbox, water_mask = FALSE){
  
  # Load files with desired layers etc.
  # NB: If run in parallel, merge() causes a crash to desktop
  if(water_mask){
    data_layers <- lapply(file_names, rast, subds = 2)
    data_merge <- do.call(merge, data_layers)
    # plot(data_merge)
    data_base <- terra::ifel(data_merge %in% c(1, 2, 3, 4, 5), NA, data_merge)
    # plot(data_base)
  } else {
    data_layers <- lapply(file_names, rast, lyrs = 3) # Blue green band width; 459-479 nm
    data_base <- do.call(merge, data_layers)
    # plot(data_base)
  }
  
  # Project to EPSG:4326
  data_base_proj <- raster::project(data_base, y = "EPSG:4326")
  # plot(data_base_proj)
  
  # Crop to bbox and exit
  data_crop <- raster::crop(data_base_proj, bbox)
  # plot(data_crop)
  return(data_crop)
}


# Statistics --------------------------------------------------------------

# Basic statistic calculations
# Expects data as two vectors of equal length sampled at the same time/space
base_stats <- function(x_vec, y_vec){
  
  # Ensure values are numeric
  if(!is.numeric(x_vec)) stop("x_vec is not numeric")
  if(!is.numeric(y_vec)) stop("y_vec is not numeric")

  # Remove paired values when either side is NA
  valid_pairs <- !is.na(x_vec) & !is.na(y_vec)
  x_valid <- x_vec[valid_pairs]
  y_valid <- y_vec[valid_pairs]

  # Check for too many negative or NA values before calculating stats
  valid_idx <- (x_valid > 0) & (y_valid > 0)
  x_clean <- x_valid[valid_idx]
  y_clean <- y_valid[valid_idx]
  n_clean <- length(x_clean)

  # Return empty data.frame if too many issues
  if(n_clean < 2){
    return(data.frame(row.names = NULL,
                         n = n_clean,
                         Slope = NA,
                         Slope_log = NA,
                         Slope_II_low = NA,
                         Slope_II = NA,
                         Slope_II_high = NA,
                         Slope_II_int_low = NA,
                         Slope_II_int = NA,
                         Slope_II_int_high = NA,
                         Slope_II_p = NA,
                         Slope_II_slope_bias_sig = NA,
                         Slope_II_int_bias_sig = NA,
                         RMSE = NA,
                         MSA = NA,
                         MAPE = NA,
                         MRD_25 = NA,
                         MRD_50 = NA,
                         MRD_75 = NA,
                         MARD_25 = NA,
                         MARD_50 = NA,
                         MARD_75 = NA,
                         Bias_25 = NA,
                         Bias_50 = NA,
                         Bias_75 = NA,
                         Error_25 = NA,
                         Error_50 = NA,
                         Error_75 = NA))
  }

  # Calculate RMSE (Root Mean Square Error)
  rmse <- sqrt(mean((y_clean - x_clean)^2, na.rm = TRUE))
  
  # Calculate MAPE (Mean Absolute Percentage Error)
  mape <- mean(abs((y_clean - x_clean) / x_clean), na.rm = TRUE) * 100
  
  # Calculate MSA (Mean Squared Adjustment)
  msa <- mean(abs(y_clean - x_clean), na.rm = TRUE)
  
  # Calculate median absolute relative difference
  mard_25 <- quantile(abs(y_clean - x_clean)/x_clean, 0.25, na.rm = TRUE)
  mard_50 <- quantile(abs(y_clean - x_clean)/x_clean, 0.50, na.rm = TRUE)
  mard_75 <- quantile(abs(y_clean - x_clean)/x_clean, 0.75, na.rm = TRUE)

  # Calculate median relative difference
  mrd_25 <- quantile((y_clean - x_clean)/x_clean, 0.25, na.rm = TRUE)
  mrd_50 <- quantile((y_clean - x_clean)/x_clean, 0.50, na.rm = TRUE)
  mrd_75 <- quantile((y_clean - x_clean)/x_clean, 0.75, na.rm = TRUE)

  # Calculate Bias
  log_ratio <- log10(y_clean / x_clean)
  log_ratio_25 <- quantile(log_ratio, 0.25, na.rm = TRUE)
  log_ratio_50 <- quantile(log_ratio, 0.50, na.rm = TRUE)
  log_ratio_75 <- quantile(log_ratio, 0.75, na.rm = TRUE)
  bias_perc_25 <- 100 * (sign(log_ratio_25) * (10^abs(log_ratio_25) - 1))
  bias_perc_50 <- 100 * (sign(log_ratio_50) * (10^abs(log_ratio_50) - 1))
  bias_perc_75 <- 100 * (sign(log_ratio_75) * (10^abs(log_ratio_75) - 1))
  
  # Calculate error
  log_ratio_25_abs <- quantile(abs(log_ratio), 0.25, na.rm = TRUE)
  log_ratio_50_abs <- quantile(abs(log_ratio), 0.50, na.rm = TRUE)
  log_ratio_75_abs <- quantile(abs(log_ratio), 0.75, na.rm = TRUE)
  error_perc_25 <- 100 * (10^log_ratio_25_abs - 1)
  error_perc_50 <- 100 * (10^log_ratio_50_abs - 1)
  error_perc_75 <- 100 * (10^log_ratio_75_abs - 1)
  
  # Calculate linear slope
  lin_fit <- lm(y_clean ~ x_clean)
  slope <- coef(lin_fit)[2]
  
  # Calculate log-log linear slope
  log_lin_fit <- lm(log10(y_clean) ~ log10(x_clean))
  log_slope <- coef(log_lin_fit)[2]
  
  # Calculate Model II weighted regression slope
  # NB: To create a weghted comparison would require the data to be weighted in advance
  # E.g.: df$w_x <- 1 / df$sd_x^2 OR df$w_xy <- 1 / (df$sd_x^2 + df$sd_y^2)
  # model II regression cannt be run on data with no variance
  # This is an issue for NIR HyperPRO data where the values are always the same
  if(length(unique(round(x_clean, 8))) == 1 | length(unique(round(y_clean, 8))) == 1){
    model_II_intercept <- NA; model_II_slope <- NA; model_II_p_perm <- NA
    model_II_slope_lo <- NA; model_II_slope_hi <- NA; model_II_int_lo <- NA
    model_II_int_hi <- NA; model_II_slope_bias_sig <- NA; model_II_int_bias_sig <- NA
  } else {
    model_II_fit <- lmodel2::lmodel2(y_clean ~ x_clean,
                                     range.y = "relative", 
                                     range.x = "relative",
                                     nperm = 99)
  # print(model_II_fit)

  # Extract results for chosen method
  model_II_method_choice <- "SMA" # Symetrical Major Axis
  model_II_results <- model_II_fit$regression.results |>
    filter(Method == model_II_method_choice)
  model_II_ci <- model_II_fit$confidence.intervals |>
    filter(Method == model_II_method_choice)

  # Get specific results
  model_II_intercept <- model_II_results$Intercept
  model_II_slope <- model_II_results$Slope
  model_II_p_perm <- model_II_results$`P-perm (1-tailed)`

  model_II_slope_lo <- model_II_ci$`2.5%-Slope`
  model_II_slope_hi <- model_II_ci$`97.5%-Slope`
  model_II_int_lo <- model_II_ci$`2.5%-Intercept`
  model_II_int_hi <- model_II_ci$`97.5%-Intercept`
  
  # Determine significance
  model_II_slope_bias_sig <- model_II_slope_lo >= 1 || model_II_slope_hi <= 1
  model_II_int_bias_sig <- model_II_int_lo >= 0 || model_II_int_hi <= 0
  }
  

  # Combine int data.frame and exit
  df_stats <- data.frame(row.names = NULL,
                         n = n_clean,
                         Slope = round(slope, 4),
                         Slope_log = round(log_slope, 4),
                         Slope_II_low = round(model_II_slope_lo, 4),
                         Slope_II = round(model_II_slope, 4),
                         Slope_II_high = round(model_II_slope_hi, 4),
                         Slope_II_int_low = round(model_II_int_lo, 6),
                         Slope_II_int = round(model_II_intercept, 6),
                         Slope_II_int_high = round(model_II_int_hi, 6),
                         Slope_II_p = round(model_II_p_perm, 4),
                         Slope_II_slope_bias_sig = model_II_slope_bias_sig,
                         Slope_II_int_bias_sig = model_II_int_bias_sig,
                         RMSE = round(rmse, 6),
                         MSA = round(msa, 6),
                         MAPE = round(mape, 4),
                         MRD_25 = round(mrd_25, 6),
                         MRD_50 = round(mrd_50, 6),
                         MRD_75 = round(mrd_75, 6),
                         MARD_25 = round(mard_25, 6),
                         MARD_50 = round(mard_50, 6),
                         MARD_75 = round(mard_75, 6),
                         Bias_25 = round(bias_perc_25, 4),
                         Bias_50 = round(bias_perc_50, 4),
                         Bias_75 = round(bias_perc_75, 4),
                         Error_25 = round(error_perc_25, 4),
                         Error_50 = round(error_perc_50, 4),
                         Error_75 = round(error_perc_75, 4))
  return(df_stats)
}


# Matchup processing ------------------------------------------------------

# get n nearest pixels
get_nearest_pixels <- function(df_data, target_lat, target_lon, n_pixels){
  
  # Extract latitude and longitude into a matrix
  df_coords <- df_data[, c("latitude", "longitude")]
  
  # Target coordinate as a data.frame
  target_coord <- data.frame(latitude = target_lat, 
                             longitude = target_lon)
  
  # Find the indices of the 5 nearest neighbors
  knn_indices <- get.knnx(df_coords, target_coord, k = n_pixels)
  
  # Extract the 5 nearest rows
  df_res <- df_data[as.vector(knn_indices$nn.index), ]
  return(df_res)
}

# Function that interrogates each matchup file to produce the needed output for all following comparisons
# file_path <- "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/MAFR/RHOW_HYPERNETS_vs_S3A/S3A_20240531T101658_vs_HYPERNETS_20240531T100000_RHOW.csv"
# file_path <- file_list[28]
process_matchup_file <- function(file_path){

  # Load the mean data
  df_mean <- load_matchup_mean(file_path)
  
  # Sensors to be compared
  sensors <- unique(df_mean$sensor)
  
  # Prep empty df for 2+ sensor comparisons
  df_results <- data.frame()
  
  # The double loop
  for(i in 1:length(sensors)){
    for(j in 1:length(sensors)){
      if(sensors[j] != sensors[i]){
        
        # Get data.frame for matchup based on the two sensors being compared
        df_sensor_sub <- df_mean |> 
          filter(sensor %in% c(sensors[i], sensors[j])) |> 
          mutate(dateTime = as.POSIXct(paste(day, time), format = "%Y%m%d %H%M%S", tz = "Europe/Paris"), .before = "latitude", .keep = "unused")
        
        # get distances
        hav_dist <- round(distHaversine(df_sensor_sub[c("longitude", "latitude")])/1000, 2) # distance in km
        
        # Time differences
        time_diff <- round(as.numeric(abs(difftime(df_sensor_sub$dateTime[[1]],
                                                   df_sensor_sub$dateTime[[2]], units = "mins"))))
        
        # Melt it for additional stats
        df_sensor_long <- df_sensor_sub |> 
          # dplyr::select(!(colnames(df_sensor_sub) %in% c("data_type", "dateTime", "longitude", "latitude", "radiometer_id", 
          #   "pixel_pos", "type", "variability_centered"))) |>
          pivot_longer(cols = matches("1|2|3|4|5|6|7|8|9"), names_to = "wavelength", values_to = "value") |> 
          mutate(wavelength = as.numeric(wavelength)) |> 
          # pivot_wider(names_from = sensor, values_from = value) |> 
          # filter(wavelength >= 380, wavelength <= 700) |> 
          na.omit()
        
        # Widen for use with stats function
        df_sensor_wide <- df_sensor_long |> 
          dplyr::select(-c(dateTime, longitude, latitude)) |>
          pivot_wider(names_from = sensor, values_from = value) |> 
          na.omit()

        # get vectors
        x_vec <- df_sensor_wide[[sensors[i]]]
        y_vec <- df_sensor_wide[[sensors[j]]]
        
        # Base stats
        df_stats <- base_stats(x_vec, y_vec)
        
        # Create data.frame of results and add them to df_results
        df_res <- df_stats |> 
          mutate(sensor_X = sensors[i],
                 sensor_Y = sensors[j],
                 lon_X = df_sensor_sub$longitude[[i]],
                 lat_X = df_sensor_sub$latitude[[i]],
                 lon_Y = df_sensor_sub$longitude[[j]],
                 lat_Y = df_sensor_sub$latitude[[j]],
                 dist = hav_dist,
                 dateTime_X = df_sensor_sub$dateTime[[i]],
                 dateTime_Y = df_sensor_sub$dateTime[[j]],
                 diff_time = time_diff, .before = "n")
        df_results <- rbind(df_results, df_res)
      }
    }
  }
  
  # For loop that cycles through the requested wavelengths and calculates stats
  df_results <- df_results |> mutate(file_name = basename(file_path), .before = sensor_X)
  return(df_results)
}

# Wrapper to be able to multicore across global stats
process_global_wavelength <- function(matchup_filt, site_name, sensor_X, sensor_Y){

  # Filter data
  # matchup_filt <- filter(match_base, wavelength == wavelength_nm)
  n_match <- nrow(matchup_filt)
  
  # Calculate stats
  # if(n_match > 0){
    
  # Correct sensor labels as necessary
  if(sensor_X == "HYPERNETS"){
    sensor_X_col <- "Hyp"
  } else {
    sensor_X_col <- sensor_X
  }
  if(sensor_Y == "HYPERNETS"){
    sensor_Y_col <- "Hyp"
  } else {
    sensor_Y_col <- sensor_Y
  }
  
  # Create vectors from filtered columns
  x_vec <- matchup_filt[[sensor_X_col]]
  y_vec <- matchup_filt[[sensor_Y_col]]
  
  # Calculate statistics
  df_stats_XY <- base_stats(x_vec, y_vec)
  # df_stats_YX <- base_stats(y_vec, x_vec)
  
  # Create named objects that differ from columna names to avoid naming bug
  sensor_X_name <- sensor_X
  sensor_Y_name <- sensor_Y

  # Create data.frame of results and add them to df_results
  df_XY <- df_stats_XY |> 
    mutate(site_name = site_name, 
           var_name = "RHOW",
           sensor_X = sensor_X_name,
           sensor_Y = sensor_Y_name,
          #  Wavelength_nm = wavelength_nm,
           n_w_nm = n_match, .before = "n")
  # df_YX <- df_stats_YX |> 
  #   mutate(site_name = site_name, 
  #          var_name = "RHOW",
  #          sensor_X = sensor_Y_name,
  #          sensor_Y = sensor_X_name,
  #         #  Wavelength_nm = wavelength_nm,
  #          n_w_nm = n_match, .before = "n")
  df_both <- df_XY |> #rbind(df_XY, df_YX) |> 
    dplyr::rename(n_w_nm_clean = n)
  return(df_both)
  # } else {
  #   print(paste0("No data for wavelength ", wavelength_nm))
  # }
}

# Global stats per matchup wavelength
# site_name = "MAFR"; sensor_Y = "S3A"
# site_name = "MAFR"; sensor_Y = "S3_all"
# site_name = "THFR"; sensor_Y = "PACE"
# site_name = "MAFR"; sensor_Y = "SNPP"
# site_name = "MAFR"; sensor_Y = "JPSS1"
# site_name = "MAFR"; sensor_Y = "AQUA"
global_stats <- function(site_name, sensor_Y, daily_average = TRUE){
  
  # Create multiple folder paths if requested
  if(sensor_Y == "S3_all"){
    folder_path <- c(file_path_build(site_name, "S3A"),
                     file_path_build(site_name, "S3B"))
  } else {
    folder_path <- file_path_build(site_name, sensor_Y)
  }
  
  # Continue with satellite versions if necessary
  if(sensor_Y  == "AQUA"){
    sensor_Z <- "MODIS"
  } else if(sensor_Y == "PACE"){
    sensor_Z <- "OCI"
  } else if(sensor_Y %in% c("SNPP", "JPSS1", "JPSS2", "VIIRS_all")){
    sensor_Z <- "VIIRS"
  } else if(sensor_Y %in% c("S3A", "S3B", "S3_all")){
    sensor_Z <- "OLCI"
  } else {
  }
  
  # Get filestub based on sensor_Y
  filestub <- paste0("_",sensor_Z,".csv")
  
  # Correct sensor_X for filtering
  sensor_X <- "HYPERNETS"
  sensor_X_filt <- "Hyp"
  
  # List all files in directory
  file_list <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)
  
  # Load individual matchup results to filter file list and for further use, restricted to this
  # site -- raw matchup filenames are not unique *across* sites (see path_site_name()), so an
  # unrestricted basename match below could wrongly include/exclude a file based on another
  # site's namesake rather than this one's own QC/outlier status
  match_base_details <- read_csv(paste0("output/matchup_stats_RHOW",filestub), show_col_types = FALSE) |>
    filter(.data$site_name == .env$site_name) |>
    dplyr::select(file_name) |> distinct()
  if(nrow(match_base_details) == 0) stop("Individual matchup file not loaded correctly.")

  # Filter accordingly
  # NB: This creates the list of valid matchups after screening for spatiotemporal range
  file_list_clean <- file_list[basename(file_list) %in% match_base_details$file_name]

  # Get outlier lists, restricted to this site (see NB above)
  outliers_sat <- read_csv("meta/satellite_outliers.csv", show_col_types = FALSE) |>
    filter(.data$site_name == .env$site_name)

  # Remove outlier files
  # NB: This creates the list of valid matchups after screening for outliers in the single matchup QC process
  file_list_no_out <- file_list_clean[!basename(file_list_clean) %in% outliers_sat$file_name]
  # file_list_no_out <- file_list_clean # Or rather do not remove outliers as this is too subjective of a process
  if(length(file_list_no_out) == 0) stop(paste("No files passed QC for", site_name, sensor_X, sensor_Y))
  
  # Load data
  match_base <- furrr::future_map_dfr(file_list_no_out, load_matchup_long, .options = furrr_options(seed = TRUE))
  
  # Melt if S3_all
  if(sensor_Y == "S3_all"){
    match_base <- match_base |> 
      pivot_longer(S3A:S3B, names_to = "name", values_to = "S3_all") |> 
      dplyr::select(-name) |> 
      filter(!is.na(S3_all))
  }
  
  # Get pre-determined wavelengths
  W_nm <- W_nm_out(sensor_Y)
  
  # Filter data.frame accordingly
  match_base_filt <- filter(match_base, wavelength %in% W_nm) #|>
    # mutate(wavelength_idx = wavelength, .before = wavelength)

  # Optional day-level temporal averaging using the site-specific time window
  if(daily_average){
    match_base_filt <- daily_average_matchups(match_base_filt, site_name)
  }

  # Get the requested wavelengths global stats, add matchup count, and exit
  df_results <- match_base_filt |>
    group_by(wavelength) |>
    group_modify(~process_global_wavelength(.x, site_name = site_name, sensor_X = sensor_X, sensor_Y = sensor_Y)) |>
    ungroup() |>
    mutate(n_clean = length(file_list_clean),
           n_no_out = length(file_list_no_out),
           .before = "n_w_nm") |>
    dplyr::select(site_name, sensor_X, sensor_Y, wavelength, everything())
  return(df_results)
}

# Function that runs this over all matchup files in a directory
# site_name = "MAFR"; sensor_Y = "S3A"
# site_name = "MAFR"; sensor_Y = "SNPP"
# site_name = "MAFR"; sensor_Y = "JPSS1"
# site_name = "MAFR"; sensor_Y = "PACE"
# site_name = "MAFR"; sensor_Y = "AQUA"
process_matchup_folder <- function(site_name, sensor_Y){
  
  # Create file path
  folder_path <- file_path_build(site_name, sensor_Y)
  
  # List all files in directory
  file_list <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)
  
  # Remove files with 'all' or 'global' in the name
  file_list <- file_list[!grepl("all|global", file_list)]
  
  # Initialise results data.frame
  df_results <- furrr::future_map_dfr(file_list, process_matchup_file, .options = furrr_options(seed = TRUE)) |>
    mutate(site_name = site_name, .after = file_name)

  # Exit
  return(df_results)
}

# Process multiple folders based on request
# sensor_Z = "MODIS"; stat_choice = "matchup"
# sensor_Z = "OLCI"; stat_choice = "global"
process_sensor <- function(sensor_Z, stat_choice = "matchup", daily_average = TRUE){
  
  # Create ply grid
  ply_grid <- sensor_grid(sensor_Z)
  
  # Add S3_all if needed
  if(sensor_Z == "OLCI" & stat_choice == "global"){
    ply_grid_bonus <- data.frame(site_name = unique(ply_grid$site_name),
                                 sensor_Y = "S3_all")
    ply_grid <- rbind(ply_grid, ply_grid_bonus)
    message("Added S3_all to sensor_Y list")
  }
  
  # Process matchups and save output
  if(stat_choice == "matchup"){
    proc_res <- furrr::future_pmap_dfr(ply_grid, process_matchup_folder, .options = furrr_options(seed = TRUE))
    # Set time and distance limits
    proc_res <- proc_res |>
      mutate(diff_time_limit = site_diff_time_limit(site_name), .after = diff_time) |>
      mutate(dist_limit = 10, .after = dist)
    # Enforce time and distance constraint
    proc_res_clean <- proc_res |>
      filter(diff_time <= diff_time_limit) |>
      filter(dist <= dist_limit)
    write_csv(proc_res_clean, paste0("output/matchup_stats_RHOW_",sensor_Z,".csv"))
    # Save the matchups removed this way
    proc_res_unclean <- proc_res[!proc_res$file_name %in% proc_res_clean$file_name,]
    write_csv(proc_res_unclean, paste0("output/matchup_noQC_stats_RHOW_",sensor_Z,".csv"))
  } else {
    proc_res <- furrr::future_pmap_dfr(ply_grid, global_stats, daily_average = daily_average, .options = furrr_options(seed = TRUE))
    write_csv(proc_res, paste0("output/global_stats_RHOW_",sensor_Z,".csv"))
  }
}


# Database-sourced matchups (SQLite) ---------------------------------------
#
# Extracts HYPERNETS-vs-satellite RHOW matchup pairs directly from a
# Hypernets_matchups SQLite database (e.g. thfr_2025.db), as an alternative
# data source to the exported .csv files read via file_path_build() /
# load_matchup_mean() above.
#
# Schema (confirmed against thfr_2025.db, 2026-08-24):
#   matchups          -- one row per HYPERNETS scan. Columns HYPERNETS/S3A/S3B/
#                        PACE/AQUA/JPSS1/JPSS2/SNPP hold a measure_info.id for
#                        that sensor's matched observation, or NULL if no
#                        matchup exists for that sensor at that time.
#   measure_info      -- one row per sensor observation: data_id (->
#                        measure_data_rhow.id, start of a contiguous block),
#                        data_count (length of that block), radiometer_id (->
#                        radiometer.id), day ("YYYYMMDD"), time ("HHMMSS"),
#                        latitude, longitude, qc.
#   measure_data_rhow -- one row per replicate/pixel spectrum: id, pixel_pos
#                        ("" for HYPERNETS scan replicates; "col,row" e.g.
#                        "0,0" for satellite pixel-box positions, "0,0" =
#                        centre pixel), pixel_lat/pixel_lon (added 2026-08-26
#                        -- true per-pixel geolocation; constant/repeated
#                        across HYPERNETS replicate rows for one scan, but
#                        varies per pixel_pos within a satellite pixel box --
#                        distinct from measure_info.latitude/longitude below,
#                        which is a single per-observation reference position,
#                        not per-pixel), then one column per integer
#                        wavelength (nm) 350-1100 holding RHOW (NULL where
#                        that sensor doesn't measure that band -- satellite
#                        rows are only populated at that sensor's actual band
#                        centres, confirmed to match W_nm_out() above).
#   radiometer        -- id -> name/full_name/radiometer_type lookup.
#   rsr               -- relative spectral response function: radiometer_id,
#                        band, wl (integer nm, same grid as measure_data_rhow),
#                        rsr (raw response), rsr_p (response pre-normalised so
#                        sum(rsr_p) == 1 per band). Present for S3A/S3B/AQUA/
#                        JPSS1/JPSS2/SNPP; NOT present for PACE (OCI is itself
#                        near-hyperspectral, so no discrete-band RSR is stored).
#
# NB: unlike the exported .csv files, the db stores raw per-pixel/per-scan
# values only -- there is no precomputed "weighted mean" row. Two distinct
# aggregation questions are handled separately here:
#  (1) Satellite side: how to combine the multiple pixels in a matchup's
#      pixel box. agg_method = "mean" (default, arithmetic mean of all stored
#      pixels) or "center" (just pixel_pos "0,0").
#  (2) HYPERNETS side: how to turn its ~1 nm hyperspectral scan into a value
#      comparable to one satellite band. This is NOT a simple same-wavelength
#      lookup -- it is reconstructed here as a proper spectral convolution
#      against the satellite's own RSR (stored in `rsr`), via
#      srf_convolve_hyp(). This is almost certainly what "weighted" means in
#      the exported .csv pipeline's data_type column (see load_matchup_mean()
#      above). PACE has no stored RSR, so its HYPERNETS side falls back to a
#      same-wavelength lookup (Hyp_method == "nearest_nm" in the output, vs
#      "srf" for the other six satellites).

# Satellite column names as they appear in `matchups` and `radiometer.name`
db_satellite_names <- c("S3A", "S3B", "PACE", "AQUA", "JPSS1", "JPSS2", "SNPP")

# Infer site_name from the db file name (e.g. "thfr_2025.db" -> "THFR"),
# mirroring the MAFR/THFR site_name convention used throughout this file
db_site_name <- function(db_path){
  toupper(str_split(basename(db_path), "_")[[1]][1])
}

# Map a radiometer's stored RSR bands onto this project's fixed nominal
# wavelength labels (W_nm_out()) by nearest distance. The RSR-weighted
# centroid of a physical band and its conventional nominal label differ by a
# few nm (e.g. OLCI's nominal 665 nm band centroids at ~665.4 nm; its nominal
# 1020 nm band centroids at ~1016 nm), so exact-equality matching would
# silently drop every band. Bands with no nominal label within tol_nm are
# dropped -- this is what naturally excludes physical channels the project
# doesn't use (e.g. OLCI's O2-absorption bands at ~762/765/768 nm, or SWIR).
srf_band_map <- function(con, radiometer_id, sensor_Y, tol_nm = 5){
  nominal_nm <- W_nm_out(sensor_Y)
  band_centres <- dbGetQuery(con, sprintf("SELECT band, wl, rsr_p FROM rsr WHERE radiometer_id = %d", radiometer_id)) |>
    summarise(centre_nm = sum(wl * rsr_p), .by = band)
  if(nrow(band_centres) == 0) return(tibble(band = integer(), wavelength = numeric()))
  band_centres |>
    rowwise() |>
    mutate(wavelength = nominal_nm[which.min(abs(nominal_nm - centre_nm))],
           dist_nm = min(abs(nominal_nm - centre_nm))) |>
    ungroup() |>
    filter(dist_nm <= tol_nm) |>
    dplyr::select(band, wavelength)
}

# Convolve a full-resolution HYPERNETS spectrum against sensor_Y's relative
# spectral response function (from the db's `rsr` table) to produce a
# satellite-band-equivalent HYPERNETS spectrum. rsr_p is already normalised
# so sum(rsr_p) == 1 per band, so the band-equivalent value is a direct
# weighted sum (no extra division needed).
# hyp_full: data.frame(matchup_id, wavelength, value) -- HYPERNETS's full
# 350-1100 nm spectrum (per-scan mean) for one or more matchups
srf_convolve_hyp <- function(con, radiometer_id, sensor_Y, hyp_full, tol_nm = 5){
  band_map <- srf_band_map(con, radiometer_id, sensor_Y, tol_nm = tol_nm) |>
    dplyr::rename(wavelength_nominal = wavelength)
  if(nrow(band_map) == 0) return(tibble())

  rsr_tbl <- dbGetQuery(con, sprintf("SELECT band, wl, rsr_p FROM rsr WHERE radiometer_id = %d", radiometer_id)) |>
    inner_join(band_map, by = "band")

  hyp_full |>
    inner_join(rsr_tbl, by = c("wavelength" = "wl"), relationship = "many-to-many") |>
    summarise(value = sum(value * rsr_p, na.rm = TRUE),
              n_wl = n(),
              .by = c(matchup_id, wavelength_nominal)) |>
    dplyr::rename(wavelength = wavelength_nominal)
}

# Aggregate the raw measure_data_rhow rows belonging to one or more
# measure_info entries into one spectrum (long format) per measure_info.id.
# info_tbl: data.frame with columns id, data_id, data_count (from measure_info)
db_load_spectra <- function(con, info_tbl, agg_method = c("mean", "center")){
  agg_method <- match.arg(agg_method)

  # Expand each measure_info's contiguous [data_id, data_id + data_count - 1]
  # block into individual measure_data_rhow ids, tagged with their parent
  info_map <- info_tbl |>
    mutate(data_id_end = data_id + data_count - 1) |>
    rowwise() |>
    reframe(measure_info_id = id, data_row_id = seq(data_id, data_id_end)) |>
    ungroup()

  # Pull the raw spectra for exactly the rows needed
  spec_query <- paste0("SELECT * FROM measure_data_rhow WHERE id IN (",
                       paste(unique(info_map$data_row_id), collapse = ","), ")")
  spec_raw <- dbGetQuery(con, spec_query)

  spec_long <- spec_raw |>
    pivot_longer(cols = -c(id, pixel_pos, pixel_lat, pixel_lon), names_to = "wavelength", values_to = "value") |>
    filter(!is.na(value)) |>
    mutate(wavelength = as.numeric(wavelength)) |>
    left_join(info_map, by = c("id" = "data_row_id"))

  # Aggregate replicates/pixels per measure_info entry x wavelength
  # NB: pixel_pos == "" identifies HYPERNETS scan replicates (always averaged);
  # "center" agg_method only applies to satellite pixel boxes (pixel_pos == "0,0")
  if(agg_method == "center"){
    spec_long <- spec_long |>
      filter(pixel_pos %in% c("", "0,0"))
  }

  # pixel_lat_mean/pixel_lon_mean: mean per-pixel geolocation (pixel_lat/
  # pixel_lon from measure_data_rhow) of whichever rows contributed to
  # value_mean -- under agg_method = "center" this is exactly the centre
  # pixel's true coordinate; under "mean" it is the pixel-box centroid
  spec_agg <- spec_long |>
    summarise(value_mean = mean(value, na.rm = TRUE),
              value_sd = sd(value, na.rm = TRUE),
              n_used = n(),
              pixel_lat_mean = mean(pixel_lat, na.rm = TRUE),
              pixel_lon_mean = mean(pixel_lon, na.rm = TRUE),
              .by = c(measure_info_id, wavelength)) |>
    mutate(cv_pct = 100 * abs(value_sd / value_mean))

  return(spec_agg)
}

# Like db_load_spectra() above but skips the mean/center aggregation step --
# returns one row per (measure_info_id, pixel_pos, wavelength), keeping every
# individual pixel's own value/pixel_lat/pixel_lon. Used by db_matchup_pixels()
# below for per-pixel spatial diagnostics (see meta/pixel_explore.R).
db_load_spectra_pixel <- function(con, info_tbl){
  info_map <- info_tbl |>
    mutate(data_id_end = data_id + data_count - 1) |>
    rowwise() |>
    reframe(measure_info_id = id, data_row_id = seq(data_id, data_id_end)) |>
    ungroup()

  spec_query <- paste0("SELECT * FROM measure_data_rhow WHERE id IN (",
                       paste(unique(info_map$data_row_id), collapse = ","), ")")
  spec_raw <- dbGetQuery(con, spec_query)

  spec_long <- spec_raw |>
    pivot_longer(cols = -c(id, pixel_pos, pixel_lat, pixel_lon), names_to = "wavelength", values_to = "value") |>
    filter(!is.na(value)) |>
    mutate(wavelength = as.numeric(wavelength)) |>
    left_join(info_map, by = c("id" = "data_row_id")) |>
    dplyr::select(measure_info_id, pixel_pos, pixel_lat, pixel_lon, wavelength, value)

  return(spec_long)
}

# Extract HYPERNETS-vs-sensor_Y matchup pairs from a matchup db, in long
# format: one row per matchup x wavelength, with a "Hyp" column and a
# sensor_Y-named column ready for use with base_stats()/process_global_wavelength()
# (which already expect a "Hyp" column name for HYPERNETS).
# db_path <- "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db"; sensor_Y <- "S3A"
db_matchup_long <- function(db_path, sensor_Y, agg_method = c("mean", "center")){
  agg_method <- match.arg(agg_method)
  sensor_Y <- match.arg(sensor_Y, db_satellite_names)
  db_path <- path.expand(db_path)
  site_name <- db_site_name(db_path)

  con <- dbConnect(RSQLite::SQLite(), db_path)
  on.exit(dbDisconnect(con), add = TRUE)

  # Matched (HYPERNETS, sensor_Y) measure_info.id pairs
  match_query <- sprintf(
    "SELECT id AS matchup_id, HYPERNETS AS hyp_info_id, %s AS sat_info_id FROM matchups WHERE HYPERNETS IS NOT NULL AND %s IS NOT NULL",
    sensor_Y, sensor_Y)
  match_ids <- dbGetQuery(con, match_query)
  if(nrow(match_ids) == 0){
    warning(paste0("No matchups found in ", basename(db_path), " for sensor_Y = ", sensor_Y))
    return(tibble())
  }

  # Metadata (day/time/lat/lon) + data block pointers for every measure_info
  # row referenced by these matchups
  info_ids <- unique(c(match_ids$hyp_info_id, match_ids$sat_info_id))
  info_query <- paste0("SELECT id, data_id, data_count, day, time, latitude, longitude, qc FROM measure_info WHERE id IN (",
                       paste(info_ids, collapse = ","), ")")
  info_tbl <- dbGetQuery(con, info_query)

  # Aggregate spectra per measure_info entry
  spec_agg <- db_load_spectra(con, info_tbl, agg_method = agg_method)

  # Attach metadata to the aggregated spectra
  info_meta <- info_tbl |>
    mutate(dateTime = as.POSIXct(paste(day, time), format = "%Y%m%d %H%M%S", tz = "Europe/Paris")) |>
    dplyr::select(id, dateTime, latitude, longitude, qc)

  spec_meta <- spec_agg |>
    left_join(info_meta, by = c("measure_info_id" = "id"))

  # HYPERNETS side: full-resolution per-scan-mean spectrum, tagged by matchup_id
  # lon_Hyp/lat_Hyp come from measure_info (single reference position per scan);
  # pixel_lon_Hyp/pixel_lat_Hyp come from measure_data_rhow.pixel_lat/pixel_lon
  # (per-pixel geolocation, added 2026-08-26 -- for HYPERNETS rows this is
  # constant across replicates, i.e. the same fixed station coordinate)
  hyp_full <- match_ids |>
    dplyr::select(matchup_id, measure_info_id = hyp_info_id) |>
    left_join(spec_meta, by = "measure_info_id", relationship = "many-to-many") |>
    dplyr::select(matchup_id, wavelength, value = value_mean, dateTime,
                  lon_Hyp = longitude, lat_Hyp = latitude,
                  pixel_lon_Hyp = pixel_lon_mean, pixel_lat_Hyp = pixel_lat_mean)

  sat_radiometer_id <- dbGetQuery(con, sprintf("SELECT id FROM radiometer WHERE name = '%s'", sensor_Y))$id
  has_rsr <- dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM rsr WHERE radiometer_id = %d", sat_radiometer_id))$n > 0

  if(has_rsr){
    # Reconstruct the satellite's band-equivalent HYPERNETS value by
    # convolving the full HYPERNETS spectrum against the satellite's stored
    # RSR -- the scientifically correct way to compare a hyperspectral
    # in-situ instrument against a discrete-band satellite sensor
    hyp_meta <- hyp_full |> dplyr::select(matchup_id, dateTime, lon_Hyp, lat_Hyp, pixel_lon_Hyp, pixel_lat_Hyp) |> distinct()
    hyp_side <- srf_convolve_hyp(con, sat_radiometer_id, sensor_Y, hyp_full) |>
      left_join(hyp_meta, by = "matchup_id") |>
      dplyr::rename(Hyp = value, Hyp_n = n_wl) |>
      mutate(Hyp_method = "srf") |>
      dplyr::select(matchup_id, wavelength, Hyp, Hyp_n, Hyp_method, dateTime_Hyp = dateTime, lon_Hyp, lat_Hyp, pixel_lon_Hyp, pixel_lat_Hyp)
  } else {
    # No stored RSR for this sensor (PACE/OCI: itself near-hyperspectral) --
    # fall back to a same-wavelength lookup
    hyp_side <- hyp_full |>
      dplyr::select(matchup_id, wavelength, Hyp = value, dateTime_Hyp = dateTime, lon_Hyp, lat_Hyp, pixel_lon_Hyp, pixel_lat_Hyp) |>
      mutate(Hyp_n = NA_integer_, Hyp_method = "nearest_nm")
  }

  # NB: relationship = "many-to-many" is expected, not a bug -- a single
  # satellite overpass (one measure_info_id) commonly matches several
  # different HYPERNETS scans, and vice versa
  # lon_/lat_<sensor_Y> come from measure_info (single reference position);
  # pixel_lon_/pixel_lat_<sensor_Y> come from measure_data_rhow.pixel_lat/
  # pixel_lon (per-pixel geolocation) -- under agg_method = "center" this is
  # the true centre-pixel coordinate, under "mean" the pixel-box centroid
  sat_side <- match_ids |>
    dplyr::select(matchup_id, measure_info_id = sat_info_id) |>
    left_join(spec_meta, by = "measure_info_id", relationship = "many-to-many") |>
    dplyr::select(matchup_id, wavelength,
                  !!sensor_Y := value_mean,
                  !!paste0(sensor_Y, "_n") := n_used,
                  !!paste0(sensor_Y, "_cv_pct") := cv_pct,
                  !!paste0("dateTime_", sensor_Y) := dateTime,
                  !!paste0("lon_", sensor_Y) := longitude,
                  !!paste0("lat_", sensor_Y) := latitude,
                  !!paste0("pixel_lon_", sensor_Y) := pixel_lon_mean,
                  !!paste0("pixel_lat_", sensor_Y) := pixel_lat_mean)

  # NB: filter on the value columns specifically, not na.omit() -- Hyp_n/cv_pct
  # can be legitimately NA (e.g. n_used == 1 under agg_method = "center") and
  # that must not drop an otherwise-valid pair
  df_pairs <- inner_join(hyp_side, sat_side, by = c("matchup_id", "wavelength")) |>
    filter(!is.na(Hyp), !is.na(.data[[sensor_Y]]))

  if(nrow(df_pairs) == 0){
    warning(paste0("Matchups found but no overlapping wavelengths for sensor_Y = ", sensor_Y))
    return(tibble())
  }

  # Per-matchup distance (km) and time difference (min), for the same QC role
  # dist/diff_time play in process_matchup_file()/process_sensor() above --
  # filtering against dist_limit / site_diff_time_limit(site_name) is left to
  # the caller
  df_pairs <- df_pairs |>
    mutate(dist_km = round(distHaversine(cbind(lon_Hyp, lat_Hyp), cbind(!!sym(paste0("lon_", sensor_Y)), !!sym(paste0("lat_", sensor_Y)))) / 1000, 3),
           diff_time_min = round(abs(as.numeric(difftime(dateTime_Hyp, !!sym(paste0("dateTime_", sensor_Y)), units = "mins")))),
           match_date = as.Date(dateTime_Hyp),
           site_name = site_name,
           sensor_Y = sensor_Y,
           .before = 1)

  return(df_pairs)
}

# Like db_matchup_long() above but keeps the satellite side at raw per-pixel
# resolution (one row per matchup x wavelength x pixel_pos) instead of
# aggregating the pixel box -- used for the spatial/directional diagnostics
# in meta/pixel_explore.R. HYPERNETS side is still averaged across scan
# replicates (pixel_pos == "" rows), matching db_matchup_long(agg_method = "mean").
# db_path <- "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db"; sensor_Y <- "S3A"
db_matchup_pixels <- function(db_path, sensor_Y){
  sensor_Y <- match.arg(sensor_Y, db_satellite_names)
  db_path <- path.expand(db_path)
  site_name <- db_site_name(db_path)

  con <- dbConnect(RSQLite::SQLite(), db_path)
  on.exit(dbDisconnect(con), add = TRUE)

  match_query <- sprintf(
    "SELECT id AS matchup_id, HYPERNETS AS hyp_info_id, %s AS sat_info_id FROM matchups WHERE HYPERNETS IS NOT NULL AND %s IS NOT NULL",
    sensor_Y, sensor_Y)
  match_ids <- dbGetQuery(con, match_query)
  if(nrow(match_ids) == 0){
    warning(paste0("No matchups found in ", basename(db_path), " for sensor_Y = ", sensor_Y))
    return(tibble())
  }

  info_ids <- unique(c(match_ids$hyp_info_id, match_ids$sat_info_id))
  info_query <- paste0("SELECT id, data_id, data_count, day, time, latitude, longitude, qc FROM measure_info WHERE id IN (",
                       paste(info_ids, collapse = ","), ")")
  info_tbl <- dbGetQuery(con, info_query)

  spec_pixel <- db_load_spectra_pixel(con, info_tbl)

  info_meta <- info_tbl |>
    mutate(dateTime = as.POSIXct(paste(day, time), format = "%Y%m%d %H%M%S", tz = "Europe/Paris")) |>
    dplyr::select(id, dateTime, latitude, longitude, qc)

  # HYPERNETS side: average across scan replicates (pixel_pos == ""), same as
  # db_load_spectra(agg_method = "mean") -- one value per matchup x wavelength
  hyp_agg <- spec_pixel |>
    filter(pixel_pos == "") |>
    summarise(value = mean(value, na.rm = TRUE),
              pixel_lat_Hyp = mean(pixel_lat, na.rm = TRUE),
              pixel_lon_Hyp = mean(pixel_lon, na.rm = TRUE),
              .by = c(measure_info_id, wavelength)) |>
    left_join(info_meta, by = c("measure_info_id" = "id"))

  hyp_full <- match_ids |>
    dplyr::select(matchup_id, measure_info_id = hyp_info_id) |>
    left_join(hyp_agg, by = "measure_info_id", relationship = "many-to-many") |>
    dplyr::select(matchup_id, wavelength, value, dateTime,
                  lon_Hyp = longitude, lat_Hyp = latitude,
                  pixel_lon_Hyp, pixel_lat_Hyp)

  sat_radiometer_id <- dbGetQuery(con, sprintf("SELECT id FROM radiometer WHERE name = '%s'", sensor_Y))$id
  has_rsr <- dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM rsr WHERE radiometer_id = %d", sat_radiometer_id))$n > 0

  if(has_rsr){
    hyp_meta <- hyp_full |> dplyr::select(matchup_id, dateTime, lon_Hyp, lat_Hyp, pixel_lon_Hyp, pixel_lat_Hyp) |> distinct()
    hyp_side <- srf_convolve_hyp(con, sat_radiometer_id, sensor_Y, hyp_full) |>
      left_join(hyp_meta, by = "matchup_id") |>
      dplyr::rename(Hyp = value, Hyp_n = n_wl) |>
      mutate(Hyp_method = "srf") |>
      dplyr::select(matchup_id, wavelength, Hyp, Hyp_n, Hyp_method, dateTime_Hyp = dateTime, lon_Hyp, lat_Hyp, pixel_lon_Hyp, pixel_lat_Hyp)
  } else {
    hyp_side <- hyp_full |>
      dplyr::select(matchup_id, wavelength, Hyp = value, dateTime_Hyp = dateTime, lon_Hyp, lat_Hyp, pixel_lon_Hyp, pixel_lat_Hyp) |>
      mutate(Hyp_n = NA_integer_, Hyp_method = "nearest_nm")
  }

  # Satellite side: raw per-pixel rows, no aggregation -- one row per
  # matchup x wavelength x pixel_pos
  sat_side <- match_ids |>
    dplyr::select(matchup_id, measure_info_id = sat_info_id) |>
    left_join(spec_pixel, by = "measure_info_id", relationship = "many-to-many") |>
    left_join(info_meta, by = c("measure_info_id" = "id")) |>
    dplyr::select(matchup_id, wavelength, pixel_pos, pixel_lat, pixel_lon,
                  !!sensor_Y := value,
                  !!paste0("dateTime_", sensor_Y) := dateTime,
                  !!paste0("lon_", sensor_Y) := longitude,
                  !!paste0("lat_", sensor_Y) := latitude)

  df_pairs <- inner_join(hyp_side, sat_side, by = c("matchup_id", "wavelength")) |>
    filter(!is.na(Hyp), !is.na(.data[[sensor_Y]]))

  if(nrow(df_pairs) == 0){
    warning(paste0("Matchups found but no overlapping wavelengths for sensor_Y = ", sensor_Y))
    return(tibble())
  }

  # Per-pixel distance (km) from the HYPERNETS station to this specific
  # satellite pixel, used downstream for bearing/quadrant diagnostics
  df_pairs <- df_pairs |>
    mutate(dist_km = round(distHaversine(cbind(pixel_lon, pixel_lat), cbind(lon_Hyp, lat_Hyp)) / 1000, 3),
           diff_time_min = round(abs(as.numeric(difftime(dateTime_Hyp, !!sym(paste0("dateTime_", sensor_Y)), units = "mins")))),
           match_date = as.Date(dateTime_Hyp),
           site_name = site_name,
           sensor_Y = sensor_Y,
           .before = 1)

  return(df_pairs)
}

# Convenience wrapper: run db_matchup_long() for every satellite present in
# the matchups table (or a subset via sensor_Y_vec), row-binding the results.
# NB: this row-binds different sensors' pair columns (S3A, AQUA, ...) using
# bind_rows(), so each sensor's satellite-value column stays separate rather
# than being merged into one column -- filter/split by sensor_Y before
# further analysis.
db_matchup_all <- function(db_path, sensor_Y_vec = db_satellite_names, agg_method = c("mean", "center")){
  agg_method <- match.arg(agg_method)
  map_dfr(sensor_Y_vec, function(sY){
    tryCatch(db_matchup_long(db_path, sY, agg_method = agg_method),
             warning = function(w){ message(conditionMessage(w)); tibble() })
  })
}


# TEMPORARY -- NE-quadrant-filtered THFR matchup CSVs -----------------------
# meta/pixel_explore.R / meta/pixel_explore_output/summary.md found that THFR's satellite
# pixel-box extraction is centered ~2.6 km SW of the true (fixed) HYPERNETS station position, and
# that there is a real, spectrally coherent quadrant-dependent bias in the resulting pixel field.
# The two functions below regenerate THFR's raw per-matchup RHOW CSV files with the satellite
# side restricted to just the NE quadrant (pixel_lat >= lat_Hyp & pixel_lon >= lon_Hyp, same
# hard-cut definition as pixel_bearing_quadrant() in meta/pixel_explore.R) before re-averaging,
# preserving the exact on-disk schema so load_matchup_mean()/load_matchup_long()/
# sat_var_check()/process_matchup_folder() can read the output completely unmodified. Written to
# a new parallel site folder ("THFR_NE") fully isolated from the real THFR/MAFR data delivered by
# Hypernets_matchups. Nothing in the existing pipeline is affected unless/until "THFR_NE" is
# separately added to available_sites()'s candidate list.

# Helper for db_export_matchups_ne() below. Writes one NE-quadrant-filtered matchup as a raw
# RHOW CSV, by locating the matching original THFR file(s) and copying everything except the
# satellite sensor's rows verbatim. The HYPERNETS rows are untouched because the SRF-convolved
# HYPERNETS value depends only on the HYPERNETS scan + the satellite's RSR, not on which
# satellite pixels get averaged.
# df_ne_matchup: one matchup_id's NE-quadrant satellite pixel rows for one sensor_Y (a subset of
# db_matchup_pixels()'s output, already filtered to pixel_lat >= lat_Hyp & pixel_lon >= lon_Hyp)
write_matchup_csv_ne <- function(df_ne_matchup, sensor_Y, out_dir, tag = "NE"){
  sat_dt_col <- paste0("dateTime_", sensor_Y)

  hyp_ts <- format(df_ne_matchup$dateTime_Hyp[1], "%Y%m%dT%H%M%S", tz = "Europe/Paris")
  sat_ts <- format(df_ne_matchup[[sat_dt_col]][1], "%Y%m%dT%H%M%S", tz = "Europe/Paris")

  # Locate the original THFR file(s) for this matchup (both _R and _C variants, if both exist)
  orig_dir <- file_path_build("THFR", sensor_Y)
  orig_prefix <- paste0(sensor_Y, "_", sat_ts, "_vs_HYPERNETS_", hyp_ts)
  orig_files <- list.files(orig_dir, pattern = paste0("^", orig_prefix, ".*\\.csv$"), full.names = TRUE)

  if(length(orig_files) == 0){
    message("  No original THFR file found for ", orig_prefix, " -- skipping")
    return(invisible(NULL))
  }

  # Per-wavelength NE-quadrant stats for this one matchup
  wl_stats <- df_ne_matchup |>
    summarise(mean_val = mean(.data[[sensor_Y]], na.rm = TRUE),
              sd_val = sd(.data[[sensor_Y]], na.rm = TRUE),
              n_used = n(),
              .by = wavelength)
  if(nrow(wl_stats) == 0) return(invisible(NULL))

  # Single-scalar CV proxy (variability_centered is one value per file, not per wavelength) --
  # exact upstream Hypernets_matchups formula unknown. This reuses db_load_spectra()'s cv_pct
  # logic (100 * sd/|mean|), pooled across wavebands with >= 2 contributing NE pixels. Uses the
  # median rather than the mean across wavebands. A handful of near-zero-mean bands (e.g.
  # S3A's 665/673/681 nm) otherwise produce cv_pct in the thousands of percent and dominate a
  # simple mean, even though most bands sit in a plausible range (confirmed empirically).
  cv_pct_by_wl <- wl_stats |> filter(n_used >= 2, mean_val != 0) |>
    mutate(cv_pct = 100 * abs(sd_val / mean_val))
  # Fall back to 0 (not NA) when every waveband has < 2 contributing NE pixels. With a single
  # pixel there is no measurable spatial spread, so "no pixel disagreement detected" is the
  # correct reading, and it keeps variability_centered numeric (sat_var_check() errors on an
  # all-NA column, since it expects exactly one non-NA value per file)
  cv_scalar <- if(nrow(cv_pct_by_wl) > 0) median(cv_pct_by_wl$cv_pct, na.rm = TRUE) else 0

  ne_lat <- mean(df_ne_matchup$pixel_lat, na.rm = TRUE)
  ne_lon <- mean(df_ne_matchup$pixel_lon, na.rm = TRUE)
  ne_pixel_pos <- paste(sort(unique(df_ne_matchup$pixel_pos)), collapse = ";")

  for(orig_file in orig_files){
    df_orig <- suppressMessages(read_delim(orig_file, delim = ";", col_types = cols(.default = "c"), show_col_types = FALSE))
    colnames(df_orig)[1] <- "sensor"

    is_sat_row <- grepl(paste0("^", sensor_Y, " "), df_orig$sensor)
    if(!any(is_sat_row)){
      message("  No ", sensor_Y, " rows found in ", basename(orig_file), " -- skipping")
      next
    }

    wl_cols <- setdiff(colnames(df_orig),
                       c("sensor", "data_type", "day", "time", "latitude", "longitude",
                         "radiometer_id", "pixel_pos", "type", "variability_centered"))

    df_new <- df_orig
    df_new$latitude[is_sat_row] <- as.character(round(ne_lat, 6))
    df_new$longitude[is_sat_row] <- as.character(round(ne_lon, 6))
    df_new$pixel_pos[is_sat_row] <- ne_pixel_pos
    df_new$variability_centered[is_sat_row] <- if(is.na(cv_scalar)) NA_character_ else as.character(round(cv_scalar, 6))

    # Recompute each satellite row's wavelength values from the NE-quadrant subset. Role
    # (mean/max/min) is read off the row's own (possibly bug-quirky, e.g. "rhow weighted")
    # data_type label so per-sensor label variants (see manuscript/upstream-data-bugs.md, Bugs
    # 1/2/7) are preserved untouched
    for(ridx in which(is_sat_row)){
      dt <- tolower(df_new$data_type[ridx])
      role <- if(grepl("std_max", dt)) "max" else if(grepl("std_min", dt)) "min" else "mean"
      for(wl in wl_cols){
        wl_num <- suppressWarnings(as.numeric(wl))
        row_stat <- wl_stats[wl_stats$wavelength == wl_num, ]
        if(nrow(row_stat) == 0 || is.na(row_stat$mean_val)){
          df_new[[wl]][ridx] <- NA_character_
          next
        }
        val <- switch(role,
                      mean = row_stat$mean_val,
                      max  = row_stat$mean_val + row_stat$sd_val,
                      min  = row_stat$mean_val - row_stat$sd_val)
        df_new[[wl]][ridx] <- if(is.na(val)) NA_character_ else as.character(val)
      }
    }

    out_name <- sub("\\.csv$", paste0("_", tag, ".csv"), basename(orig_file))
    write_delim(df_new, file.path(out_dir, out_name), delim = ";", quote = "needed", na = "")
  }
  invisible(NULL)
}

# Pixel-inclusion rules -------------------------------------------------------
# Each takes the QC-gated per-pixel df_pixel (from db_matchup_pixels(), already filtered to
# matchups passing the time/distance gate) and returns the subset of rows to keep. Passed into
# db_export_matchups_site() below -- add a new derived site by writing one of these plus a one-
# line wrapper, nothing else needs to change.

# Restrict to the inner 3x3 grid (col/row in -1..1), then to the NE quadrant relative to the true
# station position. The 3x3 restriction matches the box size Hypernets_matchups actually exports
# to .csv (confirmed by inspecting a real THFR file: pixel_pos values only span -1..1, even
# though the .db's raw measure_data_rhow block stores the wider 5x5 grid col/row -2..2) --
# excluding any pixel_pos containing "2" drops the outer ring. A hand-traced check on one S3A
# matchup found an outer-ring pixel (pixel_pos "-1,2") carrying a large negative outlier value
# that was dragging the NE-quadrant mean well below the true NE water signal; restricting to the
# inner 3x3 grid removes exactly that kind of edge/outer-ring artifact.
pixel_filter_ne_inner3x3 <- function(df_pixel){
  df_pixel |>
    filter(!grepl("2", pixel_pos)) |>
    filter(pixel_lat >= lat_Hyp, pixel_lon >= lon_Hyp)
}

# Hand-drawn "clean water" polygon, N/NE of the THFR station -- traces open water while
# deliberately excluding land and the visible oyster-bed rows. Traced by eye off
# Esri.WorldImagery at zoom 17; see meta/pixel_explore.R's "clean_water_polygon_NE.png" for the
# traced overlay and its own comments for known gaps (a couple of small oyster-table clusters not
# yet excised). Starting point for manual refinement, not a finished mask -- single-sourced here
# so meta/pixel_explore.R (which sources this file) and the pipeline never drift apart.
clean_water_polygon <- tribble(
  ~lon,   ~lat,
  3.6660, 43.4350,
  3.6625, 43.4385,
  3.6595, 43.4415,
  3.6595, 43.4440,
  3.6630, 43.4448,
  3.6660, 43.4450,
  3.6710, 43.4442,
  3.6705, 43.4400,
  3.6690, 43.4375,
  3.6685, 43.4350,
  3.6660, 43.4350  # closes the ring
)
clean_water_sf <- st_sf(geometry = st_sfc(st_polygon(list(as.matrix(clean_water_polygon))), crs = 4326))

# Keep only pixels whose true geolocation (pixel_lon/pixel_lat) falls inside clean_water_sf --
# unlike pixel_filter_ne_inner3x3(), this is a pure real-world-geometry test, independent of the
# pixel-box grid coordinates, so no inner/outer-ring restriction is needed on top of it.
pixel_filter_clean_water <- function(df_pixel){
  pts <- sf::st_as_sf(df_pixel, coords = c("pixel_lon", "pixel_lat"), crs = 4326, remove = FALSE)
  df_pixel[lengths(sf::st_within(pts, clean_water_sf)) > 0, ]
}

# Regenerates THFR's raw per-matchup RHOW CSV files for every sensor_Y in a sensor family,
# restricted to whichever pixels pixel_filter_fn keeps, and written under a new site folder
# named site_name (isolated from the real THFR/MAFR data). See the section comment above for full
# context. Matchups with zero surviving pixels for a given sensor are skipped entirely (no file
# written) -- for the coarser sensors (AQUA/VIIRS/PACE) this is expected to remove most matchups
# under the NE+3x3 filter, since the box is centered well outside the NE quadrant for those (see
# meta/pixel_explore_output/summary.md); a similarly small survival rate is expected under the
# clean-water-polygon filter, which is deliberately small and close to the station.
# sensor_Z = "OLCI"; site_name = "THFR_NE"; pixel_filter_fn = pixel_filter_ne_inner3x3; file_tag = "NE"
db_export_matchups_site <- function(sensor_Z, site_name, pixel_filter_fn, file_tag,
                                     db_path = "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db"){
  sensor_Y_list <- sensor_grid(sensor_Z)$sensor_Y |> unique()

  for(sensor_Y in sensor_Y_list){
    message("db_export_matchups_site(", site_name, "): ", sensor_Y)
    df_pixel <- tryCatch(db_matchup_pixels(db_path, sensor_Y),
                          warning = function(w){ message(conditionMessage(w)); tibble() })
    if(nrow(df_pixel) == 0) next

    # db_matchup_pixels() returns every candidate pairing in the db's raw `matchups` table,
    # time/distance QC is deliberately left to the caller (see its own docstring). The exported
    # .csv files on disk only exist for matchups that already passed this gate, so it must be
    # applied here before looking for a matching original file. diff_time_min is already
    # matchup-level (constant per matchup_id), dist_km as returned by db_matchup_pixels() is
    # per-PIXEL (station-to-pixel, used for the bearing/quadrant diagnostics), so the
    # matchup-level distance (station-to-satellite-reference-position, matching how
    # db_matchup_long()/process_sensor() gate distance) is recomputed here instead. Both derived
    # sites are reinterpretations of the same underlying THFR matchups, so this stays hardcoded
    # to THFR's own QC policy regardless of site_name.
    sat_lon_col <- paste0("lon_", sensor_Y); sat_lat_col <- paste0("lat_", sensor_Y)
    df_pixel <- df_pixel |>
      mutate(dist_km_matchup = distHaversine(cbind(lon_Hyp, lat_Hyp), cbind(.data[[sat_lon_col]], .data[[sat_lat_col]])) / 1000) |>
      filter(diff_time_min <= site_diff_time_limit("THFR"), dist_km_matchup <= 10)
    if(nrow(df_pixel) == 0){
      message("  No matchups for ", sensor_Y, " pass the time/distance QC gate -- skipping")
      next
    }

    df_filt <- pixel_filter_fn(df_pixel)
    if(nrow(df_filt) == 0){
      message("  No pixels survive the ", site_name, " filter for ", sensor_Y, " -- skipping")
      next
    }

    out_dir <- file_path_build(site_name, sensor_Y)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    for(mid in unique(df_filt$matchup_id)){
      write_matchup_csv_ne(filter(df_filt, matchup_id == mid), sensor_Y, out_dir, tag = file_tag)
    }
  }
}

# Thin, site-specific wrappers -- add a new derived site by writing one of these (plus, if the
# inclusion rule is new, a pixel_filter_*() function above)
db_export_matchups_ne <- function(sensor_Z, db_path = "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db"){
  db_export_matchups_site(sensor_Z, "THFR_NE", pixel_filter_ne_inner3x3, "NE", db_path)
}
db_export_matchups_poly <- function(sensor_Z, db_path = "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db"){
  db_export_matchups_site(sensor_Z, "THFR_poly", pixel_filter_clean_water, "poly", db_path)
}


# Plotting functions ------------------------------------------------------

# Plot data based on wavelength group
plot_matchup_nm <- function(df, x_sensor, y_sensor){
  colours_nm <- colour_nm_func(y_sensor)
  df_prep <- df |> filter(!is.na(!!sym(x_sensor)), !is.na(!!sym(y_sensor)))
  if(y_sensor == "PACE"){
    df_prep <- df_prep |>
      mutate(wavelength = cut(wavelength,
                              breaks = c(350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1050),
                              labels = names(colours_nm),
                              include.lowest = TRUE, right = TRUE))
  }
  df_prep |>
    ggplot(aes_string(x = x_sensor, y = y_sensor)) +
    geom_point(aes(colour = as.factor(wavelength), shape = site_name), size = 3) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    labs(title = paste("RHOW","-", x_sensor, "vs", y_sensor),
         x = paste("RHOW", x_sensor),
         y = paste("RHOW", y_sensor),
         colour = "Wavelength (nm)") +
    scale_colour_manual(values = colours_nm) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "bottom")
}

# Plot based on date of collection
plot_matchup_date <- function(df, x_sensor, y_sensor){
  df |>
    filter(!is.na(!!sym(x_sensor)), !is.na(!!sym(y_sensor))) |>
    mutate(date = as.factor(as.Date(dateTime_X))) |>
    ggplot(aes_string(x = x_sensor, y = y_sensor)) +
    geom_point(aes(colour = date, shape = site_name), size = 3) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    labs(title = paste("RHOW","-", x_sensor, "vs", y_sensor),
         x = paste("RHOW", x_sensor),
         y = paste("RHOW", y_sensor),
         colour = "date") +
    # scale_colour_brewer(palette = "Dark2")  +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "bottom")
}

# Plot based on dateTime of collection
plot_matchup_dateTime <- function(df, x_sensor, y_sensor, date_filter){
  df |>
    filter(!is.na(!!sym(x_sensor)), !is.na(!!sym(y_sensor))) |>
    mutate(date = as.Date(dateTime_X)) |>
    filter(date == as.Date(date_filter)) |>
    ggplot(aes_string(x = x_sensor, y = y_sensor)) +
    geom_point(aes(colour = dateTime_X, shape = site_name), size = 3) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    labs(title = paste("RHOW","-", x_sensor, "vs", y_sensor,"-", date_filter),
         x = paste("RHOW", x_sensor),
         y = paste("RHOW", y_sensor),
         colour = "time (UTC)") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "bottom")
}

# Plot scatterplot based on the MAPE values of each comparison
plot_matchup_Error_Bias <- function(df, x_sensor, y_sensor){
  pl_Error <- df |>
    filter(!is.na(!!sym(x_sensor)), !is.na(!!sym(y_sensor))) |>
    ggplot(aes_string(x = x_sensor, y = y_sensor)) +
    geom_point(aes(colour = Error_50, shape = site_name), size = 3) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    scale_colour_viridis_c(option = "D") +
    labs(title = paste(x_sensor, "vs", y_sensor,"- Error"),
         x = paste(x_sensor),
         y = paste(y_sensor),
         colour = "Error [%]") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "bottom")
  pl_Bias <- df |>
    filter(!is.na(!!sym(x_sensor)), !is.na(!!sym(y_sensor))) |>
    ggplot(aes_string(x = x_sensor, y = y_sensor)) +
    geom_point(aes(colour = Bias_50, shape = site_name), size = 3) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    scale_colour_viridis_c(option = "A") +
    labs(title = paste("RHOW","-", x_sensor, "vs", y_sensor,"- Bias"),
         x = paste("RHOW", x_sensor),
         y = paste("RHOW", y_sensor),
         colour = "Bias [%]") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "bottom")
  ggpubr::ggarrange(pl_Error, pl_Bias, nrow = 2, ncol = 1)
}

# Code that preps the labels for plotting based on a given input
pretty_label_func <- function(char_string){
  
  # Set default values
  # NB: Many names are already correct and don't need to be tuned below
  sensor_col <- char_string
  sensor_lab <- char_string
  units_lab <- char_string
  
  # Correct satellite names
  if(char_string == "HYPERNETS"){
    sensor_col <- "Hyp"
  } else if(char_string == "AQUA"){
    sensor_lab <- "MODIS-A" 
  } else if(char_string == "S3A"){ # NB: This is a special case
    sensor_lab <- "S3A" 
  } else if(char_string == "S3B"){ # NB: This is a special case
    sensor_lab <- "S3B" 
  } else if(char_string == "S3"){ # NB: This is a special case
    sensor_lab <- "S3A; S3B" 
  } else if(char_string == "PACE_V2"){
    sensor_lab <- "PACE v2.0"
  } else if(char_string == "PACE_V30"){
    sensor_lab <- "PACE v3.0"
  } else if(char_string == "PACE_V31"){
    sensor_lab <- "PACE v3.1"
  } else if(char_string == "SNPP"){
    sensor_lab <- "SNPP"
  } else if(char_string == "JPSS1"){
    sensor_lab <- "JPSS1"
  } else if(char_string == "JPSS2"){
    sensor_lab <- "JPSS2"
  }
  
  # Correct variable units
  if(char_string == "RHOW"){
    units_lab <- "<i>ρ<sub>w</sub></i>"
  }
  
  # Combine into data.frame
  if(char_string %in% c("RHOW")){
    pretty_labels <- data.frame(var_name = char_string,
                                units_lab = units_lab)
    
  } else {
    pretty_labels <- data.frame(sensor_name = char_string,
                                sensor_col = sensor_col,
                                sensor_lab = sensor_lab)
  }
  return(pretty_labels)
}

# Plot a single matchup in final format
plot_matchup_single_nm <- function(df, sensor_X, sensor_Y){

  # Prep data
  df_prep <- df |> 
    filter(sensor %in% c(sensor_X, sensor_Y)) |> 
    filter(wavelength >= 380, wavelength <= 700) |> 
    mutate(wavelength_group = cut(wavelength,
                                  breaks = c(350, 400, 450, 500, 550, 600, 650, 700),#, 750, 800),
                                  labels = labels_nm[1:7],
                                  include.lowest = TRUE, right = TRUE), .after = "wavelength") |> 
    pivot_wider(names_from = sensor, values_from = rhow:std_min)
  
  # Get max values
  if(grepl("VIIRS", sensor_Y)){
  # NB: max/min is inverted in absolute values
    max_X <- max(df_prep[paste0("std_min_",sensor_X)], na.rm = TRUE)
    max_Y <- max(df_prep[paste0("std_min_",sensor_Y)], na.rm = TRUE)
  } else {
    max_X <- max(df_prep[paste0("rhow_",sensor_X)], na.rm = TRUE)
    max_Y <- max(df_prep[paste0("rhow_",sensor_Y)], na.rm = TRUE)
  }
  max_axis <- max(max_X, max_Y)
  
  # Get global stats
  x_vec <- df_prep[[paste0("rhow_",sensor_X)]]
  y_vec <- df_prep[[paste0("rhow_",sensor_Y)]]
  
  # Calculate statistics
  df_stats <- base_stats(x_vec, y_vec)
  
  # Get pretty labels
  var_labs <- pretty_label_func("RHOW")
  
  # The first points
  if(grepl("VIIRS", sensor_Y)){
    pl_base <- ggplot(data = df_prep, 
                      aes_string(x = paste0("rhow_",sensor_X), y = paste0("`rhow_",sensor_Y,"`"))) +
      geom_errorbar(aes_string(xmin = paste0("std_min_",sensor_X), xmax = paste0("std_max_",sensor_X)), width = 0.001) +
      geom_errorbar(aes_string(ymin = paste0("`std_min_",sensor_Y,"`"), ymax = paste0("`std_max_",sensor_Y,"`")), width = 0.001) +
      geom_point(aes(colour = wavelength_group), size = 4, alpha = 0.9)
  } else {
    pl_base <- ggplot(data = df_prep, 
                      aes_string(x = paste0("rhow_",sensor_X), y = paste0("`rhow_",sensor_Y,"`"))) +
      geom_point(aes(colour = wavelength_group), size = 2, alpha = 0.7)
  }
  
  # The final plot
  pl_single <- pl_base +
    # Add 1:1 line
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "solid") +
    # Add model II linear models and 95% CI
    ## Bottom CI
    geom_abline(slope = df_stats$Slope_II_low, intercept = df_stats$Slope_II_int_low, 
                colour = "white", alpha = 0.5, linewidth = 1.5, linetype = "solid") +
    geom_abline(slope = df_stats$Slope_II_low, intercept = df_stats$Slope_II_int_low, 
                colour = "grey", linewidth = 1.0, linetype = "dashed") +
    # Mid
    geom_abline(slope = df_stats$Slope_II, intercept = df_stats$Slope_II_int, 
                colour = "white", alpha = 0.5, linewidth = 1.5, linetype = "solid") +
    geom_abline(slope = df_stats$Slope_II, intercept = df_stats$Slope_II_int, 
                colour = "black", linewidth = 1.0, linetype = "dashed") +
    ## Top CI
    geom_abline(slope = df_stats$Slope_II_high, intercept = df_stats$Slope_II_int_high, 
                colour = "white", alpha = 0.5, linewidth = 1.5, linetype = "solid") +
    geom_abline(slope = df_stats$Slope_II_high, intercept = df_stats$Slope_II_int_high, 
                colour = "grey", linewidth = 1.0, linetype = "dashed") +
    # geom_smooth(method = "lm", formula = y ~ x, colour = "white", alpha = 0.5, linewidth = 1.5, linetype = "solid", se = FALSE) +
    # geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 1, linetype = "dashed", se = FALSE) +
    # Add global stats text
    annotate(geom = "text", x = 0, y = max_axis, hjust = 0, vjust = 1, size = 4,
             label = paste0("n: ", df_stats$n,
                            "\nS: ", sprintf("%.2f", df_stats$Slope_II), "±",
                                sprintf("%.2f", abs((df_stats$Slope_II_high - df_stats$Slope_II_low)/2)), 
                            "\nβ: ", sprintf("%.1f", df_stats$Bias_50),
                            "% \nϵ: ", sprintf("%.1f", df_stats$Error_50),"%")) +
    # Make it pretty
    labs(x = paste0(sensor_X,"; ", var_labs$units_lab),
         y = paste0(sensor_Y,"; ", var_labs$units_lab),
         colour = "Wavelength (nm)") +
    scale_colour_manual(values = colour_nm_func(sensor_Y)) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1.0, size = 3))) +
    coord_fixed(xlim = c(0, max_axis), ylim = c(0, max_axis)) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          axis.title.x = element_markdown(size = 12),
          axis.title.y = element_markdown(size = 12),
          axis.text = element_text(size = 10))
  return(pl_single)
}

# Plot data based on wavelength group
# df <- match_filter; sensor_Y <- "AQUA"
plot_global_nm <- function(df, sensor_Y){
  
  # Create sensor and unit labels
  sensor_X_labs <- pretty_label_func("HYPERNETS")
  sensor_Y_labs <- pretty_label_func(sensor_Y)
  var_labs <- pretty_label_func("RHOW")
  colours_nm <- colour_nm_func(sensor_Y)
  
  # Detect Sentinel-3 data and react accordingly
  if(sensor_Y == "S3"){
    df_prep <- df |> 
      pivot_longer(cols = c(S3A, S3B), names_to = "Platform", values_to = "S3") |> 
      na.omit()
  } else {
    df_prep <- df
  }
  
  # Get max values
  max_X <- max(df_prep[sensor_X_labs$sensor_col], na.rm = TRUE)
  max_Y <- max(df_prep[sensor_Y_labs$sensor_col], na.rm = TRUE)
  max_axis <- max(max_X, max_Y)
  
  # Get unique satellite days
  unique_days <- df |> 
    dplyr::select(file_name) |> 
    mutate(date = sapply(str_split(file_name, "_"), "[[", 2)) |> 
    mutate(date = sapply(str_split(date, "T"), "[[", 1)) |> 
    distinct(date) |> 
    mutate(date = as.Date(date, format = "%Y%m%d"))
  
  # Quick for loop per site — must filter df_prep by site so each panel gets its own stats
  # (and its own max_axis, used below to give each facet an independent but still 1:1 x/y scale
  # via ggh4x::facetted_pos_scales() because coord_fixed() doesn't support facet_wrap(scales = "free")
  # directly, so per-panel limits have to be supplied this way instead)
  df_stats <- data.frame()
  for(i in 1:length(unique(df$site_name))){
    site_i <- unique(df$site_name)[i]
    df_prep_i <- df_prep |> filter(site_name == site_i)

    x_vec <- df_prep_i[[sensor_X_labs$sensor_col]]
    y_vec <- df_prep_i[[sensor_Y_labs$sensor_col]]

    df_stats_i <- base_stats(x_vec, y_vec) |>
      mutate(site_name = site_i,
             max_axis = max(x_vec, y_vec, na.rm = TRUE),
             label = paste0("n: ", n, " (", nrow(unique_days), ")",
                            "\nS: ", sprintf("%.2f", Slope_II), "±",
                            sprintf("%.2f", abs((Slope_II_high - Slope_II_low) / 2)),
                            "\nβ: ", sprintf("%.1f", Bias_50),
                            "% \nε: ", sprintf("%.1f", Error_50), "%"))

    df_stats <- rbind(df_stats, df_stats_i)
  }

  # Per-facet axis limits, in the same order facet_wrap draws panels (alphabetical by
  # site_name) so ggh4x::facetted_pos_scales() assigns each one to the right panel
  axis_by_site <- setNames(df_stats$max_axis, df_stats$site_name)[sort(df_stats$site_name)]
  pos_scales_x <- lapply(axis_by_site, function(a) scale_x_continuous(limits = c(0, a)))
  pos_scales_y <- lapply(axis_by_site, function(a) scale_y_continuous(limits = c(0, a)))
  
  # Get number of files used in matchup
  # n_files <- length(unique(df_prep$file_name))

  # Set alpha
  if(sensor_Y == "PACE"){
    point_alpha <- 0.3
  } else {
    point_alpha <- 0.9
  }
  
  # Get pre-determined wavelengths
  # NB: For the moment using all wavelengths for PACE for the plot
  # if(sensor_Y == "PACE"){
    # W_nm <- c(380, 400, 412, 443, 490, 510, 560, 620, 673, 700, 885, 1050)
  # } else {
    W_nm <- W_nm_out(sensor_Y)
  # }

  # Cut PACE colour bands to match colour labels
  if(sensor_Y == "PACE"){
    df_prep <- df_prep |> 
      mutate(wavelength = cut(wavelength,
                              breaks = c(350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1050),
                              labels = names(colours_nm),
                              include.lowest = TRUE, right = TRUE), .after = "wavelength")
  }

  # Determine rows for legend items
  if(sensor_Y %in% c("PACE", "S3A", "S3B", "S3", "S3_all")){
    legend_rows <- 3
  } else {
    legend_rows <- 2
  }

  # Filter dataframe to only plot linear models for chosen wavebands
  if(!(sensor_Y == "PACE")){
  df_sub <- filter(df_prep, wavelength %in% W_nm) |> 
    mutate(wavelength = factor(wavelength,
                               levels = sort(unique(wavelength))))
  } else {
    df_sub <- df_prep
  }

  # Plot
  # NB: facet_wrap(scales = "free") is required here for ggh4x::facetted_pos_scales() below to
  # actually take effect (confirmed empirically, without "free" it silently no-ops). That in
  # turn means coord_fixed() can't be used for the 1:1 x/y aspect ratio (ggplot2/ggh4x reject
  # free facet scales combined with a fixed coord, confirmed empirically too). theme(aspect.ratio
  # = 1) below achieves the same visual effect instead (a square panel), which is a true 1:1
  # relationship because each panel's x and y limits are set equal to each other via
  # pos_scales_x/pos_scales_y above.
  if(sensor_Y == "S3"){
    pl_base <- ggplot(data = df_sub,
                      aes_string(x = sensor_X_labs$sensor_col, y = sensor_Y_labs$sensor_col)) +
      geom_point(aes(colour = wavelength, shape = Platform), size = 2, alpha = point_alpha) +
      scale_colour_manual(values = colours_nm) +
      facet_wrap(~site_name, scales = "free", nrow = 1)
  } else {
    pl_base <- ggplot(data = df_sub,
                      aes_string(x = sensor_X_labs$sensor_col, y = sensor_Y_labs$sensor_col)) +
      geom_point(aes(colour = wavelength), alpha = point_alpha) +
      scale_colour_manual(values = colours_nm) +
      facet_wrap(~site_name, scales = "free", nrow = 1)
  }
  pl_clean <- pl_base +
    # Add 1:1 line
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "solid") +
    # Add model II linear models and 95% CI — use data = df_stats so each row
    # is routed to the matching site_name facet panel
    ## Bottom CI
    geom_abline(data = df_stats, aes(slope = Slope_II_low, intercept = Slope_II_int_low),
                colour = "white", alpha = 0.5, linewidth = 1.5, linetype = "solid") +
    geom_abline(data = df_stats, aes(slope = Slope_II_low, intercept = Slope_II_int_low),
                colour = "grey", linewidth = 1.0, linetype = "dashed") +
    # Mid
    geom_abline(data = df_stats, aes(slope = Slope_II, intercept = Slope_II_int),
                colour = "white", alpha = 0.5, linewidth = 1.5, linetype = "solid") +
    geom_abline(data = df_stats, aes(slope = Slope_II, intercept = Slope_II_int),
                colour = "black", linewidth = 1.0, linetype = "dashed") +
    ## Top CI
    geom_abline(data = df_stats, aes(slope = Slope_II_high, intercept = Slope_II_int_high),
                colour = "white", alpha = 0.5, linewidth = 1.5, linetype = "solid") +
    geom_abline(data = df_stats, aes(slope = Slope_II_high, intercept = Slope_II_int_high),
                colour = "grey", linewidth = 1.0, linetype = "dashed") +
    # Add per-panel stats text via geom_text so site_name routes it to the right facet 
    # y comes from df_stats$max_axis (per-site) rather than the single global max_axis, so the
    # label still sits at the top of each panel's own (now independent) range
    geom_text(data = df_stats, aes(label = label, y = max_axis), x = 0,
              hjust = 0, vjust = 1, size = 4, inherit.aes = FALSE) +
    # Make it pretty
    labs(x = paste0(sensor_X_labs$sensor_lab,"; ", var_labs$units_lab),
         y = paste0(sensor_Y_labs$sensor_lab,"; ", var_labs$units_lab),
         colour = "Wavelength (nm)") +
    # scale_colour_manual(values = colours_nm) +
    guides(colour = guide_legend(nrow = legend_rows, override.aes = list(alpha = 1.0, size = 3))) +
    # Per-panel axis limits (independent per site, x == y range per panel)
    ggh4x::facetted_pos_scales(x = pos_scales_x, y = pos_scales_y) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          aspect.ratio = 1, # square panels -- true 1:1 x/y since each panel's x and y limits match
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.box = "vertical",
          axis.title.x = element_markdown(size = 12),
          axis.title.y = element_markdown(size = 12),
          axis.text = element_text(size = 10))
  # pl_clean
  return(pl_clean)
}

# Takes variable and Y sensor as input to automagically create global scatterplot triptych
# cut_legend = "cut" strips the legend so several panels can be stacked under one shared legend
# sensor_Y = "S3B"; cut_legend = "cut"
# sensor_Y = "AQUA"; cut_legend = "cut"
# sensor_Y = "S3"; cut_legend = "no"
# sensor_Y = "PACE"; cut_legend = "no"
global_scatterplot <- function(sensor_Y, cut_legend = "no"){

  # Continue with satellite versions if necessary
  if(sensor_Y  == "AQUA"){
    sensor_Z <- "MODIS"
  } else if(sensor_Y %in% c("PACE")){
    sensor_Z <- "OCI"
  } else if(sensor_Y %in% c("SNPP", "JPSS1", "JPSS2")){
    sensor_Z <- "VIIRS"
  } else if(sensor_Y %in% c("S3A", "S3B", "S3")){
    sensor_Z <- "OLCI"
  } else {
    sensor_Z <- sensor_Y
  }
  
  # Get filestub based on sensor_Z
  filestub <- paste0("_",sensor_Z,".csv")
  
  # Load individual matchup results to filter file list and for further use
  match_base_details <- read_csv(paste0("output/matchup_stats_RHOW",filestub), show_col_types = FALSE) |> 
    dplyr::select(site_name, file_name) |> distinct()
  
  # Load outliers to screen them from being plotted
  outliers_all <- read_csv("meta/satellite_outliers.csv", show_col_types = FALSE) |> distinct()
  
  # Load data based on in situ comparisons or not
  # NB: site list picked up automatically via available_sites() -- THFR is included
  # once its data folder exists on disk, no code change needed here.
  print("Loading matchups")
  if(sensor_Y == "S3"){
    site_list <- unique(unlist(lapply(c("S3A", "S3B"), available_sites)))
    ply_folders <- expand.grid(site_name = site_list, sat_name = c("S3A", "S3B"))
    # match_base_1 <- bind_rows(load_matchups_folder("S3A", long = TRUE),
    #                           load_matchups_folder(e_name, "S3B", long = TRUE))
  } else {
    # match_base_1 <- load_matchups_folder(site_name, sensor_Y, long = TRUE)
    site_list <- available_sites(sensor_Y)
    ply_folders <- expand.grid(site_name = site_list, sat_name = sensor_Y)
  }
  
  # Load all folders
  match_base <- purrr::pmap_dfr(ply_folders, load_matchups_folder, long = TRUE) |> 
    right_join(match_base_details, by = join_by(site_name, file_name))

  # Filter out outliers
  match_filter <- match_base[!match_base$file_name %in% outliers_all$file_name,]

  # Remove any erroneuosly high values before plotting
  match_filter <- match_filter |> filter(Hyp <= 1)

  # Create the figure
  print("Creating figure, saving, and exiting")
  match_fig <- plot_global_nm(match_filter, sensor_Y)

  # Save the individual figure
  ggsave(paste0("figures/global_scatter_RHOW_",sensor_Y,".png"), match_fig, width = 16, height = 5)

  # Remove the legend when this panel will be stacked beneath another with a shared legend
  # Then return it (invisibly) for reuse by global_scatterplot_stack()
  if(cut_legend == "cut") match_fig <- match_fig + theme(legend.position = "none")
  invisible(match_fig)
}

# Stack the per-sensor global scatterplots for a sensor family into one composite figure
# NB: Requires that the matchup and global stats CSVs for sensor_Z already exist (see process_sensor())
# sensor_Z = "MODIS"
global_scatterplot_stack <- function(sensor_Z){

  # Get the sensor_Y platforms that belong to this sensor family
  sensor_Y_list <- unique(sensor_grid(sensor_Z)$sensor_Y)
  # NB: Disabling S3 all for the moment
  # if(sensor_Z == "OLCI"){ 
  #   sensor_Y_list <- c("S3A", "S3B", "S3")
  #   print("Added S3 to Sat names")
  # }
  sensor_count <- length(sensor_Y_list)

  # Build one panel per platform
  # NB: The legend is cut from every panel but the last so it appears once, at the bottom of the stack
  fig_list <- vector("list", sensor_count)
  for(i in seq_len(sensor_count)){
    fig_list[[i]] <- global_scatterplot(sensor_Y_list[i], cut_legend = ifelse(i < sensor_count, "cut", "no"))
  }

  # Give the legend-bearing (last) panel extra relative height to fit the legend
  panel_heights <- c(rep(1, sensor_count - 1), 1.15)

  # Stack panels vertically and exit
  fig_stack <- ggpubr::ggarrange(plotlist = fig_list, ncol = 1, nrow = sensor_count, heights = panel_heights) +
    ggpubr::bgcolor("white") + ggpubr::border("white", size = 2)
  ggsave(paste0("figures/global_scatter_RHOW_",sensor_Z,".png"), fig_stack, width = 16, height = 5 * sensor_count)
}

