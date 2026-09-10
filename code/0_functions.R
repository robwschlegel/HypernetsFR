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
library(ggh4x) # For independent per-facet axis limits
library(sf) # For the clean_water_sf point-in-polygon pixel filter


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

# Bucket PACE/OCI's continuous wavelength (350-1150 nm) into the same broad bands used for its
# colour legend above, since it has no small fixed set of nominal bands like the other sensors
# (W_nm_out("PACE") intentionally stays a continuous 350:1150 range -- it's consumed elsewhere,
# e.g. srf_band_map(), as a plain nominal-wavelength lookup, not a bucketed set). Shared by
# plot_global_nm(), global_scatterplot_waveband(), and code/5_figures.R's Figure 11 matrix so the
# binning stays consistent everywhere it's used.
pace_waveband_bucket <- function(wavelength){
  breaks <- c(350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1050)
  as.character(cut(wavelength, breaks = breaks, labels = names(colour_nm_func("PACE")),
                    include.lowest = TRUE, right = TRUE))
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

# Load all files in a given folder
load_matchups_folder <- function(site_name, sat_name, long = FALSE){
  
  # Create file path
  folder_path <- file_path_build(site_name, sat_name)
  
  # List all files in directory
  file_list <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)
  
  # Remove stats output files
  file_list_clean <- file_list[!grepl("all|global", file_list)]

  # Graceful skip when a site/sensor combination has no exported matchup files at all -- without
  # this, the mutate(.before = wavelength/sensor) below fails on the resulting zero-column empty
  # tibble ("Column `wavelength` doesn't exist"), aborting the whole calling pmap_dfr()/figure.
  if(length(file_list_clean) == 0){
    message("No matchup files for ", site_name, " ", sat_name, " -- skipping")
    return(tibble())
  }

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

# Check the amount of variance in satellite files and return a message if there is an issue
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
# THFR_NE (added 2026-08-28), THFR_poly (added 2026-08-28), and THFR_pixel (added 2026-09-03) are
# independent derived sites alongside MAFR/THFR -- NE-quadrant, clean-water-polygon-filtered, and
# no-spatial-filter (per-pixel-QC-only) THFR re-analyses respectively (see the TEMPORARY section in
# this file and meta/pixel_explore.R), each intentionally kept as its own site for direct
# comparison rather than replacing THFR anywhere in the pipeline. MAFR_pixel (added 2026-09-07) is
# the analogous no-spatial-filter, .db-direct reconstruction of MAFR (db_export_matchups_mafr_pixel()),
# built from mafr_2024.db/mafr_2025.db rather than the real Hypernets_matchups MAFR files --
# likewise its own site, not a replacement for MAFR anywhere in the pipeline. THFR_raw/MAFR_raw
# (added 2026-09-07) go one step further: the same 3x3-box reconstruction with NO pixel-level QC
# gates at all (db_export_matchups_thfr_raw()/db_export_matchups_mafr_raw(), apply_pixel_qc =
# FALSE on db_export_matchups_site()), meant to reproduce Hypernets_matchups' own unfiltered
# aggregation as closely as possible -- also their own sites, not yet a replacement for THFR/MAFR
# anywhere, pending validate_derived_site() review.
available_sites <- function(sat_name){
  candidate_sites <- c("MAFR", "MAFR_pixel", "MAFR_raw", "THFR", "THFR_pixel", "THFR_NE", "THFR_poly", "THFR_raw")
  site_present <- vapply(candidate_sites, function(s) dir.exists(file_path_build(s, sat_name)), logical(1))
  sites_found <- candidate_sites[site_present]
  if(length(sites_found) == 0) stop(paste0("No site data found on disk for sensor: ", sat_name))
  return(sites_found)
}

# Site-specific matchup time-window limit (minutes), added 2026-07-10.
# NB: unlike Doxaran et al. 2024 (who additionally varied the *spatial* matchup criterion per site --
# a 3x3-pixel box at Berre vs. nearest-pixel-only at Gironde), this pipeline keeps ONE spatial rule
# (nearest pixel + a per-sensor dist_limit = 2x sensor_resolution_km(sensor_Y), see process_sensor();
# a flat 10 km ceiling was used until 2026-09-03, and a 3x multiplier until 2026-09-04) for every site.
# MAFR raised from 15 to 30 min on 2026-09-10, after an empirical check (pooled across sensors,
# dist-QC-passed matchups only, from output/matchup_noQC_stats_RHOW_*.csv) found the 15->30 min
# widening nearly doubles usable MAFR matchups (571 -> 1062) for a small cost in typical error
# (median Error_50 31.96% -> 33.03%, mean 48.59% -> 50.02%). MAFR and THFR now share the same 30 min
# window, mirroring Doxaran et al. 2024's Berre value rather than their tighter Gironde value -- see
# manuscript/roadmap.md's "site-specific matchup criteria" item for the full comparison table and
# literature context (Doxaran et al. 2024 used Gironde's own tighter +/-15 min, justified by >20%
# variability beyond +/-30 min there, so this is a considered departure from that precursor's choice,
# not an oversight). THFR's own 30 min ceiling was separately found to already be the maximum
# diff_time present anywhere in the raw matchup pool for either site -- Hypernets_matchups itself
# appears to only ever generate candidate pairs within a +/-30 min window upstream of this gate.
site_diff_time_limit <- function(site_name){
  dplyr::case_when(
    grepl("MAFR", site_name) ~ 30,
    grepl("THFR", site_name) ~ 30,
    TRUE ~ 30 # fallback for any site name containing neither "MAFR" nor "THFR"
  )
}

# Site-specific RHOW plausibility ceiling, added 2026-09-03. Derived from the empirical max
# HYPERNETS RHOW observed at each site across every sensor family/waveband (a full scan of every
# raw matchup .csv for MAFR/THFR found a natural max of 0.174 at MAFR, OLCI 681 nm, and 0.0294 at
# THFR, VIIRS SNPP 551 nm -- see chat log 2026-09-03), with headroom applied (~45% at MAFR, ~70% at
# THFR) so the ceiling doesn't clip genuine turbid-water extremes while still catching retrieval
# failures (sun glint, land/adjacency contamination, thin cloud) that push RHOW well outside the
# physically plausible range for each site. Used as a QC gate on BOTH the Hyp and satellite RHOW
# values (not just Hyp -- unlike the pre-existing Hyp <= 1 gate in global_scatterplot(), which only
# catches gross parsing failures like netCDF fill values, e.g. Bug: three matchup files found with
# Hyp == 9.969209968386869e+36, the standard netCDF double fill value).
site_rhow_limit <- function(site_name){
  dplyr::case_when(
    grepl("MAFR", site_name) ~ 0.25,
    grepl("THFR", site_name) ~ 0.05,
    TRUE ~ 0.05 # fallback for any site name containing neither "MAFR" nor "THFR"
  )
}

# Nominal per-sensor pixel resolution (km at nadir), added 2026-09-03. Used to derive a per-pixel
# distance QC gate (2x the sensor's own resolution, tightened from 3x on 2026-09-04 -- e.g. PACE's
# 1.2 km resolution gives a 2.4 km ceiling) that is much tighter than the fixed 10 km dist_limit used
# at the whole-matchup level (process_sensor()/db_export_matchups_site()), since a pixel several resolution-cells away from
# the station is increasingly unlikely to represent the water actually seen by HYPERNETS.
sensor_resolution_km <- function(sensor_Y){
  if(sensor_Y %in% c("S3A", "S3B")){
    0.3   # Sentinel-3 OLCI Full Resolution, ~300 m at nadir
  } else if(sensor_Y == "AQUA"){
    1     # MODIS-Aqua ocean-colour bands, 1 km at nadir
  } else if(sensor_Y %in% c("SNPP", "JPSS1", "JPSS2")){
    0.75  # VIIRS ocean-colour EDR, ~750 m at nadir
  } else if(sensor_Y == "PACE"){
    1.2   # PACE OCI, 1.2 km x 1.2 km ground sample footprint at nadir (corrected 2026-09-10 from
          # 1 km, which was OCI's design-goal GSD, not its specified/delivered resolution -- see
          # https://pace.oceansciences.org/oci.htm ("1.2 km x 1.2 km ground sample footprint at
          # the center of the scan") vs. https://pace.oceansciences.org/requirements.htm ("a GSD
          # of 1km enables global science", a design rationale, not the OCI-60 GSD requirement
          # itself). meta/pixel_explore_output/pixel_spacing_summary.csv's empirical PACE nearest-
          # neighbour spacing (~1.68-1.71 km) is 1.40-1.42x this corrected value, not 1.68-1.71x.
  } else {
    stop(paste0("Incorrect value for 'sensor_Y' : ", sensor_Y))
  }
}

# Minimum valid-pixel count per matchup-wavelength, added 2026-09-04. Mirrors Hypernets_matchups'
# own rule for the real MAFR/THFR data (>= 6 of the 9 pixels in its native 3x3 box must be valid
# before it will compute a daily pixel-box mean) -- the DB-derived sites (THFR_NE/THFR_poly/
# THFR_pixel) recompute pixel-box aggregates themselves in write_matchup_csv_ne(), so this closes
# that gap for them. Different minimums per site reflect how many pixels are ever physically
# available under each site's own spatial restriction (THFR_pixel: full inner 3x3, no further
# restriction; THFR_NE: NE-quadrant of the inner 3x3; THFR_poly: hand-drawn clean-water polygon --
# S3A/S3B's finer 300m resolution puts more candidate pixels in a given area than the coarser
# sensors, hence the higher S3A/S3B minimums at both derived sites).
site_pixel_min <- function(site_name, sensor_Y){
  if(site_name == "THFR_pixel"){
    6
  } else if(site_name == "THFR_NE"){
    if(sensor_Y %in% c("S3A", "S3B")) 3 else 1
  } else if(site_name == "THFR_poly"){
    1 # TEMPORARY (2026-09-04): reduced from 5 (S3A/S3B) / 3 (others) -- the clean-water polygon is
      # too small for MODIS/VIIRS/OCI to ever clear 3 pixels, leaving those sensors at zero data.
      # Revisit alongside the other temporary QC relaxations next week.
  } else if(site_name == "MAFR_pixel"){
    6
  } else if(site_name %in% c("THFR_raw", "MAFR_raw")){
    6 # mirrors the real tool's own documented ">= 6 of 9" rule (see this function's own docstring
      # above) -- with no other pixel QC applied under apply_pixel_qc = FALSE, this is the only gate
      # standing between "raw" and a matchup built from a single stray pixel
  } else {
    stop(paste0("No pixel-count minimum defined for site_name: ", site_name))
  }
}

# Select, per calendar day, the single closest-in-time HYPERNETS-scan-vs-satellite-overpass
# matchup, for use at the global-stats stage. Added 2026-07-10 as daily_average_matchups()
# (averaging same-day matchups together); replaced 2026-09-03 with closest-match selection instead
# -- HYPERNETS scans and satellite overpasses on the same day are not independent measurements of
# the same water, so rather than blend them, only the temporally-nearest pairing is retained per
# day. This still addresses the "not all matchups are independent of one another" caveat raised in
# the Tara "in review" paper's Conclusion, just via selection rather than averaging.
#
# Date-extraction logic reuses the exact convention used in global_scatterplot() -- i.e. split
# file_name on "_", take the 2nd element, split that on "T", take the 1st element -- confirmed to
# work for S3A/S3B/JPSS1/JPSS2/SNPP/AQUA/PACE naming patterns at both MAFR and THFR.
#
# "Within the site time-difference limit" is NOT re-enforced here -- df is already restricted (via
# global_stats()'s file_list_no_out, built from matchup_stats_RHOW_<sensor_Z>.csv, which only ever
# contains files that already passed process_sensor()'s diff_time <= diff_time_limit gate) to
# already-qualifying candidates, so "closest" selection only ever chooses among matchups that
# already satisfy site_diff_time_limit().
#
# Selection happens once per (day, file_name), never independently per (day, wavelength) row --
# df has one row per (file_name, wavelength), and diff_time is constant across all wavelength-rows
# of a given file_name, but resolving a tie independently per wavelength could in principle pick a
# DIFFERENT winning file_name for different wavelengths on the same day, silently splicing together
# two different matchup events. This picks exactly one winning file_name per day first (with a
# deterministic tie-break: smallest diff_time, then smallest dist, then alphabetically-first
# file_name), then keeps every row of that one file.
#
# df: expects the long-format data.frame produced by load_matchup_long()/load_matchups_folder(long =
#     TRUE), i.e. one row per file_name x wavelength, already filtered to wavelength %in% W_nm, with
#     diff_time and dist already joined on by file_name (see global_stats()).

daily_closest_matchup <- function(df, site_name){
  df <- df |>
    mutate(match_date = sapply(str_split(file_name, "_"), "[[", 2),
           match_date = sapply(str_split(match_date, "T"), "[[", 1),
           match_date = as.Date(match_date, format = "%Y%m%d"))

  # One row per (day, file_name); pick the single winning file per day
  winners <- df |>
    dplyr::select(match_date, file_name, diff_time, dist) |>
    distinct() |>
    arrange(match_date, diff_time, dist, file_name) |>
    slice_head(n = 1, by = match_date)

  df |>
    filter(file_name %in% winners$file_name) |>
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
  # n_clean < 3 guards the clearest failure mode (an exactly-2-point fit has zero residual degrees
  # of freedom), but lmodel2()'s nperm = 99 permutation test can still fail internally
  # (summary(<internal fit>)$coefficients[2, 1]: subscript out of bounds) on small-n data with tied
  # values, where some permutation happens to produce a rank-deficient fit -- not something worth
  # enumerating in advance, so the call itself is wrapped in tryCatch() below and falls back to the
  # same all-NA Model II result used for the zero-variance/n<3 cases. First actually reached
  # 2026-09-07 by the unfiltered "_raw" sites (apply_pixel_qc = FALSE): without the pixel-level
  # RHOW-ceiling/negative-value pre-filtering, more implausible values and smaller/more repetitive
  # n reach this function than the QC'd pipeline ever produced.
  if(n_clean < 3 | length(unique(round(x_clean, 8))) == 1 | length(unique(round(y_clean, 8))) == 1){
    model_II_intercept <- NA; model_II_slope <- NA; model_II_p_perm <- NA
    model_II_slope_lo <- NA; model_II_slope_hi <- NA; model_II_int_lo <- NA
    model_II_int_hi <- NA; model_II_slope_bias_sig <- NA; model_II_int_bias_sig <- NA
  } else {
    model_II_vals <- tryCatch({
      model_II_fit <- lmodel2::lmodel2(y_clean ~ x_clean,
                                       range.y = "relative",
                                       range.x = "relative",
                                       nperm = 99)
      # Extract results for chosen method
      model_II_method_choice <- "SMA" # Symetrical Major Axis
      model_II_results <- model_II_fit$regression.results |>
        filter(Method == model_II_method_choice)
      model_II_ci <- model_II_fit$confidence.intervals |>
        filter(Method == model_II_method_choice)
      list(intercept = model_II_results$Intercept, slope = model_II_results$Slope,
           p_perm = model_II_results$`P-perm (1-tailed)`,
           slope_lo = model_II_ci$`2.5%-Slope`, slope_hi = model_II_ci$`97.5%-Slope`,
           int_lo = model_II_ci$`2.5%-Intercept`, int_hi = model_II_ci$`97.5%-Intercept`)
    }, error = function(e){
      message("base_stats(): lmodel2() failed for n = ", n_clean, " (", conditionMessage(e), ") -- Model II fields set to NA")
      NULL
    })
    if(is.null(model_II_vals)){
      model_II_intercept <- NA; model_II_slope <- NA; model_II_p_perm <- NA
      model_II_slope_lo <- NA; model_II_slope_hi <- NA; model_II_int_lo <- NA
      model_II_int_hi <- NA; model_II_slope_bias_sig <- NA; model_II_int_bias_sig <- NA
    } else {
      model_II_intercept <- model_II_vals$intercept
      model_II_slope <- model_II_vals$slope
      model_II_p_perm <- model_II_vals$p_perm
      model_II_slope_lo <- model_II_vals$slope_lo
      model_II_slope_hi <- model_II_vals$slope_hi
      model_II_int_lo <- model_II_vals$int_lo
      model_II_int_hi <- model_II_vals$int_hi
      # Determine significance
      model_II_slope_bias_sig <- model_II_slope_lo >= 1 || model_II_slope_hi <= 1
      model_II_int_bias_sig <- model_II_int_lo >= 0 || model_II_int_hi <= 0
    }
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

# Function that interrogates each matchup file to produce the needed output for all following comparisons
process_matchup_file <- function(file_path){

  # Load the mean data
  df_mean <- load_matchup_mean(file_path)

  # Sensors to be compared -- always exactly 2 in practice (Hyp + one satellite; Hyp_nosc is
  # already dropped upstream in load_matchup_mean()).
  sensors <- unique(df_mean$sensor)
  if(length(sensors) != 2 || !("Hyp" %in% sensors)){
    stop(paste0("process_matchup_file() expects exactly 2 sensors incl. 'Hyp', got: ",
                paste(sensors, collapse = ", "), " in file: ", basename(file_path)))
  }
  sensor_X_name <- "Hyp"
  sensor_Y_name <- setdiff(sensors, "Hyp")

  # Get data.frame for the matchup
  df_sensor_sub <- df_mean |>
    mutate(dateTime = as.POSIXct(paste(day, time), format = "%Y%m%d %H%M%S", tz = "Europe/Paris"), .before = "latitude", .keep = "unused")

  # get distances
  hav_dist <- round(distHaversine(df_sensor_sub[c("longitude", "latitude")])/1000, 2) # distance in km

  # Time differences
  time_diff <- round(as.numeric(abs(difftime(df_sensor_sub$dateTime[[1]],
                                             df_sensor_sub$dateTime[[2]], units = "mins"))))

  # Row indices for the two sensors in df_sensor_sub, needed for the per-row lon/lat/dateTime below
  idx_X <- which(df_sensor_sub$sensor == sensor_X_name)
  idx_Y <- which(df_sensor_sub$sensor == sensor_Y_name)

  # Melt it for additional stats
  df_sensor_long <- df_sensor_sub |>
    pivot_longer(cols = matches("1|2|3|4|5|6|7|8|9"), names_to = "wavelength", values_to = "value") |>
    mutate(wavelength = as.numeric(wavelength)) |>
    na.omit()

  # Widen for use with stats function
  df_sensor_wide <- df_sensor_long |>
    dplyr::select(-c(dateTime, longitude, latitude)) |>
    pivot_wider(names_from = sensor, values_from = value) |>
    na.omit()

  # get vectors
  x_vec <- df_sensor_wide[[sensor_X_name]]
  y_vec <- df_sensor_wide[[sensor_Y_name]]

  # Base stats
  df_stats <- base_stats(x_vec, y_vec)

  # Create data.frame of results
  df_results <- df_stats |>
    mutate(sensor_X = sensor_X_name,
           sensor_Y = sensor_Y_name,
           lon_X = df_sensor_sub$longitude[[idx_X]],
           lat_X = df_sensor_sub$latitude[[idx_X]],
           lon_Y = df_sensor_sub$longitude[[idx_Y]],
           lat_Y = df_sensor_sub$latitude[[idx_Y]],
           dist = hav_dist,
           dateTime_X = df_sensor_sub$dateTime[[idx_X]],
           dateTime_Y = df_sensor_sub$dateTime[[idx_Y]],
           diff_time = time_diff, .before = "n") |>
    mutate(file_name = basename(file_path), .before = sensor_X)
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
           n_w_nm = n_match, .before = "n") |> 
    dplyr::rename(n_w_nm_clean = n)
  return(df_XY)
}

# Global stats per matchup wavelength
# site_name = "MAFR"; sensor_Y = "S3A"
# site_name = "MAFR"; sensor_Y = "S3_all"
# site_name = "THFR"; sensor_Y = "PACE"
# site_name = "MAFR"; sensor_Y = "SNPP"
# site_name = "MAFR"; sensor_Y = "JPSS1"
# site_name = "MAFR"; sensor_Y = "AQUA"
#
# select_daily controls whether daily_closest_matchup() collapses same-day matchups to the single
# closest-in-time one (the only behaviour global_stats() itself ever uses, see its thin wrapper
# below) or leaves every QC-passed matchup as an independent point. Kept as a private toggle on
# this impl function -- not on global_stats()'s own public signature -- purely so
# code/3_sensitivity.R's before/after comparison can still get the "before" (no selection) state
# without a stale daily_average-style parameter re-appearing on the production entry point. There
# is no longer any averaging behind either setting: "select_daily = FALSE" was never
# daily_average_matchups() (removed 2026-09-03), it is simply "no day-level selection at all".
global_stats_impl <- function(site_name, sensor_Y, select_daily = TRUE){
  
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
  
  # Load individual matchup results to filter file list and for further use,
  # use the single closest-in-time matchup per day without recomputing anything
  match_base_details <- read_csv(paste0("output/matchup_stats_RHOW",filestub), show_col_types = FALSE) |>
    filter(.data$site_name == .env$site_name) |>
    dplyr::select(file_name, diff_time, dist) |> distinct()
  if(nrow(match_base_details) == 0){
    message("No QC-passed files for site ", site_name, ", sensor_Y ", sensor_Y, " -- skipping")
    return(tibble())
  }

  # Filter accordingly
  # NB: This creates the list of valid matchups after screening for spatiotemporal range
  file_list_clean <- file_list[basename(file_list) %in% match_base_details$file_name]

  # Get outlier lists, restricted to this site (see NB above)
  outliers_sat <- read_csv("meta/satellite_outliers.csv", show_col_types = FALSE) |>
    filter(.data$site_name == .env$site_name)

  # Remove outlier files
  # NB: This creates the list of valid matchups after screening for outliers in the single matchup QC process
  file_list_no_out <- file_list_clean[!basename(file_list_clean) %in% outliers_sat$file_name]
  if(length(file_list_no_out) == 0){
    message("No files passed QC (post-outlier-screen) for ", site_name, " ", sensor_X, " ", sensor_Y, " -- skipping")
    return(tibble())
  }
  
  # Load data
  match_base <- furrr::future_map_dfr(file_list_no_out, load_matchup_long, .options = furrr_options(seed = TRUE))
  
  # Melt if S3_all -- only pivot whichever of S3A/S3B columns actually survived QC for this site.
  # Both are normally present, but a site can end up with only one satellite's files passing QC
  # (e.g. every S3B file failing the outlier/distance gate for a given site), in which case
  # match_base never gets an "S3B" column at all and pivot_longer(S3A:S3B, ...) errors.
  if(sensor_Y == "S3_all"){
    s3_cols <- intersect(c("S3A", "S3B"), colnames(match_base))
    match_base <- match_base |>
      pivot_longer(all_of(s3_cols), names_to = "name", values_to = "S3_all") |>
      dplyr::select(-name) |>
      filter(!is.na(S3_all))
  }
  
  # Get pre-determined wavelengths
  W_nm <- W_nm_out(sensor_Y)
  
  # Filter data.frame accordingly
  match_base_filt <- filter(match_base, wavelength %in% W_nm) #|>
    # mutate(wavelength_idx = wavelength, .before = wavelength)

  # Remove any erroneously high (physically impossible) Hyp values before computing statistics --
  # mirrors the equivalent filter already applied in global_scatterplot() (see Bug 3,
  # manuscript/upstream-data-bugs.md), previously missing here (roadmap.md open item 11), which let
  # figures and statistics CSVs disagree whenever such a value was present.
  match_base_filt <- match_base_filt |> filter(Hyp <= 1)

  # Attach diff_time/dist (per file_name) needed by daily_closest_matchup() below
  match_base_filt <- match_base_filt |> left_join(match_base_details, by = "file_name")

  # Day-level collapsing: keep only the single closest-in-time matchup per day. See select_daily's
  # own docstring above -- global_stats() itself always wants this; only the sensitivity comparison
  # ever asks for the un-collapsed alternative.
  if(select_daily){
    match_base_filt <- daily_closest_matchup(match_base_filt, site_name)
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

# Public entry point: always the single closest-in-time matchup per day (daily_closest_matchup()).
# There is no averaging option any more -- see global_stats_impl()'s select_daily docstring.
global_stats <- function(site_name, sensor_Y){
  global_stats_impl(site_name, sensor_Y, select_daily = TRUE)
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

  # Graceful skip when a site/sensor combination has no exported matchup files at all (e.g. a
  # derived site's spatial + per-pixel QC gates left nothing for a given sensor) -- without this,
  # the mutate() below fails with "Column `file_name` doesn't exist" on the resulting zero-column
  # empty tibble, aborting the whole process_sensor() run for every other site/sensor combination.
  if(length(file_list) == 0){
    message("No exported matchup files for ", site_name, " ", sensor_Y, " -- skipping")
    return(tibble())
  }

  # Initialise results data.frame
  df_results <- furrr::future_map_dfr(file_list, process_matchup_file, .options = furrr_options(seed = TRUE)) |>
    mutate(site_name = site_name, .after = file_name)

  # Exit
  return(df_results)
}

# Process multiple folders based on request
# sensor_Z = "MODIS"; stat_choice = "matchup"
# sensor_Z = "OLCI"; stat_choice = "global"
process_sensor <- function(sensor_Z, stat_choice = "matchup"){
  
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
    # Set time and distance limits.
    # TEMPORARY (2026-09-04): the resolution-based distance ceiling (2x sensor_resolution_km(),
    # itself tightened from 3x earlier the same day) turned out to eliminate essentially all MAFR
    # matchups -- reverted to a flat 5 km ceiling as a stopgap so results can ship today. The
    # resolution-based approach needs revisiting (properly tuned, not just 2x/3x guesses) --
    # see manuscript/upstream-data-bugs.md and revisit next week.
    proc_res <- proc_res |>
      mutate(diff_time_limit = site_diff_time_limit(site_name), .after = diff_time) |>
      mutate(dist_limit = 5, .after = dist)
    # Enforce time and distance constraint
    proc_res_clean <- proc_res |>
      filter(diff_time <= diff_time_limit) |>
      filter(dist <= dist_limit)
    write_csv(proc_res_clean, paste0("output/matchup_stats_RHOW_",sensor_Z,".csv"))
    # Save the matchups removed this way
    proc_res_unclean <- proc_res[!proc_res$file_name %in% proc_res_clean$file_name,]
    write_csv(proc_res_unclean, paste0("output/matchup_noQC_stats_RHOW_",sensor_Z,".csv"))
  } else {
    proc_res <- furrr::future_pmap_dfr(ply_grid, global_stats, .options = furrr_options(seed = TRUE))
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

# Expands each measure_info row's [data_id, data_id+data_count-1] block of measure_data_rhow ids
# into one row per (measure_info_id, data_row_id) -- shared by db_load_spectra() and
# db_load_spectra_pixel(). HYPERNETS scans (radiometer_id == 1) are capped to their first 3 raw
# replicate rows (by ascending id): data_count is always 6, but replicates 4-6 have been confirmed
# (manuscript/upstream-data-bugs.md Bug 12) to sometimes belong to a second, different measurement
# condition bundled under the same measure_info id, and the real Hypernets_matchups tool's own
# Hyp 1/2/3 values match only replicates 1-3 to full floating-point precision. Satellite pixel-box
# rows are never truncated -- their whole data_count block IS the pixel box.
# info_tbl: data.frame with columns id, data_id, data_count, radiometer_id (from measure_info)
db_expand_measure_info_ids <- function(info_tbl){
  info_tbl |>
    mutate(data_count_used = if_else(radiometer_id == 1L, pmin(data_count, 3L), data_count),
           data_id_end = data_id + data_count_used - 1) |>
    rowwise() |>
    reframe(measure_info_id = id, data_row_id = seq(data_id, data_id_end)) |>
    ungroup()
}

# Aggregate the raw measure_data_rhow rows belonging to one or more
# measure_info entries into one spectrum (long format) per measure_info.id.
# info_tbl: data.frame with columns id, data_id, data_count, radiometer_id (from measure_info)
db_load_spectra <- function(con, info_tbl, agg_method = c("mean", "center")){
  agg_method <- match.arg(agg_method)

  info_map <- db_expand_measure_info_ids(info_tbl)

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
  info_map <- db_expand_measure_info_ids(info_tbl)

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
  info_query <- paste0("SELECT id, data_id, data_count, radiometer_id, day, time, latitude, longitude, qc FROM measure_info WHERE id IN (",
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
  info_query <- paste0("SELECT id, data_id, data_count, radiometer_id, day, time, latitude, longitude, qc FROM measure_info WHERE id IN (",
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

# Helper for db_export_matchups_site() below. Builds and writes one matchup's raw RHOW CSV
# entirely from db_matchup_pixels() output -- no dependency on a real on-disk Hypernets_matchups
# file existing for this matchup. Writes only the two rows anything downstream actually reads
# (load_matchup_mean()'s "weighted"/"rhow" pair, sat_var_check()'s variability_centered): the
# Hyp value is already correct by the time it reaches df_matchup (SRF-convolved for RSR sensors,
# raw hyperspectral for PACE -- both computed by db_matchup_pixels(), dispatched on Hyp_method
# below), and the satellite value is the mean of whichever pixels survived every upstream QC gate
# and the site's pixel_filter_fn(). Full hyperspectral raw rows, Hyp_nosc, and separate std_max/
# std_min rows present in real files are NOT reproduced -- confirmed unused by every consumer in
# this pipeline (load_matchup_mean(), load_matchup_long(), sat_var_check()) and Hyp_nosc isn't
# derivable from the .db at all (see 0_functions.R's db_matchup_pixels() docs).
# df_matchup: one matchup_id's surviving satellite pixel rows for one sensor_Y (post QC gates and
#     pixel_filter_fn(), from db_export_matchups_site())
# sat_radiometer_id: the satellite's radiometer.id (db_radiometer_id()), written into the
#     satellite row's radiometer_id column
# pixel_min: minimum number of contributing pixels required for a wavelength's mean/SD/CV to be
#     written at all (site_pixel_min()) -- mirrors Hypernets_matchups' own >=6-of-9-pixels rule.
# Returns list(written, removal_log) -- removal_log carries the wavebands dropped for falling
# below pixel_min, in the same schema as db_export_matchups_site()'s own gate log (see
# log_gate_drop()), for the caller to fold into the site's combined pixel-removal CSV.

# Computes the exact filename write_matchup_csv_db() writes a matchup's CSV to. Used internally
# by write_matchup_csv_db() itself, and available for joining meta/<site>_pixel_removals.csv's
# (matchup_id, db_source, dateTime_Hyp, dateTime_sat) key against today's still-file-based
# output/matchup_stats_RHOW_*.csv (file_name = basename(file_path), set in process_matchup_file())
# -- not needed once/if matchup stats themselves become .db-native and drop the intermediate CSV.
db_matchup_filename <- function(sensor_Y, dateTime_Hyp, dateTime_sat, tag){
  hyp_ts <- format(dateTime_Hyp, "%Y%m%dT%H%M%S", tz = "Europe/Paris")
  sat_ts <- format(dateTime_sat, "%Y%m%dT%H%M%S", tz = "Europe/Paris")
  paste0(sensor_Y, "_", sat_ts, "_vs_HYPERNETS_", hyp_ts, "_RHOW_", tag, ".csv")
}

write_matchup_csv_db <- function(df_matchup, sensor_Y, sat_radiometer_id, out_dir, pixel_min, tag,
                                  sensor_Z, site_name){
  # NB: df_matchup's own `site_name` column (from db_matchup_pixels()) is always the .db's base
  # site (e.g. "THFR"), not the derived site (e.g. "THFR_NE") -- site_name/sensor_Z must be passed
  # in explicitly by the caller rather than read off df_matchup.
  matchup_id <- df_matchup$matchup_id[1]
  match_date <- df_matchup$match_date[1]
  db_source <- df_matchup$db_source[1]

  is_pace_style <- unique(df_matchup$Hyp_method) == "nearest_nm"
  data_type_label <- if(is_pace_style) "rhow" else "weighted"

  # Per-wavelength satellite pixel-box stats for this one matchup
  wl_stats <- df_matchup |>
    summarise(mean_val = mean(.data[[sensor_Y]], na.rm = TRUE),
              sd_val = sd(.data[[sensor_Y]], na.rm = TRUE),
              n_used = n(),
              pixel_pos_wl = paste(sort(unique(pixel_pos)), collapse = ";"),
              .by = wavelength)

  removal_log <- tibble()
  below_min <- wl_stats |> filter(n_used < pixel_min)
  if(nrow(below_min) > 0){
    removal_log <- below_min |>
      transmute(site_name = site_name, sensor_Z = sensor_Z, sensor_Y = sensor_Y, db_source = db_source,
                matchup_id = matchup_id, match_date = match_date,
                gate = "pixel_min_count", pixel_pos = pixel_pos_wl, wavelength = wavelength,
                value = n_used, threshold = pixel_min)
  }
  wl_stats <- wl_stats |> filter(n_used >= pixel_min)   # mirrors Hypernets_matchups' own minimum-valid-pixel rule
  if(nrow(wl_stats) == 0) return(list(written = FALSE, removal_log = removal_log))

  hyp_stats <- df_matchup |> distinct(wavelength, Hyp) |> filter(wavelength %in% wl_stats$wavelength)

  # Single-scalar CV proxy (variability_centered is one value per file, not per wavelength) --
  # exact upstream Hypernets_matchups formula unknown. This reuses db_load_spectra()'s cv_pct
  # logic (100 * sd/|mean|), pooled across wavebands with >= 2 contributing pixels. Uses the
  # median rather than the mean across wavebands. A handful of near-zero-mean bands (e.g.
  # S3A's 665/673/681 nm) otherwise produce cv_pct in the thousands of percent and dominate a
  # simple mean, even though most bands sit in a plausible range (confirmed empirically).
  cv_pct_by_wl <- wl_stats |> filter(n_used >= 2, mean_val != 0) |>
    mutate(cv_pct = 100 * abs(sd_val / mean_val))
  # Fall back to 0 (not NA) when every waveband has < 2 contributing pixels. With a single
  # pixel there is no measurable spatial spread, so "no pixel disagreement detected" is the
  # correct reading, and it keeps variability_centered numeric (sat_var_check() errors on an
  # all-NA column, since it expects exactly one non-NA value per file)
  cv_scalar <- if(nrow(cv_pct_by_wl) > 0) median(cv_pct_by_wl$cv_pct, na.rm = TRUE) else 0

  sat_lat <- mean(df_matchup$pixel_lat, na.rm = TRUE)
  sat_lon <- mean(df_matchup$pixel_lon, na.rm = TRUE)
  sat_pixel_pos <- paste(sort(unique(df_matchup$pixel_pos)), collapse = ";")

  out_name <- db_matchup_filename(sensor_Y, df_matchup$dateTime_Hyp[1],
                                   df_matchup[[paste0("dateTime_", sensor_Y)]][1], tag)
  hyp_day <- format(df_matchup$dateTime_Hyp[1], "%Y%m%d", tz = "Europe/Paris")
  hyp_time <- format(df_matchup$dateTime_Hyp[1], "%H%M%S", tz = "Europe/Paris")
  sat_day <- format(df_matchup[[paste0("dateTime_", sensor_Y)]][1], "%Y%m%d", tz = "Europe/Paris")
  sat_time <- format(df_matchup[[paste0("dateTime_", sensor_Y)]][1], "%H%M%S", tz = "Europe/Paris")

  wl_chr <- as.character(wl_stats$wavelength)
  hyp_row <- c(list(sensor = "Hyp 1", data_type = data_type_label, day = hyp_day, time = hyp_time,
                     latitude = round(df_matchup$lat_Hyp[1], 6), longitude = round(df_matchup$lon_Hyp[1], 6),
                     radiometer_id = 1, pixel_pos = "", type = "in-situ"),
               setNames(as.list(hyp_stats$Hyp[match(wl_stats$wavelength, hyp_stats$wavelength)]), wl_chr),
               list(variability_centered = NA_real_))
  sat_row <- c(list(sensor = paste0(sensor_Y, " 1"), data_type = data_type_label, day = sat_day, time = sat_time,
                     latitude = round(sat_lat, 6), longitude = round(sat_lon, 6),
                     radiometer_id = sat_radiometer_id, pixel_pos = sat_pixel_pos, type = "satellite"),
               setNames(as.list(wl_stats$mean_val), wl_chr),
               list(variability_centered = round(cv_scalar, 6)))

  df_new <- bind_rows(as_tibble(hyp_row), as_tibble(sat_row))
  colnames(df_new)[1] <- ""

  write_delim(df_new, file.path(out_dir, out_name), delim = ";", quote = "needed", na = "")

  list(written = TRUE, removal_log = removal_log)
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

# Restrict to the inner 3x3 grid (col/row in -1..1) -- matches the box size Hypernets_matchups
# actually exports to .csv (see pixel_filter_ne_inner3x3()'s comment for the raw-5x5-vs-exported-3x3
# discrepancy) -- but apply no further spatial restriction on top of that, unlike
# pixel_filter_ne_inner3x3() (NE quadrant) or pixel_filter_clean_water() (hand-drawn polygon).
# Isolates the effect of the shared per-pixel QC gates (RHOW ceiling, negative-value, distance,
# minimum-valid-pixel-count) from any spatial subsetting. Renamed 2026-09-04 from
# pixel_filter_none() (which applied literally no restriction, so THFR_pixel was drawing from the
# database's raw 5x5/25-pixel grid rather than Hypernets_matchups' real 3x3/9-pixel box).
pixel_filter_inner3x3 <- function(df_pixel){
  df_pixel |>
    filter(!grepl("2", pixel_pos))
}

# Looks up a radiometer's own id in the .db, e.g. for stamping the correct radiometer_id into a
# from-scratch-written matchup row. Small helper factoring out a query already inlined identically
# in db_matchup_long()/db_matchup_pixels().
db_radiometer_id <- function(db_path, sensor_Y){
  con <- dbConnect(RSQLite::SQLite(), path.expand(db_path))
  on.exit(dbDisconnect(con), add = TRUE)
  dbGetQuery(con, sprintf("SELECT id FROM radiometer WHERE name = '%s'", sensor_Y))$id
}

# Captures the rows one QC/spatial gate removes from df, at a given grouping grain, before they're
# actually filtered out -- used by db_export_matchups_site() below to build the per-site
# pixel-removal audit CSV instead of letting dropped pixels vanish with no record. group_vars
# controls granularity: "wavelength" in group_vars means the reason genuinely depends on which
# waveband (kept separate); "pixel_pos" in group_vars means the reason depends on which specific
# pixel (kept separate) -- when absent, pixel_pos is written as the "ALL" sentinel because the
# gate's value/reason doesn't vary by pixel (e.g. Hyp-side gates, matchup-level time/distance).
# value_col is the name of the column in df holding the value that triggered the gate (NULL for
# gates with no single numeric trigger, e.g. a spatial-filter exclusion).
# Structural key (matchup_id, db_source, dateTime_Hyp, dateTime_sat) rather than a per-matchup
# CSV filename: file_name only means anything while a physical per-matchup CSV exists, whereas
# these come straight from measure_info/matchups and stay meaningful even once ingestion no
# longer produces a per-matchup file. db_matchup_filename() can derive a filename on demand from
# these fields for joining against today's still-file-based output/matchup_stats_RHOW_*.csv.
log_gate_drop <- function(df, keep_lgl, gate, group_vars, value_col = NULL, threshold = NA_real_,
                           sensor_Z, sensor_Y, site_name){
  dropped <- df[!keep_lgl, , drop = FALSE]
  if(nrow(dropped) == 0) return(tibble())
  has_wavelength <- "wavelength" %in% group_vars
  has_pixel <- "pixel_pos" %in% group_vars
  sat_dt_col <- paste0("dateTime_", sensor_Y)
  distinct_cols <- c("matchup_id", "match_date", "db_source", "dateTime_Hyp", sat_dt_col, group_vars)
  if(!is.null(value_col)) distinct_cols <- c(distinct_cols, value_col)
  dropped |>
    distinct(across(all_of(distinct_cols))) |>
    transmute(site_name = site_name, sensor_Z = sensor_Z, sensor_Y = sensor_Y,
              matchup_id = matchup_id, match_date = match_date, db_source = db_source,
              dateTime_Hyp = dateTime_Hyp, dateTime_sat = .data[[sat_dt_col]],
              gate = gate,
              pixel_pos = if(has_pixel) pixel_pos else "ALL",
              wavelength = if(has_wavelength) wavelength else NA_real_,
              value = if(!is.null(value_col)) .data[[value_col]] else NA_real_,
              threshold = threshold)
}

# Regenerates THFR's raw per-matchup RHOW CSV files for every sensor_Y in a sensor family,
# restricted to whichever pixels pixel_filter_fn keeps, and written under a new site folder
# named site_name (isolated from the real THFR/MAFR data). See the section comment above for full
# context. Matchups with zero surviving pixels for a given sensor are skipped entirely (no file
# written) -- for the coarser sensors (AQUA/VIIRS/PACE) this is expected to remove most matchups
# under the NE+3x3 filter, since the box is centered well outside the NE quadrant for those (see
# meta/pixel_explore_output/summary.md); a similarly small survival rate is expected under the
# clean-water-polygon filter, which is deliberately small and close to the station.
# Returns a tibble accumulating every pixel/waveband dropped along the way (QC gates, the site's
# own spatial filter, and write_matchup_csv_db()'s per-wavelength minimum-pixel-count drop), in
# the schema produced by log_gate_drop() -- callers (the db_export_matchups_ne/_poly/_pixel()
# wrappers below) combine this across all 4 sensor families into one CSV per site.
# sensor_Z = "OLCI"; site_name = "THFR_NE"; pixel_filter_fn = pixel_filter_ne_inner3x3; file_tag = "NE"; qc_site_name = "THFR"
# apply_pixel_qc = FALSE skips the four pixel-level plausibility gates below (RHOW ceiling and
# negative-value on both Hyp and the satellite pixel, plus the per-pixel distance gate), keeping
# only the matchup-level diff_time gate, the spatial filter (box shape, not QC), and
# write_matchup_csv_db()'s own minimum-valid-pixel-count gate. Used by the "_raw" sites
# (db_export_matchups_thfr_raw()/db_export_matchups_mafr_raw() below) to reconstruct MAFR/THFR as a
# deliberately unfiltered 3x3-box average -- as close as this repo can get to reproducing whatever
# Hypernets_matchups' own opaque aggregation currently does, so validate_derived_site() has a fair
# like-for-like comparison to check before anything reading the real tool's CSVs is retargeted.
db_export_matchups_site <- function(sensor_Z, site_name, pixel_filter_fn, file_tag, qc_site_name,
                                     db_path = "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db",
                                     apply_pixel_qc = TRUE){
  removal_log_list <- list()

  sensor_Y_list <- sensor_grid(sensor_Z)$sensor_Y |> unique()

  for(sensor_Y in sensor_Y_list){
    message("db_export_matchups_site(", site_name, "): ", sensor_Y)
    df_pixel <- tryCatch(db_matchup_pixels(db_path, sensor_Y),
                          warning = function(w){ message(conditionMessage(w)); tibble() },
                          error = function(e){ message(conditionMessage(e)); tibble() })
    if(nrow(df_pixel) == 0) next

    # db_matchup_pixels() returns every candidate pairing in the db's raw `matchups` table,
    # time/distance QC is deliberately left to the caller (see its own docstring). The exported
    # .csv files on disk only exist for matchups that already passed this gate, so it must be
    # applied here before looking for a matching original file. diff_time_min is already
    # matchup-level (constant per matchup_id), dist_km as returned by db_matchup_pixels() is
    # per-PIXEL (station-to-pixel, used for the bearing/quadrant diagnostics), so the
    # matchup-level distance (station-to-satellite-reference-position, matching how
    # db_matchup_long()/process_sensor() gate distance) is recomputed here instead. qc_site_name
    # states which real site's QC policy applies -- "THFR" for the three THFR-derived sites (all
    # reinterpretations of the same underlying THFR matchups), "MAFR" for a MAFR-sourced site.
    # The error= handler above copes with mafr_2025.db's missing SNPP column (Bug 11,
    # manuscript/upstream-data-bugs.md) the same way code/5_figures.R's Figure 1 loop already does.
    #
    # Four further per-pixel QC gates in total, before write_matchup_csv_db() aggregates the pixel
    # box into one matchup value (and long before the separate day-level daily_closest_matchup()
    # step further downstream in the main CSV pipeline). Three are applied here, directly on
    # df_pixel: a site-specific RHOW plausibility ceiling (site_rhow_limit()) applied to both Hyp
    # and the satellite pixel value, a per-pixel distance gate (TEMPORARY flat 5 km, 2026-09-04 --
    # the resolution-based version, 2x/3x sensor_resolution_km(), turned out to eliminate
    # essentially all MAFR matchups when applied to the real-data pipeline; reverted here too for
    # consistency until it's properly revisited next week), and a per-pixel-per-waveband
    # negative-value gate (RHOW < 0 is optically impossible; satellite pixels commonly go negative
    # after an over-aggressive atmospheric/glint correction, most often in the blue bands) applied
    # to both Hyp and the satellite pixel value. The fourth -- a site-specific minimum-valid-pixel-
    # count per wavelength (site_pixel_min()) -- is applied inside write_matchup_csv_db(), since it
    # has to be evaluated per-wavelength on the pixel_filter_fn()-restricted candidate pool, not on
    # df_pixel directly.
    #
    # Applied as a sequential waterfall (rather than one combined filter()) so each dropped row can
    # be attributed to exactly the one gate that removed it, for the pixel-removal audit log --
    # the final surviving set is identical to a combined AND-filter either way.
    sat_lon_col <- paste0("lon_", sensor_Y); sat_lat_col <- paste0("lat_", sensor_Y)
    rhow_limit <- site_rhow_limit(qc_site_name)
    pixel_dist_limit <- 5
    df_pixel <- df_pixel |>
      mutate(db_source = basename(db_path),
             dist_km_matchup = distHaversine(cbind(lon_Hyp, lat_Hyp), cbind(.data[[sat_lon_col]], .data[[sat_lat_col]])) / 1000)

    keep_time <- df_pixel$diff_time_min <= site_diff_time_limit(qc_site_name)
    removal_log_list[[paste0(sensor_Y, "_diff_time")]] <- log_gate_drop(
      df_pixel, keep_time, "diff_time", "matchup_id", "diff_time_min", site_diff_time_limit(qc_site_name),
      sensor_Z, sensor_Y, site_name)
    df_pixel <- df_pixel[keep_time, , drop = FALSE]

    if(apply_pixel_qc){
      keep_hyp_rhow <- df_pixel$Hyp <= rhow_limit
      removal_log_list[[paste0(sensor_Y, "_hyp_rhow")]] <- log_gate_drop(
        df_pixel, keep_hyp_rhow, "hyp_rhow_ceiling", c("matchup_id", "wavelength"), "Hyp", rhow_limit,
        sensor_Z, sensor_Y, site_name)
      df_pixel <- df_pixel[keep_hyp_rhow, , drop = FALSE]

      keep_sat_rhow <- df_pixel[[sensor_Y]] <= rhow_limit
      removal_log_list[[paste0(sensor_Y, "_sat_rhow")]] <- log_gate_drop(
        df_pixel, keep_sat_rhow, "sat_rhow_ceiling", c("matchup_id", "wavelength", "pixel_pos"), sensor_Y, rhow_limit,
        sensor_Z, sensor_Y, site_name)
      df_pixel <- df_pixel[keep_sat_rhow, , drop = FALSE]

      keep_hyp_neg <- df_pixel$Hyp >= 0
      removal_log_list[[paste0(sensor_Y, "_hyp_neg")]] <- log_gate_drop(
        df_pixel, keep_hyp_neg, "hyp_negative", c("matchup_id", "wavelength"), "Hyp", 0,
        sensor_Z, sensor_Y, site_name)
      df_pixel <- df_pixel[keep_hyp_neg, , drop = FALSE]

      keep_sat_neg <- df_pixel[[sensor_Y]] >= 0
      removal_log_list[[paste0(sensor_Y, "_sat_neg")]] <- log_gate_drop(
        df_pixel, keep_sat_neg, "sat_negative", c("matchup_id", "wavelength", "pixel_pos"), sensor_Y, 0,
        sensor_Z, sensor_Y, site_name)
      df_pixel <- df_pixel[keep_sat_neg, , drop = FALSE]

      keep_dist <- df_pixel$dist_km_matchup <= pixel_dist_limit
      removal_log_list[[paste0(sensor_Y, "_dist")]] <- log_gate_drop(
        df_pixel, keep_dist, "dist_matchup", "matchup_id", "dist_km_matchup", pixel_dist_limit,
        sensor_Z, sensor_Y, site_name)
      df_pixel <- df_pixel[keep_dist, , drop = FALSE]
    }

    if(nrow(df_pixel) == 0){
      message("  No matchups for ", sensor_Y, " pass the time/distance", if(apply_pixel_qc) "/RHOW/pixel-distance", " QC gate -- skipping")
      next
    }

    df_pixel <- df_pixel |> mutate(row_uid = row_number())
    df_filt <- pixel_filter_fn(df_pixel)
    spatial_gate <- switch(file_tag, NE = "spatial_ne_quadrant", poly = "spatial_clean_water",
                           pixel = "spatial_inner3x3", paste0("spatial_", file_tag))
    keep_spatial <- df_pixel$row_uid %in% df_filt$row_uid
    removal_log_list[[paste0(sensor_Y, "_spatial")]] <- log_gate_drop(
      df_pixel, keep_spatial, spatial_gate, "pixel_pos", value_col = NULL, threshold = NA_real_,
      sensor_Z = sensor_Z, sensor_Y = sensor_Y, site_name = site_name)

    if(nrow(df_filt) == 0){
      message("  No pixels survive the ", site_name, " filter for ", sensor_Y, " -- skipping")
      next
    }

    pixel_min <- site_pixel_min(site_name, sensor_Y)
    sat_radiometer_id <- db_radiometer_id(db_path, sensor_Y)

    out_dir <- file_path_build(site_name, sensor_Y)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    for(mid in unique(df_filt$matchup_id)){
      csv_result <- write_matchup_csv_db(filter(df_filt, matchup_id == mid), sensor_Y, sat_radiometer_id,
                                          out_dir, pixel_min, tag = file_tag,
                                          sensor_Z = sensor_Z, site_name = site_name)
      removal_log_list[[paste0(sensor_Y, "_pixel_min_", mid)]] <- csv_result$removal_log
    }
  }

  bind_rows(removal_log_list)
}

# Thin, site-specific wrappers -- add a new derived site by writing one of these (plus, if the
# inclusion rule is new, a pixel_filter_*() function above). Each loops internally over every
# sensor family (rather than being called once per sensor_Z externally) so the pixel-removal log
# can be accumulated in memory across all 4 sensor families and written once per site -- mirroring
# how meta/satellite_outliers.csv is already built (code/2_outliers.R), rather than read-modify-
# appending to a file across separate calls.
db_export_matchups_ne <- function(sensor_Z_vec = c("MODIS", "VIIRS", "OLCI", "OCI"),
                                   db_path = "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db"){
  removal_log <- purrr::map_dfr(sensor_Z_vec, db_export_matchups_site,
                                 site_name = "THFR_NE", pixel_filter_fn = pixel_filter_ne_inner3x3,
                                 file_tag = "NE", qc_site_name = "THFR", db_path = db_path)
  write_csv(removal_log, "meta/THFR_NE_pixel_removals.csv")
  invisible(removal_log)
}
db_export_matchups_poly <- function(sensor_Z_vec = c("MODIS", "VIIRS", "OLCI", "OCI"),
                                     db_path = "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db"){
  removal_log <- purrr::map_dfr(sensor_Z_vec, db_export_matchups_site,
                                 site_name = "THFR_poly", pixel_filter_fn = pixel_filter_clean_water,
                                 file_tag = "poly", qc_site_name = "THFR", db_path = db_path)
  write_csv(removal_log, "meta/THFR_poly_pixel_removals.csv")
  invisible(removal_log)
}
db_export_matchups_pixel <- function(sensor_Z_vec = c("MODIS", "VIIRS", "OLCI", "OCI"),
                                      db_path = "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db"){
  removal_log <- purrr::map_dfr(sensor_Z_vec, db_export_matchups_site,
                                 site_name = "THFR_pixel", pixel_filter_fn = pixel_filter_inner3x3,
                                 file_tag = "pixel", qc_site_name = "THFR", db_path = db_path)
  write_csv(removal_log, "meta/THFR_pixel_pixel_removals.csv")
  invisible(removal_log)
}

# Loops db_export_matchups_site() (unchanged) over every (db_path, sensor_Z) combination for a
# site drawing from more than one .db file -- needed for MAFR, whose two files (mafr_2024.db/
# mafr_2025.db) have independent auto-increment matchups.id values that collide across the two.
# db_export_matchups_site() already stamps db_source = basename(db_path) onto every removal-log
# row itself, so this just combines and writes once per site. Output CSV filenames themselves
# never collide across db_paths since they're built from timestamps, not matchup_id.
db_export_matchups_multi <- function(sensor_Z_vec, db_paths, site_name, pixel_filter_fn,
                                      file_tag, qc_site_name, apply_pixel_qc = TRUE){
  removal_log <- purrr::map_dfr(db_paths, function(db_path){
    purrr::map_dfr(sensor_Z_vec, db_export_matchups_site,
                    site_name = site_name, pixel_filter_fn = pixel_filter_fn,
                    file_tag = file_tag, qc_site_name = qc_site_name, db_path = db_path,
                    apply_pixel_qc = apply_pixel_qc)
  })
  write_csv(removal_log, paste0("meta/", site_name, "_pixel_removals.csv"))
  invisible(removal_log)
}

# Regenerates MAFR's own matchup data directly from its .db files, mirroring THFR_pixel exactly:
# same pixel_filter_inner3x3 (3x3 box, no extra spatial restriction -- what real Hypernets_matchups
# itself exports), same "pixel" tag, same per-matchup-file convention under file_path_build()'s
# existing FR/ tree -- just sourced from MAFR's two .db files with MAFR's own QC policy
# (site_rhow_limit("MAFR")/site_diff_time_limit("MAFR")) instead of THFR's. Wired into
# available_sites()'s candidate list and called from code/1_matchups_single.R.
db_export_matchups_mafr_pixel <- function(sensor_Z_vec = c("MODIS", "VIIRS", "OLCI", "OCI"),
                                           db_paths = c("~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/mafr_2024.db",
                                                        "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/mafr_2025.db")){
  db_export_matchups_multi(sensor_Z_vec, db_paths, site_name = "MAFR_pixel",
                            pixel_filter_fn = pixel_filter_inner3x3, file_tag = "pixel",
                            qc_site_name = "MAFR")
}

# Deliberately unfiltered reconstruction of THFR/MAFR themselves (apply_pixel_qc = FALSE -- see its
# docstring on db_export_matchups_site() above): same 3x3 box as the real tool exports
# (pixel_filter_inner3x3), same per-matchup-file convention, but none of the RHOW-ceiling/negative-
# value/pixel-distance gates, so this is as close to Hypernets_matchups' own (opaque) box-averaging
# as this repo can reconstruct. Written to "THFR_raw"/"MAFR_raw" -- a new site, not a replacement
# for THFR/MAFR -- specifically so validate_derived_site() can check it against the real tool
# output before anything is retargeted onto it.
db_export_matchups_thfr_raw <- function(sensor_Z_vec = c("MODIS", "VIIRS", "OLCI", "OCI"),
                                         db_path = "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db"){
  removal_log <- purrr::map_dfr(sensor_Z_vec, db_export_matchups_site,
                                 site_name = "THFR_raw", pixel_filter_fn = pixel_filter_inner3x3,
                                 file_tag = "raw", qc_site_name = "THFR", db_path = db_path,
                                 apply_pixel_qc = FALSE)
  write_csv(removal_log, "meta/THFR_raw_pixel_removals.csv")
  invisible(removal_log)
}
db_export_matchups_mafr_raw <- function(sensor_Z_vec = c("MODIS", "VIIRS", "OLCI", "OCI"),
                                         db_paths = c("~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/mafr_2024.db",
                                                      "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/mafr_2025.db")){
  db_export_matchups_multi(sensor_Z_vec, db_paths, site_name = "MAFR_raw",
                            pixel_filter_fn = pixel_filter_inner3x3, file_tag = "raw",
                            qc_site_name = "MAFR", apply_pixel_qc = FALSE)
}

# Validation diagnostic (non-blocking) -------------------------------------
# Since THFR_NE/THFR_poly/THFR_pixel/THFR_raw are reinterpretations of the same underlying THFR
# matchups (and MAFR_pixel/MAFR_raw of the same underlying MAFR matchups), a real Hypernets_matchups
# file exists for many of the same matchup timestamps -- giving a ground truth to numerically check
# write_matchup_csv_db()'s from-scratch reconstruction against, since the real tool's own
# box-averaging/variability_centered formula is undocumented anywhere.
# Not run as part of the main pipeline (see the "Validation Hypernets_matchups vs direct .db"
# section of code/3_sensitivity.R, which is where these calls actually live) -- purely diagnostic,
# flags mismatches for human review rather than asserting/blocking.
validate_derived_site <- function(site_name, tag, db_path = "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db",
                                   real_site_name = "THFR"){
  results <- purrr::map_dfr(db_satellite_names, function(sensor_Y){
    derived_dir <- file_path_build(site_name, sensor_Y)
    if(!dir.exists(derived_dir)) return(tibble())
    derived_files <- list.files(derived_dir, pattern = paste0("_RHOW_", tag, "\\.csv$"), full.names = TRUE)
    if(length(derived_files) == 0) return(tibble())

    real_dir <- file_path_build(real_site_name, sensor_Y)

    purrr::map_dfr(derived_files, function(f){
      fname <- basename(f)
      m <- regmatches(fname, regexec(
        paste0("^", sensor_Y, "_(\\d{8}T\\d{6})_vs_HYPERNETS_(\\d{8}T\\d{6})_RHOW_", tag, "\\.csv$"), fname))[[1]]
      if(length(m) != 3) return(tibble())
      sat_ts <- m[2]; hyp_ts <- m[3]

      real_prefix <- paste0(sensor_Y, "_", sat_ts, "_vs_HYPERNETS_", hyp_ts)
      real_files <- list.files(real_dir, pattern = paste0("^", real_prefix, ".*\\.csv$"), full.names = TRUE)
      if(length(real_files) == 0){
        message("  validate_derived_site(): no real ", real_site_name, " file for ", real_prefix, " -- skipping")
        return(tibble())
      }

      df_derived <- tryCatch(load_matchup_mean(f), error = function(e) NULL)
      df_real <- tryCatch(load_matchup_mean(real_files[1]), error = function(e) NULL)
      if(is.null(df_derived) || is.null(df_real)) return(tibble())

      to_wide <- function(df, suffix){
        df |>
          # day/time/latitude/longitude differ between the Hyp and satellite rows (each sensor's
          # own timestamp/position) -- must be dropped before pivot_wider() the same way
          # load_matchup_long() already does, otherwise they act as extra id columns and pivot_wider
          # can't merge the two rows into one per wavelength
          dplyr::select(-day, -time, -latitude, -longitude) |>
          pivot_longer(cols = matches("^[0-9]+$"), names_to = "wavelength", values_to = "value") |>
          mutate(wavelength = as.numeric(wavelength)) |>
          na.omit() |>
          pivot_wider(names_from = sensor, values_from = value) |>
          rename_with(~paste0(.x, suffix), -wavelength)
      }

      joined <- tryCatch(
        inner_join(to_wide(df_derived, "_derived"), to_wide(df_real, "_real"), by = "wavelength"),
        error = function(e) tibble())
      if(nrow(joined) == 0) return(tibble())

      sat_derived_col <- paste0(sensor_Y, "_derived"); sat_real_col <- paste0(sensor_Y, "_real")
      if(!all(c(sat_derived_col, sat_real_col) %in% colnames(joined))) return(tibble())

      joined |>
        transmute(site_name = site_name, sensor_Y = sensor_Y, file_name = fname, wavelength,
                  Hyp_real, Hyp_derived,
                  Hyp_diff = Hyp_derived - Hyp_real,
                  Hyp_pct_diff = 100 * Hyp_diff / Hyp_real,
                  Sat_real = .data[[sat_real_col]], Sat_derived = .data[[sat_derived_col]],
                  Sat_diff = Sat_derived - Sat_real,
                  Sat_pct_diff = 100 * Sat_diff / Sat_real)
    })
  })

  if(nrow(results) == 0){
    message("validate_derived_site(", site_name, "): no overlapping real/derived matchups found to compare")
    return(invisible(results))
  }

  write_csv(results, paste0("meta/validate_", site_name, ".csv"))

  summary_tbl <- results |>
    summarise(median_abs_Hyp_pct = median(abs(Hyp_pct_diff), na.rm = TRUE),
              iqr_abs_Hyp_pct = IQR(abs(Hyp_pct_diff), na.rm = TRUE),
              median_abs_Sat_pct = median(abs(Sat_pct_diff), na.rm = TRUE),
              iqr_abs_Sat_pct = IQR(abs(Sat_pct_diff), na.rm = TRUE),
              n_flagged_over_5pct = sum(abs(Hyp_pct_diff) > 5 | abs(Sat_pct_diff) > 5, na.rm = TRUE),
              n = n(),
              .by = sensor_Y)
  print(summary_tbl)

  invisible(results)
}


# Plotting functions ------------------------------------------------------

# Plot data based on wavelength group
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

    # Unique satellite days for this site only
    unique_days_i <- df_prep_i |>
      dplyr::select(file_name) |>
      mutate(date = sapply(str_split(file_name, "_"), "[[", 2)) |>
      mutate(date = sapply(str_split(date, "T"), "[[", 1)) |>
      distinct(date) |>
      mutate(date = as.Date(date, format = "%Y%m%d"))

    df_stats_i <- base_stats(x_vec, y_vec) |>
      mutate(site_name = site_i,
             max_axis = max(x_vec, y_vec, na.rm = TRUE),
             label = paste0("n: ", n, " (", nrow(unique_days_i), ")",
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
      mutate(wavelength = pace_waveband_bucket(wavelength), .after = "wavelength")
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

# Load QC-passed, outlier-screened, long-format HYPERNETS-vs-satellite matchups for one
# sensor_Y platform. Shared by global_scatterplot() and global_scatterplot_waveband().
# Returns columns: site_name, file_name, wavelength, Hyp, <sensor_Y column>
# sites: optional character vector restricting which site(s) to load (must be valid
# available_sites() candidates); NULL (default) preserves the original behaviour of loading
# every site found on disk for this sensor.
load_global_matchup_data <- function(sensor_Y, sites = NULL){

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
  # once its data folder exists on disk, no code change needed here. If `sites` is supplied,
  # it overrides this auto-discovery (still validated against what's actually on disk).
  print("Loading matchups")
  if(sensor_Y == "S3"){
    site_list <- unique(unlist(lapply(c("S3A", "S3B"), available_sites)))
    if(!is.null(sites)) site_list <- intersect(site_list, sites)
    ply_folders <- expand.grid(site_name = site_list, sat_name = c("S3A", "S3B"))
    # match_base_1 <- bind_rows(load_matchups_folder("S3A", long = TRUE),
    #                           load_matchups_folder(e_name, "S3B", long = TRUE))
  } else {
    # match_base_1 <- load_matchups_folder(site_name, sensor_Y, long = TRUE)
    site_list <- available_sites(sensor_Y)
    if(!is.null(sites)) site_list <- intersect(site_list, sites)
    ply_folders <- expand.grid(site_name = site_list, sat_name = sensor_Y)
  }

  # Load all folders
  match_base <- purrr::pmap_dfr(ply_folders, load_matchups_folder, long = TRUE) |>
    right_join(match_base_details, by = join_by(site_name, file_name))

  # Filter out outliers
  match_filter <- match_base[!match_base$file_name %in% outliers_all$file_name,]

  # Remove any erroneuosly high values before plotting
  match_filter <- match_filter |> filter(Hyp <= 1)

  match_filter
}

# Takes variable and Y sensor as input to automagically create global scatterplot triptych
# cut_legend = "cut" strips the legend so several panels can be stacked under one shared legend
# sensor_Y = "S3B"; cut_legend = "cut"
# sensor_Y = "AQUA"; cut_legend = "cut"
# sensor_Y = "S3"; cut_legend = "no"
# sensor_Y = "PACE"; cut_legend = "no"
global_scatterplot <- function(sensor_Y, cut_legend = "no", sites = NULL){

  # Load QC-passed, outlier-screened, long-format matchups
  print("Loading matchups")
  match_filter <- load_global_matchup_data(sensor_Y, sites = sites)

  # Create the figure
  print("Creating figure, saving, and exiting")
  match_fig <- plot_global_nm(match_filter, sensor_Y)

  # Save the individual figure -- width scales with the number of site columns actually plotted
  # (plot_global_nm() facets ~site_name, nrow = 1), rather than a fixed value tuned for the original
  # 8-candidate-site case, so a restricted `sites` list (e.g. the 2-site manuscript_sites filter in
  # code/5_figures.R) doesn't leave panels stretched across a canvas sized for far more columns.
  n_sites <- length(unique(match_filter$site_name))
  ggsave(paste0("figures/global_scatter_RHOW_",sensor_Y,".png"), match_fig, width = max(2 * n_sites, 6), height = 5)

  # Remove the legend when this panel will be stacked beneath another with a shared legend
  # Then return it (invisibly) for reuse by global_scatterplot_stack()
  if(cut_legend == "cut") match_fig <- match_fig + theme(legend.position = "none")
  invisible(match_fig)
}

# Stack the per-sensor global scatterplots for a sensor family into one composite figure
# NB: Requires that the matchup and global stats CSVs for sensor_Z already exist (see process_sensor())
# sensor_Z = "MODIS"
global_scatterplot_stack <- function(sensor_Z, sites = NULL){

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
    fig_list[[i]] <- global_scatterplot(sensor_Y_list[i], cut_legend = ifelse(i < sensor_count, "cut", "no"), sites = sites)
  }

  # Give the legend-bearing (last) panel extra relative height to fit the legend
  panel_heights <- c(rep(1, sensor_count - 1), 1.15)

  # Stack panels vertically and exit
  fig_stack <- ggpubr::ggarrange(plotlist = fig_list, ncol = 1, nrow = sensor_count, heights = panel_heights) +
    ggpubr::bgcolor("white") + ggpubr::border("white", size = 2)
  # Width scales with the number of site columns actually plotted (see global_scatterplot()'s
  # equivalent note) rather than a fixed value tuned for the original 8-candidate-site case.
  n_sites <- if(!is.null(sites)) length(sites) else length(available_sites(sensor_Y_list[1]))
  ggsave(paste0("figures/global_scatter_RHOW_",sensor_Z,".png"), fig_stack, width = max(2.5 * n_sites, 8), height = 5 * sensor_count)
}

# Single-platform scatterplot faceted by (site, waveband) -- complements global_scatterplot_stack()
# (which facets by site and colours by waveband) by making a single systematically bad waveband
# immediately visible as its own panel, rather than scattered across a wavelength-coloured legend.
# One figure per platform (sensor_Y), not per sensor family (2026-09-10 rewrite): with only the
# 2-site manuscript_sites filter typically in play, combining a whole family's platforms into one
# figure (as the previous cross-platform-band-aligned version did) is no longer necessary, and a
# true per-(site, waveband) grid -- each panel scaled independently to its own data, both axes --
# is clearer built one platform at a time. Colour is mapped to waveband, matching plot_global_nm()'s
# convention in the "_RHOW_<sensor>" figures; site is not separately encoded by shape/colour since
# each site already has its own row of panels, identified by its facet label.
# sensor_Y = "SNPP"
global_scatterplot_waveband <- function(sensor_Y, sites = NULL){

  # Load QC-passed, outlier-screened, long-format matchups for this one platform
  match_filter <- load_global_matchup_data(sensor_Y, sites = sites) |>
    filter(wavelength %in% W_nm_out(sensor_Y)) |>
    rename(Sat = all_of(sensor_Y))

  # PACE/OCI is continuous (350-1150 nm) -- bucket into the same broad bands used for its
  # colour legend elsewhere rather than one facet per nm
  if(sensor_Y == "PACE"){
    match_filter <- match_filter |> mutate(waveband_label = pace_waveband_bucket(wavelength))
    waveband_levels <- names(colour_nm_func("PACE"))
  } else {
    match_filter <- match_filter |> mutate(waveband_label = as.character(wavelength))
    waveband_levels <- as.character(sort(unique(as.numeric(match_filter$waveband_label))))
  }
  match_filter <- match_filter |> mutate(waveband_label = factor(waveband_label, levels = waveband_levels))

  # When `sites` is supplied, use it directly (in that order) for the row layout below; otherwise
  # fall back to whatever available_sites() finds on disk for this platform.
  site_levels <- if(!is.null(sites)) sites else available_sites(sensor_Y)
  match_filter <- match_filter |> mutate(site_name = factor(site_name, levels = site_levels))
  n_sites <- length(site_levels)

  # Distinct satellite overpass days per (site, waveband) -- same date-extraction convention as
  # plot_global_nm()'s parenthetical day count, computed before base_stats()'s NA/positive-value
  # filtering (i.e. counts distinct days contributing a matchup, not distinct valid point pairs)
  df_days <- match_filter |>
    dplyr::select(site_name, waveband_label, file_name) |>
    mutate(date = sapply(str_split(file_name, "_"), "[[", 2),
           date = sapply(str_split(date, "T"), "[[", 1),
           date = as.Date(date, format = "%Y%m%d")) |>
    distinct(site_name, waveband_label, date) |>
    count(site_name, waveband_label, name = "n_days")

  # Per-(site, waveband) stats and per-panel max (for 1:1 axis limits below), formatted the same
  # multi-line way as plot_global_nm()'s per-site label (one stat per row: n, S, beta, epsilon).
  # tidyr::complete() fills in every (site_name, waveband_label) combination, including ones with
  # no surviving match-ups at all (e.g. THFR_pixel's PACE data doesn't reach every waveband a site
  # with fuller coverage does) -- without this, facet_wrap() below only draws panels for combinations
  # that actually have data, which silently breaks the row = site / column = waveband alignment
  # whenever one site is missing a waveband the other has (confirmed 2026-09-10: PACE's THFR_pixel
  # panels were shifted left relative to MAFR_pixel's without this fix).
  df_stats <- match_filter |>
    group_by(site_name, waveband_label) |>
    group_modify(~ base_stats(.x$Hyp, .x$Sat) |>
                   mutate(max_axis = suppressWarnings(max(c(.x$Hyp, .x$Sat), na.rm = TRUE)))) |>
    ungroup() |>
    tidyr::complete(site_name, waveband_label) |>
    left_join(df_days, by = c("site_name", "waveband_label")) |>
    arrange(site_name, waveband_label) |>
    mutate(label = if_else(is.na(n) | is.na(max_axis) | !is.finite(max_axis),
                           "no match-ups",
                           paste0("n: ", n, " (", n_days, ")",
                                  "\nS: ", sprintf("%.2f", Slope_II),
                                  "\nβ: ", sprintf("%.1f", Bias_50), "%",
                                  "\nε: ", sprintf("%.1f", Error_50), "%")),
           # Fallback text-placement y for empty panels (where max_axis is NA/non-finite and the
           # panel's own y-scale is left to auto-compute, limits = c(0, NA)) -- doesn't affect the
           # axis itself, just keeps the "no match-ups" label from being silently dropped.
           label_y = if_else(is.na(max_axis) | !is.finite(max_axis), 1, max_axis))

  # True per-cell (not just per-row or per-column) independent 1:1 axis limits, each panel scaled
  # to its own data -- facet_grid's native free scales only vary per row/column, not per cell, so
  # this uses the same facet_wrap() + ggh4x::facetted_pos_scales() technique as plot_global_nm()
  # (see that function's own note on why coord_fixed()/facet_grid can't achieve this), with an
  # explicit nrow = n_sites so the wrap still lays out as site rows x waveband columns: facet_wrap()
  # fills panels in the order of the *sorted* combinations of its facetting variables (site_name
  # then waveband_label here, matching df_stats' arrange() above), so the first n_wavebands panels
  # are all of site 1's wavebands, the next n_wavebands are site 2's, and so on -- one full row per
  # site. A missing/non-finite max_axis (an empty cell after complete()) falls back to limits =
  # c(0, NA), i.e. let ggplot auto-scale that one (empty) panel rather than erroring.
  axis_limit <- function(a) if(is.na(a) || !is.finite(a)) NA_real_ else a
  pos_scales_x <- lapply(df_stats$max_axis, function(a) scale_x_continuous(limits = c(0, axis_limit(a))))
  pos_scales_y <- lapply(df_stats$max_axis, function(a) scale_y_continuous(limits = c(0, axis_limit(a))))

  colours_nm <- colour_nm_func(sensor_Y)
  var_labs <- pretty_label_func("RHOW")

  fig <- ggplot(match_filter, aes(x = Hyp, y = Sat)) +
    geom_point(aes(colour = waveband_label), size = 2, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, colour = "black") +
    geom_text(data = df_stats, aes(label = label, y = label_y), x = 0,
              hjust = 0, vjust = 1, size = 2.5, inherit.aes = FALSE) +
    facet_wrap(site_name ~ waveband_label, nrow = n_sites, scales = "free", drop = FALSE) +
    scale_colour_manual(values = colours_nm, drop = TRUE) +
    ggh4x::facetted_pos_scales(x = pos_scales_x, y = pos_scales_y) +
    labs(x = paste0("HYPERNETS; ", var_labs$units_lab),
         y = paste0(sensor_Y, "; ", var_labs$units_lab),
         colour = "Wavelength (nm)") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          aspect.ratio = 1,
          legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.title.x = element_markdown(size = 12),
          axis.title.y = element_markdown(size = 12),
          axis.text = element_text(size = 9))

  n_col <- length(waveband_levels)
  ggsave(paste0("figures/global_scatter_waveband_RHOW_",sensor_Y,".png"), fig,
         width = 3 * n_col, height = 4 * n_sites + 1, limitsize = FALSE)
}

