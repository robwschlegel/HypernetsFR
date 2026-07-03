# code/outliers.R
# Convenience functions to rapidly visualise outliers


# Setup -------------------------------------------------------------------

source("code/functions.R")


# Satellite Outliers -------------------------------------------------------

## MODIS -------------------------------------------------------------------

# # Load processed in situ matchups
# matchup_MODIS <- read_csv("output/matchup_stats_RHOW_MODIS.csv") |> 
#   filter(sensor_X %in% c("Hyp", "HYPERPRO", "TRIOS")) |> 
#   mutate(comp_sensors = paste0(sensor_X," vs ",sensor_Y))

# # File list
# file_list_MODIS <- list.files(dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260203", 
#                                   pattern = "AQUA", full.names = TRUE), pattern = "*.csv", full.names = TRUE)

# # Load base W_nm matchup values
# base_MODIS <- plyr::ldply(file_list_MODIS, load_matchup_long, .parallel = TRUE)

# # Join for full range of stats
# join_MODIS <- right_join(base_MODIS, matchup_MODIS, by = join_by(file_name))

# # Plot matchup by Error + Bias
# plot_matchup_Error_Bias(join_MODIS, "Rhow", "Hyp", "AQUA") # Error > 50
# plot_matchup_Error_Bias(join_MODIS, "Rhow", "TRIOS", "AQUA") # OK
# plot_matchup_Error_Bias(join_MODIS, "Rhow", "HYPERPRO", "AQUA") # OK

# # Check satellite variance in files
# sat_var_MODIS <- plyr::ldply(file_list_MODIS, sat_var_check, .parallel = TRUE) # OK

# # Filter all by Error or Bias to get an initial idea of the issues
# # No need to filter by satellite variance, there is one bad AQUA matchup across all sensors and dateTimes
# # 2024-08-15 10:45:01
# filter_MODIS <- filter(matchup_MODIS, Error >= 50) |> mutate(val_filter = "Error >= 50%") |> 
#   filter((dateTime_Y == "2024-08-15 10:45:01" & sensor_Y == "AQUA")) # Manually checked, not an outlier, just a poor matchup
# filter_join_MODIS <- right_join(join_MODIS, filter_MODIS)
# clean_join_MODIS <- anti_join(join_MODIS, filter_MODIS)

# # Plot matchups by date
# plot_matchup_date(filter_join_MODIS, "Rhow", "Hyp", "AQUA") # OK
# plot_matchup_date(filter_join_MODIS, "Rhow", "TRIOS", "AQUA") # OK
# plot_matchup_date(filter_join_MODIS, "Rhow", "HYPERPRO", "AQUA") # OK

# # Plot all wavelength matchups
# plot_matchup_nm(join_MODIS, "Rhow", "Hyp", "AQUA")
# plot_matchup_nm(filter_join_MODIS, "Rhow", "Hyp", "AQUA")
# plot_matchup_nm(clean_join_MODIS, "Rhow", "Hyp", "AQUA")
# plot_matchup_nm(join_MODIS, "Rhow", "TRIOS", "AQUA")
# plot_matchup_nm(filter_join_MODIS, "Rhow", "TRIOS", "AQUA")
# plot_matchup_nm(clean_join_MODIS, "Rhow", "TRIOS", "AQUA")
# plot_matchup_nm(join_MODIS, "Rhow", "HYPERPRO", "AQUA")
# plot_matchup_nm(filter_join_MODIS, "Rhow", "HYPERPRO", "AQUA") # OK
# plot_matchup_nm(clean_join_MODIS, "Rhow", "HYPERPRO", "AQUA")


## VIIRS -------------------------------------------------------------------

# # Load processed in situ matchups
# matchup_VIIRS <- read_csv("output/matchup_stats_RHOW_VIIRS.csv") |> 
#   filter(sensor_X %in% c("Hyp", "HYPERPRO", "TRIOS")) |> 
#   mutate(comp_sensors = paste0(sensor_X," vs ",sensor_Y))

# # File list
# file_list_VIIRS <- list.files(dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260203", 
#                                   pattern = "VIIRS", full.names = TRUE), pattern = "*.csv", full.names = TRUE)

# # Load base W_nm matchup values
# base_VIIRS <- plyr::ldply(file_list_VIIRS, load_matchup_long, .parallel = TRUE)

# # Join for full range of stats
# join_VIIRS <- right_join(base_VIIRS, matchup_VIIRS, by = join_by(file_name))
                           
# # Plot matchup by Error + Bias
# # NB: VIIRS_N is visually the least similar, so using this for base reference
# plot_matchup_Error_Bias(join_VIIRS, "Rhow", "Hyp", "VIIRS_N") # Error > 50
# plot_matchup_Error_Bias(join_VIIRS, "Rhow", "TRIOS", "VIIRS_N") # OK
# plot_matchup_Error_Bias(join_VIIRS, "Rhow", "HYPERPRO", "VIIRS_N") # OK

# # Check satellite variance in files
# sat_var_VIIRS <- plyr::ldply(file_list_VIIRS, sat_var_check, .parallel = TRUE)
# filter_var_VIIRS <- filter(matchup_VIIRS, file_name %in% sat_var_VIIRS$file_name) |> mutate(val_filter = "CV >= 20%")

# # Filter all by Error or Bias to get an initial idea of the issues
# filter_VIIRS <- matchup_VIIRS |> 
#   filter(!file_name %in% filter_var_VIIRS$file_name) |> 
#   filter(Error >= 50) |> mutate(val_filter = "Error >= 50%") |> 
#     filter(!file_name %in% c("HYPERNETS_vs_VIIRS_N_vs_20240813T102100_RHOW.csv", 
#                              "TRIOS_vs_VIIRS_N_vs_20240813T100254_RHOW.csv",
#                              "TRIOS_vs_VIIRS_N_vs_20240813T101805_RHOW.csv",
#                              "TRIOS_vs_VIIRS_N_vs_20240813T102853_RHOW.csv")) |> # Manually checked, not an outlier
#   bind_rows(filter_var_VIIRS)
# filter_join_VIIRS <- right_join(join_VIIRS, filter_VIIRS)
# clean_join_VIIRS <- anti_join(join_VIIRS, filter_VIIRS)

# # Plot matchups by date
# plot_matchup_date(filter_join_VIIRS, "Rhow", "Hyp", "VIIRS_N")
# plot_matchup_date(filter_join_VIIRS, "Rhow", "TRIOS", "VIIRS_N")
# plot_matchup_date(filter_join_VIIRS, "Rhow", "HYPERPRO", "VIIRS_N") # OK

# # Plot all wavelength matchups
# plot_matchup_nm(join_VIIRS, "Rhow", "Hyp", "VIIRS_N")
# plot_matchup_nm(filter_join_VIIRS, "Rhow", "Hyp", "VIIRS_N")
# plot_matchup_nm(clean_join_VIIRS, "Rhow", "Hyp", "VIIRS_N")
# plot_matchup_nm(join_VIIRS, "Rhow", "TRIOS", "VIIRS_N")
# plot_matchup_nm(filter_join_VIIRS, "Rhow", "TRIOS", "VIIRS_N")
# plot_matchup_nm(clean_join_VIIRS, "Rhow", "TRIOS", "VIIRS_N")
# plot_matchup_nm(join_VIIRS, "Rhow", "HYPERPRO", "VIIRS_N")
# plot_matchup_nm(filter_join_VIIRS, "Rhow", "HYPERPRO", "VIIRS_N") # OK
# plot_matchup_nm(clean_join_VIIRS, "Rhow", "HYPERPRO", "VIIRS_N")


## OLCI --------------------------------------------------------------------

# Load processed in situ matchups
matchup_OLCI <- read_csv("output/matchup_stats_RHOW_OLCI.csv") |> 
  filter(sensor_X %in% c("Hyp")) |> 
  mutate(comp_sensors = paste0(sensor_X," vs ",sensor_Y),
         file_name = gsub("_V4", "", file_name)) # Rename for easier joining with base data

# OLCI files
file_list_OLCI <- stringr::str_subset(string = dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/", 
                                         pattern = "S3", full.names = TRUE, recursive = TRUE), pattern = "csv")

# Load base W_nm matchup values
base_OLCI <- plyr::ldply(file_list_OLCI, load_matchup_long, .parallel = TRUE)

# Join for full range of stats
join_OLCI <- right_join(base_OLCI, matchup_OLCI, by = join_by(file_name))

# Check satellite variance in files
sat_var_OLCI <- plyr::ldply(file_list_OLCI, sat_var_check, .parallel = TRUE) # Many high CV obs
sat_var_filt_OLCI <- filter(sat_var_OLCI, cv > cv_limit)
filter_var_OLCI <- filter(matchup_OLCI, file_name %in% basename(sat_var_filt_OLCI$file_name)) |> 
  mutate(val_filter = "CV >= 30%")

# Plot matchup by Error + Bias
# NB: There are more S3A matchups, so using that sensor for analysis
plot_matchup_Error_Bias(join_OLCI, "Rhow", "Hyp", "S3B") # Error > 50

# Filter all by Error or Bias to get an initial idea of the issues
filter_OLCI <- matchup_OLCI |> 
  filter(!file_name %in% filter_var_OLCI$file_name) |> 
  filter(Error_50 >= 100) |> mutate(val_filter = "Error >= 100%") |> 
  # filter(!file_name %in% c("HYPERNETS_vs_S3A_20240814T123100_RHOW.csv"),
  #        !as.character(dateTime_Y) %in% c("2024-08-14 10:05:12", "2024-08-15 09:06:24")) |> # Manually checked, not an outlier
  bind_rows(filter_var_OLCI)
filter_join_OLCI <- right_join(join_OLCI, filter_OLCI)
clean_join_OLCI <- anti_join(join_OLCI, filter_OLCI)

# Plot matchups by date
plot_matchup_date(filter_join_OLCI, "Rhow", "Hyp", "S3B")

# Plot all wavelength matchups
plot_matchup_nm(join_OLCI, "Rhow", "Hyp", "S3B")
plot_matchup_nm(filter_join_OLCI, "Rhow", "Hyp", "S3B")
plot_matchup_nm(clean_join_OLCI, "Rhow", "Hyp", "S3B")


## OCI ---------------------------------------------------------------------

# # Load processed in situ matchups
# matchup_OCI <- read_csv("output/matchup_stats_RHOW_OCI.csv") |> 
#   filter(sensor_X %in% c("Hyp", "HYPERPRO", "TRIOS")) |> 
#   mutate(comp_sensors = paste0(sensor_X," vs ",sensor_Y))

# # File list
# file_list_OCI <- list.files(dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260203", 
#                                   pattern = "PACE", full.names = TRUE), pattern = "*.csv", full.names = TRUE)
# file_list_OCI <- file_list_OCI[!grepl("RHOW_PACE_V2_vs_PACE_V30_vs_PACE_V31", file_list_OCI)] # Remove self-comparisons

# # Load base W_nm matchup values
# base_OCI <- plyr::ldply(file_list_OCI, load_matchup_long, .parallel = TRUE)

# # Join for full range of stats
# join_OCI <- right_join(base_OCI, matchup_OCI, by = join_by(file_name))

# # Plot matchup by Error + Bias
# # NB: PACE_v30 is visually the least similar, so using this for base reference
# # NB: There are many PACE files with negative values
# plot_matchup_Error_Bias(join_OCI, "Rhow", "Hyp", "PACE_V31") # Bias < -50
# plot_matchup_Error_Bias(join_OCI, "Rhow", "TRIOS", "PACE_V31") # Bias < -50
# plot_matchup_Error_Bias(join_OCI, "Rhow", "HYPERPRO", "PACE_V31") # Bias < -50

# # Check satellite variance in files
# sat_var_OCI <- plyr::ldply(file_list_OCI, sat_var_check, .parallel = TRUE)
# filter_var_OCI <- filter(matchup_OCI, file_name %in% sat_var_OCI$file_name) |> mutate(val_filter = "CV >= 20%")

# # Filter all by Error or Bias to get an initial idea of the issues
# filter_OCI <- matchup_OCI |> 
#   filter(!file_name %in% filter_var_OCI$file_name) |> 
#   filter(Error >= 50) |> mutate(val_filter = "Error >= 50%") |> 
#   filter(!file_name %in% c("HYPERNETS_vs_PACE_V31_vs_20240814T123100_RHOW.csv"),
#          !as.character(dateTime_Y) %in% c("2024-08-14 10:05:12", "2024-08-15 09:06:24")) |> # Manually checked, not an outlier
#   bind_rows(filter_var_OCI)
# filter_join_OCI <- right_join(base_OCI, filter_OCI, by = join_by(file_name))
# clean_join_OCI <- anti_join(base_OCI, filter_OCI, by = join_by(file_name))

# # Plot matchups by date
# plot_matchup_date(filter_join_OCI, "Rhow", "Hyp", "PACE_V2")
# plot_matchup_date(filter_join_OCI, "Rhow", "TRIOS", "PACE_V2")
# plot_matchup_date(filter_join_OCI, "Rhow", "HYPERPRO", "PACE_V2")

# # Plot all wavelength matchups
# plot_matchup_nm(join_OCI, "Rhow", "Hyp", "PACE_V30")
# plot_matchup_nm(filter_join_OCI, "Rhow", "Hyp", "PACE_V30")
# plot_matchup_nm(clean_join_OCI, "Rhow", "Hyp", "PACE_V30")
# plot_matchup_nm(join_OCI, "Rhow", "TRIOS", "PACE_V30")
# plot_matchup_nm(filter_join_OCI, "Rhow", "TRIOS", "PACE_V30")
# plot_matchup_nm(clean_join_OCI, "Rhow", "TRIOS", "PACE_V30")
# plot_matchup_nm(join_OCI, "Rhow", "HYPERPRO", "PACE_V30")
# plot_matchup_nm(filter_join_OCI, "Rhow", "HYPERPRO", "PACE_V30")
# plot_matchup_nm(clean_join_OCI, "Rhow", "HYPERPRO", "PACE_V30")


## Combine satellite outliers ----------------------------------------------

# TODO: Reactivate once all satellite matchups are sorted
# Stack all filtered data.frames with file names that appear to be outliers
# satellite_outliers <- rbind(filter_OLCI, filter_VIIRS, filter_MODIS, filter_OCI) |> 
#   dplyr::select(file_name, sensor_X, sensor_Y, var_name, comp_sensors, 
#                 dateTime_X, dateTime_Y, Slope, Error, Bias, val_filter) |> 
#   distinct()
satellite_outliers <- filter_OLCI |> 
  dplyr::select(file_name, sensor_X, sensor_Y, var_name, comp_sensors, 
                dateTime_X, dateTime_Y, Slope_II, Error_50, Bias_50, val_filter) |> 
  distinct()
write_csv(satellite_outliers, "meta/satellite_outliers.csv")

