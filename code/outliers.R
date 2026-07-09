# code/outliers.R
# Convenience functions to rapidly visualise outliers


# Setup -------------------------------------------------------------------

source("code/functions.R")


# Satellite Outliers -------------------------------------------------------

## MODIS -------------------------------------------------------------------

# # Load processed in situ matchups
matchup_MODIS <- read_csv("output/matchup_stats_RHOW_MODIS.csv") |>
  filter(sensor_X %in% c("Hyp")) |>
  mutate(comp_sensors = paste0(sensor_X," vs ",sensor_Y))

# # MODIS files
file_list_MODIS <- stringr::str_subset(string = dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/",
                                         pattern = "AQUA", full.names = TRUE, recursive = TRUE), pattern = "csv")

# # Load base W_nm matchup values
base_MODIS <- plyr::ldply(file_list_MODIS, load_matchup_long, .parallel = TRUE)

# # Join for full range of stats
join_MODIS <- right_join(base_MODIS, matchup_MODIS, by = join_by(file_name))

# Check satellite variance in files
sat_var_MODIS <- plyr::ldply(file_list_MODIS, sat_var_check, .parallel = TRUE)
sat_var_filt_MODIS <- filter(sat_var_MODIS, cv > cv_limit)
filter_var_MODIS <- filter(matchup_MODIS, file_name %in% basename(sat_var_filt_MODIS$file_name)) |>
  mutate(val_filter = "CV >= 30%")

# Plot matchup by Error + Bias
plot_matchup_Error_Bias(join_MODIS, "Hyp", "AQUA")

# Filter all by Error or Bias to get an initial idea of the issues
filter_MODIS <- matchup_MODIS |>
  filter(!file_name %in% filter_var_MODIS$file_name) |>
  filter(Error_50 >= 100) |> mutate(val_filter = "Error >= 100%") |>
  bind_rows(filter_var_MODIS)
filter_join_MODIS <- right_join(join_MODIS, filter_MODIS)
clean_join_MODIS <- anti_join(join_MODIS, filter_MODIS)

# Plot matchups by date
plot_matchup_date(filter_join_MODIS, "Hyp", "AQUA")

# Plot all wavelength matchups
plot_matchup_nm(join_MODIS, "Hyp", "AQUA")
plot_matchup_nm(filter_join_MODIS, "Hyp", "AQUA")
plot_matchup_nm(clean_join_MODIS, "Hyp", "AQUA")


## VIIRS -------------------------------------------------------------------

# Load processed in situ matchups
matchup_VIIRS <- read_csv("output/matchup_stats_RHOW_VIIRS.csv") |>
  filter(sensor_X %in% c("Hyp")) |>
  mutate(comp_sensors = paste0(sensor_X," vs ",sensor_Y))

# # VIIRS files
file_list_VIIRS <- stringr::str_subset(string = dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/",
                                         pattern = "VIIRS", full.names = TRUE, recursive = TRUE), pattern = "csv")

# Load base W_nm matchup values
base_VIIRS <- plyr::ldply(file_list_VIIRS, load_matchup_long, .parallel = TRUE)

# Join for full range of stats
join_VIIRS <- right_join(base_VIIRS, matchup_VIIRS, by = join_by(file_name))

# # Check satellite variance in files
# sat_var_VIIRS <- plyr::ldply(file_list_VIIRS, sat_var_check, .parallel = TRUE)
# sat_var_filt_VIIRS <- filter(sat_var_VIIRS, cv > cv_limit)
# filter_var_VIIRS <- filter(matchup_VIIRS, file_name %in% basename(sat_var_filt_VIIRS$file_name)) |>
#   mutate(val_filter = "CV >= 30%")

# # Plot matchup by Error + Bias
# # NB: SNPP is used as the reference platform for these diagnostic plots; the filtering below
# # still screens outliers across all three VIIRS platforms (SNPP, JPSS1, JPSS2)
# plot_matchup_Error_Bias(join_VIIRS, "Hyp", "SNPP")

# # Filter all by Error or Bias to get an initial idea of the issues
# filter_VIIRS <- matchup_VIIRS |>
#   filter(!file_name %in% filter_var_VIIRS$file_name) |>
#   filter(Error_50 >= 100) |> mutate(val_filter = "Error >= 100%") |>
#   bind_rows(filter_var_VIIRS)
# filter_join_VIIRS <- right_join(join_VIIRS, filter_VIIRS)
# clean_join_VIIRS <- anti_join(join_VIIRS, filter_VIIRS)

# # Plot matchups by date
# plot_matchup_date(filter_join_VIIRS, "Hyp", "SNPP")

# # Plot all wavelength matchups
# plot_matchup_nm(join_VIIRS, "Hyp", "SNPP")
# plot_matchup_nm(filter_join_VIIRS, "Hyp", "SNPP")
# plot_matchup_nm(clean_join_VIIRS, "Hyp", "SNPP")


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
join_OLCI <- right_join(base_OLCI, matchup_OLCI, by = join_by(file_name)) |> 
  filter(Hyp <= 1)

# Check satellite variance in files
sat_var_OLCI <- plyr::ldply(file_list_OLCI, sat_var_check, .parallel = TRUE) # Many high CV obs
sat_var_filt_OLCI <- filter(sat_var_OLCI, cv > cv_limit)
filter_var_OLCI <- filter(matchup_OLCI, file_name %in% basename(sat_var_filt_OLCI$file_name)) |> 
  mutate(val_filter = "CV >= 30%")

# Plot matchup by Error + Bias
# NB: There are more S3A matchups, so using that sensor for analysis
plot_matchup_Error_Bias(join_OLCI, "Hyp", "S3B") # Error > 50

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
plot_matchup_date(filter_join_OLCI, "Hyp", "S3B")

# Plot all wavelength matchups
plot_matchup_nm(join_OLCI, "Hyp", "S3B")
plot_matchup_nm(filter_join_OLCI, "Hyp", "S3B")
plot_matchup_nm(clean_join_OLCI, "Hyp", "S3B")

# Extract the extreme outliers for further investigation
# filter_extreme_OLCI <- matchup_OLCI |> 
#   filter(as.Date(dateTime_X) %in% c("2024-06-17", "2024-07-29", "2024-09-02",
#                                     "2024-11-10", "2024-11-11", "2024-11-14",
#                                     "2025-07-15", "2025-08-16", "2025-10-27", "2025-12-05"))
#   # filter(Bias_50 >= 100 | Bias_50 <= -200)
# filter_join_extreme_OLCI <- right_join(join_OLCI, filter_extreme_OLCI)
# plot_matchup_date(filter_join_extreme_OLCI, "Hyp", "S3B")

# Save as a list
# filter_join_extreme_OLCI_list <- filter_join_extreme_OLCI |> 
#   dplyr::select(file_name) |> 
#   distinct() |> 
#   pull(file_name)
# write_csv(data.frame(file_name = filter_join_extreme_OLCI_list), "meta/filter_join_extreme_OLCI_list.csv")


## OCI ---------------------------------------------------------------------

# # Load processed in situ matchups
# matchup_OCI <- read_csv("output/matchup_stats_RHOW_OCI.csv") |>
#   filter(sensor_X %in% c("Hyp")) |>
#   mutate(comp_sensors = paste0(sensor_X," vs ",sensor_Y))

# # OCI files
# # NB: Unlike the Tara matchup results, the FR folder layout only ever pairs a satellite against
# # HYPERNETS (see file_path_build()), so there are no PACE-vs-PACE self-comparison folders to exclude
# file_list_OCI <- stringr::str_subset(string = dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/",
#                                          pattern = "PACE", full.names = TRUE, recursive = TRUE), pattern = "csv")

# # Load base W_nm matchup values
# base_OCI <- plyr::ldply(file_list_OCI, load_matchup_long, .parallel = TRUE)

# # Join for full range of stats
# join_OCI <- right_join(base_OCI, matchup_OCI, by = join_by(file_name))

# # Check satellite variance in files
# sat_var_OCI <- plyr::ldply(file_list_OCI, sat_var_check, .parallel = TRUE)
# sat_var_filt_OCI <- filter(sat_var_OCI, cv > cv_limit)
# filter_var_OCI <- filter(matchup_OCI, file_name %in% basename(sat_var_filt_OCI$file_name)) |>
#   mutate(val_filter = "CV >= 30%")

# # Plot matchup by Error + Bias
# # NB: PACE_V30 is visually the least similar, so using this platform as the reference; the
# # filtering below still screens outliers across all three PACE versions (V2, V30, V31)
# plot_matchup_Error_Bias(join_OCI, "Hyp", "PACE_V30")

# # Filter all by Error or Bias to get an initial idea of the issues
# filter_OCI <- matchup_OCI |>
#   filter(!file_name %in% filter_var_OCI$file_name) |>
#   filter(Error_50 >= 100) |> mutate(val_filter = "Error >= 100%") |>
#   bind_rows(filter_var_OCI)
# filter_join_OCI <- right_join(join_OCI, filter_OCI)
# clean_join_OCI <- anti_join(join_OCI, filter_OCI)

# # Plot matchups by date
# plot_matchup_date(filter_join_OCI, "Hyp", "PACE_V30")

# # Plot all wavelength matchups
# plot_matchup_nm(join_OCI, "Hyp", "PACE_V30")
# plot_matchup_nm(filter_join_OCI, "Hyp", "PACE_V30")
# plot_matchup_nm(clean_join_OCI, "Hyp", "PACE_V30")


## Combine satellite outliers ----------------------------------------------

# TODO: Reactivate once the MODIS, VIIRS, and OCI matchups have been processed (see code/matchups.R)
# Stack all filtered data.frames with file names that appear to be outliers
# satellite_outliers <- rbind(filter_OLCI, filter_VIIRS, filter_MODIS, filter_OCI) |>
#   dplyr::select(file_name, sensor_X, sensor_Y, var_name, comp_sensors,
#                 dateTime_X, dateTime_Y, Slope_II, Error_50, Bias_50, val_filter) |>
#   distinct()
satellite_outliers <- filter_OLCI |>
  dplyr::select(file_name, sensor_X, sensor_Y, var_name, comp_sensors, 
                dateTime_X, dateTime_Y, Slope_II, Error_50, Bias_50, val_filter) |> 
  distinct()
write_csv(satellite_outliers, "meta/satellite_outliers.csv")

