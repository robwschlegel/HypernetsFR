# code/2_outliers.R
# Convenience functions to rapidly visualise outliers
#
# Pipeline run order: 0_functions.r -> 1_matchups_single.R -> 2_outliers.R ->
#                      3_sensitivity.R -> 4_matchups_global.R -> 5_figures.R
# (renamed 2026-07-10 from code/outliers.R -- see manuscript/track-changes.md)


# Setup -------------------------------------------------------------------

source("code/0_functions.r")

# Outlier gating (standardised 2026-07-10; see manuscript/track-changes.md)
# --------------------------------------------------------------------------
# For now every sensor uses ONE gate only: the satellite pixel-variance CV
# check (sat_var_check(), cv_limit = 30). The previously-used Error_50 >= 100%
# gate, and a prospective RMSE-based gate, are left in below as commented-out
# placeholders in case either is reinstated later -- do not delete them.
cv_limit_choice <- 30 # TODO: revisit; this is a placeholder value, not yet validated
error_50_limit <- 100 # TODO (placeholder, currently unused): may reinstate an Error_50-based gate later
rmse_limit <- NA # TODO (placeholder, currently unused): determine a principled RMSE-based threshold (units: RHOW) before enabling


# Satellite Outliers -------------------------------------------------------

## MODIS -------------------------------------------------------------------

# Load processed in situ matchups
matchup_MODIS <- read_csv("output/matchup_stats_RHOW_MODIS.csv") |>
  filter(sensor_X %in% c("Hyp")) |>
  mutate(comp_sensors = paste0(sensor_X," vs ",sensor_Y))

# MODIS files
file_list_MODIS <- stringr::str_subset(string = dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/",
                                         pattern = "AQUA", full.names = TRUE, recursive = TRUE), pattern = "csv")

# Load base W_nm matchup values
base_MODIS <- plyr::ldply(file_list_MODIS, load_matchup_long, .parallel = TRUE)

# Join for full range of stats
join_MODIS <- right_join(base_MODIS, matchup_MODIS, by = join_by(file_name))

# Check satellite variance in files
sat_var_MODIS <- plyr::ldply(file_list_MODIS, sat_var_check, cv_limit = cv_limit_choice, .parallel = TRUE)
sat_var_filt_MODIS <- filter(sat_var_MODIS, cv > cv_limit)
filter_var_MODIS <- filter(matchup_MODIS, file_name %in% basename(sat_var_filt_MODIS$file_name)) |>
  mutate(val_filter = paste0("CV >= ", cv_limit_choice, "%"))

# Plot matchup by Error + Bias
plot_matchup_Error_Bias(join_MODIS, "Hyp", "AQUA")

# Standardised outlier gate: CV only (see cv_limit_choice above)
filter_MODIS <- filter_var_MODIS

# Placeholder: previously-used Error_50 gate, kept for reference / possible future reactivation
# filter_MODIS <- matchup_MODIS |>
#   filter(!file_name %in% filter_var_MODIS$file_name) |>
#   filter(Error_50 >= error_50_limit) |> mutate(val_filter = paste0("Error >= ", error_50_limit, "%")) |>
#   bind_rows(filter_var_MODIS)

# Placeholder: possible future RMSE-based gate (threshold TBD)
# filter_RMSE_MODIS <- matchup_MODIS |>
#   filter(RMSE >= rmse_limit) |> mutate(val_filter = paste0("RMSE >= ", rmse_limit))
# filter_MODIS <- bind_rows(filter_MODIS, filter_RMSE_MODIS)

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

# VIIRS files
file_list_VIIRS <- stringr::str_subset(string = dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/",
                                         pattern = "VIIRS", full.names = TRUE, recursive = TRUE), pattern = "csv")

# Load base W_nm matchup values
base_VIIRS <- plyr::ldply(file_list_VIIRS, load_matchup_long, .parallel = TRUE)

# Join for full range of stats
join_VIIRS <- right_join(base_VIIRS, matchup_VIIRS, by = join_by(file_name))

# Check satellite variance in files
sat_var_VIIRS <- plyr::ldply(file_list_VIIRS, sat_var_check, cv_limit = cv_limit_choice, .parallel = TRUE)
sat_var_filt_VIIRS <- filter(sat_var_VIIRS, cv > cv_limit)
filter_var_VIIRS <- filter(matchup_VIIRS, file_name %in% basename(sat_var_filt_VIIRS$file_name)) |>
  mutate(val_filter = paste0("CV >= ", cv_limit_choice, "%"))

# Plot matchup by Error + Bias
# NB: SNPP is used as the reference platform for these diagnostic plots; the filtering below
# still screens outliers across all three VIIRS platforms (SNPP, JPSS1, JPSS2)
plot_matchup_Error_Bias(join_VIIRS, "Hyp", "SNPP")

# Standardised outlier gate: CV only (see cv_limit_choice above)
filter_VIIRS <- filter_var_VIIRS

# Placeholder: previously-used Error_50 gate, kept for reference / possible future reactivation
# filter_VIIRS <- matchup_VIIRS |>
#   filter(!file_name %in% filter_var_VIIRS$file_name) |>
#   filter(Error_50 >= error_50_limit) |> mutate(val_filter = paste0("Error >= ", error_50_limit, "%")) |>
#   bind_rows(filter_var_VIIRS)

# Placeholder: possible future RMSE-based gate (threshold TBD)
# filter_RMSE_VIIRS <- matchup_VIIRS |>
#   filter(RMSE >= rmse_limit) |> mutate(val_filter = paste0("RMSE >= ", rmse_limit))
# filter_VIIRS <- bind_rows(filter_VIIRS, filter_RMSE_VIIRS)

filter_join_VIIRS <- right_join(join_VIIRS, filter_VIIRS)
clean_join_VIIRS <- anti_join(join_VIIRS, filter_VIIRS)

# Plot matchups by date
plot_matchup_date(filter_join_VIIRS, "Hyp", "SNPP")

# Plot all wavelength matchups
plot_matchup_nm(join_VIIRS, "Hyp", "SNPP")
plot_matchup_nm(filter_join_VIIRS, "Hyp", "SNPP")
plot_matchup_nm(clean_join_VIIRS, "Hyp", "SNPP")


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
sat_var_OLCI <- plyr::ldply(file_list_OLCI, sat_var_check, cv_limit = cv_limit_choice, .parallel = TRUE) # Many high CV obs
sat_var_filt_OLCI <- filter(sat_var_OLCI, cv > cv_limit)
filter_var_OLCI <- filter(matchup_OLCI, file_name %in% basename(sat_var_filt_OLCI$file_name)) |>
  mutate(val_filter = paste0("CV >= ", cv_limit_choice, "%"))

# Plot matchup by Error + Bias
# NB: There are more S3A matchups, so using that sensor for analysis
plot_matchup_Error_Bias(join_OLCI, "Hyp", "S3B") # Error > 50

# Standardised outlier gate: CV only (see cv_limit_choice above)
filter_OLCI <- filter_var_OLCI

# Placeholder: previously-used Error_50 gate, kept for reference / possible future reactivation
# filter_OLCI <- matchup_OLCI |>
#   filter(!file_name %in% filter_var_OLCI$file_name) |>
#   filter(Error_50 >= error_50_limit) |> mutate(val_filter = paste0("Error >= ", error_50_limit, "%")) |>
#   bind_rows(filter_var_OLCI)

# Placeholder: possible future RMSE-based gate (threshold TBD)
# filter_RMSE_OLCI <- matchup_OLCI |>
#   filter(RMSE >= rmse_limit) |> mutate(val_filter = paste0("RMSE >= ", rmse_limit))
# filter_OLCI <- bind_rows(filter_OLCI, filter_RMSE_OLCI)

filter_join_OLCI <- right_join(join_OLCI, filter_OLCI)
clean_join_OLCI <- anti_join(join_OLCI, filter_OLCI)

# Plot matchups by date
plot_matchup_date(filter_join_OLCI, "Hyp", "S3B")

# Plot all wavelength matchups
plot_matchup_nm(join_OLCI, "Hyp", "S3B")
plot_matchup_nm(filter_join_OLCI, "Hyp", "S3B")
plot_matchup_nm(clean_join_OLCI, "Hyp", "S3B")


## OCI ---------------------------------------------------------------------

# Load processed in situ matchups
matchup_OCI <- read_csv("output/matchup_stats_RHOW_OCI.csv") |>
  filter(sensor_X %in% c("Hyp")) |>
  mutate(comp_sensors = paste0(sensor_X," vs ",sensor_Y))

# OCI files
# NB: Unlike the Tara matchup results, the FR folder layout only ever pairs a satellite against
# HYPERNETS (see file_path_build()), so there are no PACE-vs-PACE self-comparison folders to exclude
file_list_OCI <- stringr::str_subset(string = dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/",
                                         pattern = "PACE", full.names = TRUE, recursive = TRUE), pattern = "csv")

# Load base W_nm matchup values
base_OCI <- plyr::ldply(file_list_OCI, load_matchup_long, .parallel = TRUE)

# Join for full range of stats
join_OCI <- right_join(base_OCI, matchup_OCI, by = join_by(file_name))

# Check satellite variance in files
sat_var_OCI <- plyr::ldply(file_list_OCI, sat_var_check, cv_limit = cv_limit_choice, .parallel = TRUE)
sat_var_filt_OCI <- filter(sat_var_OCI, cv > cv_limit)
filter_var_OCI <- filter(matchup_OCI, file_name %in% basename(sat_var_filt_OCI$file_name)) |>
  mutate(val_filter = paste0("CV >= ", cv_limit_choice, "%"))

# Plot matchup by Error + Bias
# NB: PACE_V30 is visually the least similar, so using this platform as the reference; the
# filtering below still screens outliers across all three PACE versions (V2, V30, V31)
plot_matchup_Error_Bias(join_OCI, "Hyp", "PACE_V30")

# Standardised outlier gate: CV only (see cv_limit_choice above)
filter_OCI <- filter_var_OCI

# Placeholder: previously-used Error_50 gate, kept for reference / possible future reactivation
# filter_OCI <- matchup_OCI |>
#   filter(!file_name %in% filter_var_OCI$file_name) |>
#   filter(Error_50 >= error_50_limit) |> mutate(val_filter = paste0("Error >= ", error_50_limit, "%")) |>
#   bind_rows(filter_var_OCI)

# Placeholder: possible future RMSE-based gate (threshold TBD)
# filter_RMSE_OCI <- matchup_OCI |>
#   filter(RMSE >= rmse_limit) |> mutate(val_filter = paste0("RMSE >= ", rmse_limit))
# filter_OCI <- bind_rows(filter_OCI, filter_RMSE_OCI)

filter_join_OCI <- right_join(join_OCI, filter_OCI)
clean_join_OCI <- anti_join(join_OCI, filter_OCI)

# Plot matchups by date
plot_matchup_date(filter_join_OCI, "Hyp", "PACE_V30")

# Plot all wavelength matchups
plot_matchup_nm(join_OCI, "Hyp", "PACE_V30")
plot_matchup_nm(filter_join_OCI, "Hyp", "PACE_V30")
plot_matchup_nm(clean_join_OCI, "Hyp", "PACE_V30")


## Combine satellite outliers ----------------------------------------------

# Stack all filtered data.frames with file names that appear to be outliers
satellite_outliers <- rbind(filter_OLCI, filter_VIIRS, filter_MODIS, filter_OCI) |>
  dplyr::select(file_name, sensor_X, sensor_Y, var_name, comp_sensors,
                dateTime_X, dateTime_Y, Slope_II, Error_50, Bias_50, val_filter) |>
  distinct()
write_csv(satellite_outliers, "meta/satellite_outliers.csv")
