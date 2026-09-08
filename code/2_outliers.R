# code/2_outliers.R
# Convenience functions to rapidly visualise outliers


# Setup -------------------------------------------------------------------

source("code/0_functions.R")

# Outlier gating
# Every sensor uses ONE live gate: the satellite pixel-variance CV check (sat_var_check()).
# MODIS/VIIRS/OLCI use cv_limit_choice below; OCI uses a separate hardcoded 50% threshold (see
# the OCI branch below) as a temporary workaround while all THFR PACE files have anomalously high
# CVs -- see Bug 8 in manuscript/upstream-data-bugs.md.
cv_limit_choice <- 30


# Per-sensor CV screening ---------------------------------------------------

sensor_Z_list <- c("MODIS", "VIIRS", "OLCI", "OCI")
outlier_list <- list()

for(sensor_Z in sensor_Z_list){

  print(paste0("Beginning ", sensor_Z, " outlier screening"))

  # Load processed in situ matchups
  matchup_sensor <- read_csv(paste0("output/matchup_stats_RHOW_", sensor_Z, ".csv"), show_col_types = FALSE) |>
    filter(sensor_X %in% c("Hyp")) |>
    mutate(comp_sensors = paste0(sensor_X," vs ",sensor_Y))

  # Sensor files -- pattern derived from sensor_grid() so it always matches the same sensor_Y
  # names used everywhere else in the pipeline, rather than a separately hardcoded regex per sensor
  sat_name_pattern <- paste(unique(sensor_grid(sensor_Z)$sensor_Y), collapse = "|")
  file_list_sensor <- stringr::str_subset(string = dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/",
                                           pattern = sat_name_pattern, full.names = TRUE, recursive = TRUE), pattern = "csv")

  # Check satellite variance in files -- the only read of each raw matchup CSV needed for
  # outlier screening (NB: High THFR AQUA/OCI CVs are genuine failed AC retrievals / a known
  # pixel-extraction offset -- see upstream-data-bugs.md)
  cv_limit_sensor <- if(sensor_Z == "OCI") 50 else cv_limit_choice # TODO reimplement OCI at cv_limit_choice once its high CV values are understood
  sat_var_sensor <- furrr::future_map_dfr(file_list_sensor, sat_var_check, cv_limit = cv_limit_sensor, .options = furrr_options(seed = TRUE)) |>
    mutate(site_name = path_site_name(file_name), file_name = basename(file_name))
  sat_var_filt_sensor <- filter(sat_var_sensor, cv > cv_limit)
  val_filter_sensor <- if(sensor_Z == "OCI") paste0("CV >= ", cv_limit_sensor) else paste0("CV >= ", cv_limit_sensor, "%")
  filter_var_sensor <- inner_join(matchup_sensor, dplyr::select(sat_var_filt_sensor, file_name, site_name),
                                  by = c("file_name", "site_name")) |>
    mutate(val_filter = val_filter_sensor)

  # TODO: meta/<site>_pixel_removals.csv (code/0_functions.R) now carries a structural key --
  # matchup_id, db_source, dateTime_Hyp, dateTime_sat -- that survives even once ingestion no
  # longer produces a per-matchup CSV (added 2026-09-07). Joining it against matchup_sensor here
  # (still file_name-keyed today) needs either matching on the derived filename
  # (db_matchup_filename(), code/0_functions.R) or -- once this script's refactor reaches it --
  # switching matchup_sensor's own key to the same structural (matchup_id, db_source) pair.

  outlier_list[[sensor_Z]] <- filter_var_sensor
}


# Combine satellite outliers ----------------------------------------------

print("Combining results and exiting")

# Stack all filtered data.frames with file names that appear to be outliers. site_name is kept
# alongside file_name (not just file_name) because raw matchup filenames are not unique across
# sites (see path_site_name() in code/0_functions.R) -- global_stats() filters this file by both.
satellite_outliers <- bind_rows(outlier_list) |>
  dplyr::select(file_name, site_name, sensor_X, sensor_Y, comp_sensors,
                dateTime_X, dateTime_Y, Slope_II, Error_50, Bias_50, val_filter) |>
  distinct()
write_csv(satellite_outliers, "meta/satellite_outliers.csv")

