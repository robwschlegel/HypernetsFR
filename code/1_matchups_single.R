# code/1_matchups_single.R
# Compute per-file (individual) match-up statistics for every sensor family.


# Setup -------------------------------------------------------------------

source("code/0_functions.R")


# TEMPORARY: NE-quadrant-filtered THFR matchups ---------------------------
# Addresses the off-center satellite pixel-box extraction found in meta/pixel_explore.R (see
# meta/pixel_explore_output/summary.md). Regenerates raw per-matchup RHOW CSVs under a new
# THFR_NE site folder (db_export_matchups_ne() / write_matchup_csv_ne() in code/0_functions.R),
# kept fully independent of THFR/MAFR (not a replacement for THFR anywhere). THFR_NE is now a
# first-class site alongside MAFR/THFR (added to available_sites()'s candidate list) so it
# must be generated here, before process_sensor() below, which picks it up automatically via
# sensor_grid()/available_sites() the same way it already does for MAFR/THFR.

for(sZ in c("MODIS", "VIIRS", "OLCI", "OCI")) db_export_matchups_ne(sZ)


# Individual matchup stats ------------------------------------------------

process_sensor("MODIS")
process_sensor("VIIRS")
process_sensor("OLCI")
process_sensor("OCI")


# Summary single matchups stats ------------------------------------------

# Re-load all single matchups
matchup_single_all <- map_dfr(dir("output", pattern = "matchup_stats_", full.names = TRUE), read_csv)

# Date and time range of samples per sensor
matchup_date_time_range <- matchup_single_all |>
  dplyr::select(sensor_X, dateTime_X) |>
  distinct() |>
  mutate(date = as.Date(dateTime_X),
         time = format(dateTime_X, format = "%H:%M:%S")) |>
  summarise(date_min = min(date), date_max = max(date),
            time_min = min(time), time_max = max(time), .by = "sensor_X")

# Unique number of satellite passes available for each platform+sensor/version
matchup_sat_uniq <- matchup_single_all |>
  dplyr::select(sensor_X, dateTime_X) |>
  filter(!(sensor_X %in% c("Hyp", "TRIOS", "HYPERPRO"))) |>
  distinct() |>
  summarise(sat_count = n(), .by = "sensor_X")

