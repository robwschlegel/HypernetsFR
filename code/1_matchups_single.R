# code/1_matchups_single.R
# Compute per-file (individual) match-up statistics for every sensor family.


# Setup -------------------------------------------------------------------

source("code/0_functions.R")


# TEMPORARY: derived THFR/MAFR matchup sites -------------------------------
# Addresses the off-center satellite pixel-box extraction found in meta/pixel_explore.R (see
# meta/pixel_explore_output/summary.md). Regenerates raw per-matchup RHOW CSVs under new site
# folders (db_export_matchups_site() / write_matchup_csv_db() in code/0_functions.R), kept fully
# independent of THFR/MAFR (not a replacement for THFR anywhere): THFR_NE restricts to the NE
# quadrant of the inner 3x3 pixel grid, THFR_poly restricts to the hand-drawn clean-water polygon,
# and THFR_pixel restricts only to the inner 3x3 pixel grid with no further spatial subsetting
# (isolates the effect of the shared per-pixel QC gates in db_export_matchups_site() -- RHOW
# ceiling, negative-value, distance, and minimum-valid-pixel-count -- from the spatial filters
# used by THFR_NE/THFR_poly). MAFR_pixel (added 2026-09-07) is the analogous no-spatial-filter
# reconstruction of MAFR, sourced from mafr_2024.db/mafr_2025.db via db_export_matchups_multi().
# THFR_raw/MAFR_raw (added 2026-09-07) go one step further again: the same 3x3-box reconstruction
# with apply_pixel_qc = FALSE, i.e. no pixel-level QC gates at all -- meant to reproduce
# Hypernets_matchups' own unfiltered aggregation, so validate_derived_site() (see
# code/3_sensitivity.R) has a fair like-for-like check before MAFR/THFR themselves are ever
# retargeted onto .db-native data. All six are first-class sites alongside MAFR/THFR (added to
# available_sites()'s candidate list) so they must be generated here, before process_sensor()
# below, which picks them up automatically via sensor_grid()/available_sites() the same way it
# already does for MAFR/THFR. Each call loops internally over all 4 sensor families (rather than
# being called once per sensor_Z here) so it can accumulate one pixel-removal audit CSV per site
# (meta/<site>_pixel_removals.csv) across every sensor family in a single write.

db_export_matchups_ne()
db_export_matchups_poly()
db_export_matchups_pixel()
db_export_matchups_mafr_pixel()
db_export_matchups_thfr_raw()
db_export_matchups_mafr_raw()


# Individual matchup stats ------------------------------------------------

process_sensor("MODIS")
process_sensor("VIIRS")
process_sensor("OLCI")
process_sensor("OCI")


# Summary single matchups stats ------------------------------------------

# Re-load all single matchups
matchup_single_all <- map_dfr(dir("output", pattern = "matchup_stats_", full.names = TRUE), read_csv)

# Date and time range of samples per sensor
# NB: sensor_Y/dateTime_Y (not sensor_X/dateTime_X) since process_matchup_file() now computes
# stats exactly once per file with sensor_X fixed to "Hyp" and sensor_Y the satellite (2026-09-03 --
# previously also computed the reverse direction, which is what this used to key off of).
matchup_date_time_range <- matchup_single_all |>
  dplyr::select(sensor_Y, dateTime_Y) |>
  distinct() |>
  mutate(date = as.Date(dateTime_Y),
         time = format(dateTime_Y, format = "%H:%M:%S")) |>
  summarise(date_min = min(date), date_max = max(date),
            time_min = min(time), time_max = max(time), .by = "sensor_Y")

# Unique number of satellite passes available for each platform+sensor/version
matchup_sat_uniq <- matchup_single_all |>
  dplyr::select(sensor_Y, dateTime_Y) |>
  distinct() |>
  summarise(sat_count = n(), .by = "sensor_Y")

