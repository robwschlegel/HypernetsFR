# meta/pixel_explore.R
# Exploratory, line-by-line script for investigating why THFR's global
# scatterplots are so much worse than MAFR's, for every sensor family.
#
# Working hypothesis: oyster beds to the north-west, west, and south-west of 
# the THFR station contaminate some of the pixels in the extracted 5x5 satellite 
# pixel box (measure_data_rhow, pixel_pos = "col,row", col/row in -2..2), and 
# whichever pixels get averaged into the matchup value are what's dragging THFR's
# Hyp-vs-Sat comparison off the 1:1 line.
#
# This is sourced by any pipeline script and is meant to berun interactively, 
# section by section, and edit the "EDIT ME" variables in section 4 to inspect 
# different dates/sensors by hand. Uses db_matchup_pixels()/db_load_spectra_pixel()
# (code/0_functions.R) added alongside this script, which expose the raw,
# unaggregated per-pixel rows (db_matchup_long()/db_load_spectra() only ever
# return the pixel-box mean or center pixel).


# Setup ---------------------------------------------------------------------

source("code/0_functions.R")
library(sf)
library(maptiles)
library(tidyterra)

db_path <- "~/pCloudDrive/Documents/OMTAB/HYPERNETS/FR/thfr_2025.db"
out_dir <- "meta/pixel_explore_output"
dir.create(out_dir, showWarnings = FALSE)

# Sensors actually present for THFR in this db (db_matchup_pixels() warns and
# returns an empty tibble for any sensor_Y with no matchups, so it's safe to
# just try every db_satellite_names entry rather than pre-filtering)
sensor_Y_thfr <- db_satellite_names


# Bearing/quadrant helpers ---------------------------------------------------

# Classify each satellite pixel's position relative to the HYPERNETS station
# (lon_Hyp/lat_Hyp, the per-scan GPS position stored in measure_info used
# in preference to meta/station_in_situ.csv's fixed THFR coordinate, which is
# ~1.9 km off from the db's actual per-scan position, see section 3 below).
# quadrant is the primary (hard-cut) grouping. Bearing_deg (0-360, clockwise
# from north) is kept as a continuous robustness check since "north and west"
# doesn't cleanly bisect the circle into two equal groups.
pixel_bearing_quadrant <- function(df){
  df |>
    mutate(bearing_deg = (geosphere::bearing(cbind(lon_Hyp, lat_Hyp), cbind(pixel_lon, pixel_lat)) + 360) %% 360,
           quadrant = case_when(
             pixel_lat >= lat_Hyp & pixel_lon <  lon_Hyp ~ "NW",
             pixel_lat >= lat_Hyp & pixel_lon >= lon_Hyp ~ "NE",
             pixel_lat <  lat_Hyp & pixel_lon <  lon_Hyp ~ "SW",
             TRUE ~ "SE"))
}


# Station coordinate sanity check --------------------------------------------
# One-off diagnostic: does the fixed meta/station_in_situ.csv THFR coordinate
# agree with the actual per-scan HYPERNETS GPS logged in the db? A large,
# consistent offset here could itself be feeding into Bug 9 (THFR PACE pixel
# extraction offset, see manuscript/upstream-data-bugs.md) if the same fixed
# coordinate is what Hypernets_matchups queries against upstream.

station_csv <- read_csv("meta/station_in_situ.csv", show_col_types = FALSE) |>
  dplyr::filter(site == "THFR")

hyp_positions <- db_matchup_pixels(db_path, "S3A") |>
  dplyr::select(match_date, lon_Hyp, lat_Hyp) |>
  distinct()

station_offset <- hyp_positions |>
  mutate(dist_km = round(distHaversine(cbind(lon_Hyp, lat_Hyp), cbind(station_csv$lon, station_csv$lat)) / 1000, 3),
         bearing_from_csv_deg = round((geosphere::bearing(cbind(station_csv$lon, station_csv$lat), cbind(lon_Hyp, lat_Hyp)) + 360) %% 360, 1))

print(station_offset)
# All rows should show the same dist_km/bearing (one fixed station, GPS noise
# aside). If so, meta/station_in_situ.csv is just stale/rounded and the
# db's per-scan position is the one to trust (as used throughout this script).

# Map of the two coordinate sources against a real basemap tile. The GADM
# shapefile in meta/FRANCE_shapefile/ is a national-boundary (level 0) outline
# only -- it does not resolve the Etang de Thau lagoon as water, so it can't
# show local context here. OpenStreetMap tiles (via maptiles) are used instead;
# Esri's satellite imagery has a coverage gap over this part of the lagoon and
# returns broken placeholder tiles for it.
db_centroid <- hyp_positions |> summarise(lon = mean(lon_Hyp), lat = mean(lat_Hyp))

pad <- 0.05
map_bbox <- st_bbox(c(xmin = min(c(hyp_positions$lon_Hyp, station_csv$lon)) - pad,
                       xmax = max(c(hyp_positions$lon_Hyp, station_csv$lon)) + pad,
                       ymin = min(c(hyp_positions$lat_Hyp, station_csv$lat)) - pad,
                       ymax = max(c(hyp_positions$lat_Hyp, station_csv$lat)) + pad),
                     crs = 4326) |> st_as_sfc() |> st_sf()
thau_tile <- get_tiles(map_bbox, provider = "Esri.WorldImagery", crop = TRUE)

pts_sf <- bind_rows(
  station_csv |> transmute(lon, lat, source = "CSV (station_in_situ.csv)"),
  hyp_positions |> transmute(lon = lon_Hyp, lat = lat_Hyp, source = "DB (per-scan GPS)")
) |> st_as_sf(coords = c("lon", "lat"), crs = 4326)

p_station_map <- ggplot() +
  geom_spatraster_rgb(data = thau_tile) +
  geom_segment(aes(x = station_csv$lon, y = station_csv$lat, xend = db_centroid$lon, yend = db_centroid$lat),
               colour = "white", linetype = "dashed", linewidth = 0.6) +
  geom_sf(data = pts_sf, aes(colour = source, shape = source), size = 3, stroke = 1.2, alpha = 0.7) +
  scale_colour_manual(values = c("CSV (station_in_situ.csv)" = "red", "DB (per-scan GPS)" = "yellow")) +
  scale_shape_manual(values = c("CSV (station_in_situ.csv)" = 17, "DB (per-scan GPS)" = 16)) +
  coord_sf(crs = st_crs(thau_tile), expand = FALSE) +
  labs(title = paste0("THFR station position -- CSV vs. DB GPS (offset ~", round(mean(station_offset$dist_km), 2), " km)"),
       x = NULL, y = NULL, colour = NULL, shape = NULL) +
  theme_bw()
ggsave(file.path(out_dir, "thau_station_position_map.png"), p_station_map, width = 8, height = 7)


# All per-pixel geographic positions, per sensor ------------------------------
# Every individual matchup's raw pixel_lon/pixel_lat per pixel_pos (col,row
# within the 5x5 box), rather than one average point per cell, coloured by
# pixel_pos on a discrete colour scale (ordered numerically by col then row,
# not alphabetically, so nearby grid cells get nearby colours). This shows
# whether a sensor's pixel grid sits in a consistent spatial location/
# orientation across passes (S3A/OLCI: a clean, tight, smoothly-graded
# rotated 5x5 cloud) or is scattered/overlapping because swath geometry
# rotates pass-to-pass (SNPP/VIIRS and PACE both show heavy inter-cell mixing).
# Esri.WorldImagery is used so real water-surface detail (e.g. the oyster-bed
# rows visible NW/W of the station) is visible under the points; note Esri has
# a confirmed, persistent tile gap (renders black) over the western half of
# the lagoon, which will show through for any sensor whose scatter reaches
# that far, accepted as a known limitation rather than switching provider.
# Facetting this on one shared basemap isn't possible: coord_sf() doesn't
# support facet_wrap(scales = "free") in this ggplot2 version, and pixel
# spacing differs ~5x between sensors, so each sensor gets its own tile/bbox/
# ggsave() call instead, matching the other per-sensor figures in this script.

for(sY in sensor_Y_thfr){
  pixel_positions_all <- tryCatch({
    db_matchup_pixels(db_path, sY) |>
      distinct(matchup_id, pixel_pos, pixel_lon, pixel_lat, lon_Hyp, lat_Hyp)
  }, warning = function(w){ message(conditionMessage(w)); tibble() })
  if(nrow(pixel_positions_all) == 0) next

  # Numeric col/row ordering (not alphabetical) drives the discrete legend --
  # alphabetical order would put "-1,*" before "-2,*" (ASCII sorts 1 before 2).
  pixel_pos_key <- pixel_positions_all |>
    distinct(pixel_pos) |>
    separate(pixel_pos, into = c("col", "row"), sep = ",", convert = TRUE, remove = FALSE) |>
    arrange(col, row)

  pixel_positions_all <- pixel_positions_all |>
    mutate(pixel_pos = factor(pixel_pos, levels = pixel_pos_key$pixel_pos))

  station_pos <- pixel_positions_all |> summarise(lon_Hyp = mean(lon_Hyp), lat_Hyp = mean(lat_Hyp))

  pad3 <- 0.003
  map_bbox3 <- st_bbox(c(xmin = min(c(pixel_positions_all$pixel_lon, station_pos$lon_Hyp)) - pad3,
                          xmax = max(c(pixel_positions_all$pixel_lon, station_pos$lon_Hyp)) + pad3,
                          ymin = min(c(pixel_positions_all$pixel_lat, station_pos$lat_Hyp)) - pad3,
                          ymax = max(c(pixel_positions_all$pixel_lat, station_pos$lat_Hyp)) + pad3),
                        crs = 4326) |> st_as_sfc() |> st_sf()
  pixel_pos_tile <- get_tiles(map_bbox3, provider = "Esri.WorldImagery", crop = TRUE)

  p_pixel_all_pos <- ggplot() +
    geom_spatraster_rgb(data = pixel_pos_tile) +
    geom_point(data = pixel_positions_all, aes(x = pixel_lon, y = pixel_lat, colour = pixel_pos),
               size = 1.2, alpha = 0.5) +
    geom_point(data = station_pos, aes(x = lon_Hyp, y = lat_Hyp), colour = "red", shape = 4, size = 3, stroke = 1.5) +
    scale_colour_viridis_d(name = "pixel_pos", guide = guide_legend(ncol = 5, override.aes = list(size = 3, alpha = 1))) +
    coord_sf(crs = st_crs(pixel_pos_tile), expand = FALSE) +
    labs(title = paste0("THFR ", sY, " -- all satellite pixel positions, coloured by grid cell"),
         x = NULL, y = NULL) +
    theme_bw()
  ggsave(file.path(out_dir, paste0("pixel_position_all_", sY, ".png")), p_pixel_all_pos, width = 9.5, height = 7)
}


# Hand-drawn "clean water" polygon, N/NE of station --------------------------
# A polygon that traces the open water north/north-east of THFR,
# deliberately excluding land and the visible oyster-bed rows.
# Traced by eye off Esri.WorldImagery at zoom 17
# using a fine (0.001 deg) coordinate grid overlay for calibration.

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

pad4 <- 0.012
map_bbox4 <- st_bbox(c(xmin = db_centroid$lon - pad4, xmax = db_centroid$lon + pad4,
                        ymin = db_centroid$lat - pad4 * 0.6, ymax = db_centroid$lat + pad4),
                      crs = 4326) |> st_as_sfc() |> st_sf()
clean_water_tile <- get_tiles(map_bbox4, provider = "Esri.WorldImagery", crop = TRUE, zoom = 17)

p_clean_water <- ggplot() +
  geom_spatraster_rgb(data = clean_water_tile) +
  geom_sf(data = clean_water_sf, fill = "cyan", alpha = 0.25, colour = "cyan", linewidth = 1.2) +
  geom_point(aes(x = db_centroid$lon, y = db_centroid$lat), colour = "red", shape = 4, size = 4, stroke = 2) +
  coord_sf(crs = st_crs(clean_water_tile), expand = FALSE) +
  labs(title = "THFR clean-water polygon N/NE of station", x = NULL, y = NULL) +
  theme_bw()
ggsave(file.path(out_dir, "clean_water_polygon_NE.png"), p_clean_water, width = 9, height = 7.5)


# Single-matchup inspection (EDIT ME) ----------------------------------------
# Re-run this section for different sensors/dates to look at one matchup's
# raw 5x5 pixel grid by hand.

# Load all data for one sensor
sensor_Y_pick <- "S3A"                 # EDIT ME
df_pick_all <- db_matchup_pixels(db_path, sensor_Y_pick) |> pixel_bearing_quadrant()

# Pick one date for testing
unique(df_pick_all$match_date)
date_pick <- as.Date("2025-11-30")#NA    # EDIT ME or NA to just take the first available matchup
matchup_id_pick <- if(is.na(date_pick)){
  df_pick_all$matchup_id[1]
} else {
  df_pick_all |> dplyr::filter(match_date == date_pick) |> pull(matchup_id) |> head(1)
}

# Check for differences
df_pick <- df_pick_all |> dplyr::filter(matchup_id == matchup_id_pick) |>
  mutate(deviation = .data[[sensor_Y_pick]] - Hyp)

# Spatial map of one wavelength's pixel box: raw satellite RHOW (left) and
# deviation from the paired HYPERNETS value (right), with the HYPERNETS
# position marked and the quadrant split drawn in
wavelength_pick <- df_pick$wavelength[which.min(abs(df_pick$wavelength - 560))] # nearest to 560 nm (green)

p_raw <- df_pick |> dplyr::filter(wavelength == wavelength_pick) |>
  ggplot(aes(x = pixel_lon, y = pixel_lat)) +
  geom_point(aes(colour = .data[[sensor_Y_pick]]), size = 8) +
  geom_point(data = df_pick[1,], aes(x = lon_Hyp, y = lat_Hyp), shape = 4, size = 5, stroke = 2) +
  geom_vline(xintercept = df_pick$lon_Hyp[1], linetype = "dashed") +
  geom_hline(yintercept = df_pick$lat_Hyp[1], linetype = "dashed") +
  scale_colour_viridis_c(name = paste0(sensor_Y_pick, "\nRHOW")) +
  coord_fixed() +
  labs(title = paste0(sensor_Y_pick, " raw pixel RHOW -- ", df_pick$match_date[1], " -- ", wavelength_pick, " nm"),
       x = "Longitude", y = "Latitude") +
  theme_bw()

p_dev <- df_pick |> dplyr::filter(wavelength == wavelength_pick) |>
  ggplot(aes(x = pixel_lon, y = pixel_lat)) +
  geom_point(aes(colour = deviation), size = 8) +
  geom_point(data = df_pick[1,], aes(x = lon_Hyp, y = lat_Hyp), shape = 4, size = 5, stroke = 2) +
  geom_vline(xintercept = df_pick$lon_Hyp[1], linetype = "dashed") +
  geom_hline(yintercept = df_pick$lat_Hyp[1], linetype = "dashed") +
  scale_colour_gradient2(name = "Sat - Hyp", low = "blue", mid = "white", high = "red") +
  coord_fixed() +
  labs(title = paste0(sensor_Y_pick, " deviation from HYPERNETS -- ", df_pick$match_date[1], " -- ", wavelength_pick, " nm"),
       x = "Longitude", y = "Latitude") +
  theme_bw()

p_raw + p_dev
ggsave(file.path(out_dir, paste0("single_matchup_", sensor_Y_pick, "_", df_pick$match_date[1], ".png")), p_raw + p_dev, width = 12, height = 6)


# Aggregate diagnostic across all THFR matchups ------------------------------
# Builds the full per-pixel deviation table across every sensor family, then
# summarises Sat-Hyp deviation and raw Sat RHOW by quadrant.
# See meta/pixel_explore_output/ for the written CSVs/figures/summary.

# pixel_all <- map_dfr(sensor_Y_thfr, function(sY){
#   tryCatch({
#     db_matchup_pixels(db_path, sY) |>
#       pixel_bearing_quadrant() |>
#       mutate(sat_value = .data[[sY]], deviation = .data[[sY]] - Hyp) |>
#       dplyr::select(sensor_Y, matchup_id, match_date, wavelength, pixel_pos,
#                     bearing_deg, quadrant, dist_km, diff_time_min, Hyp, sat_value, deviation)
#   }, warning = function(w){ message(conditionMessage(w)); tibble() })
# })
# write_csv(pixel_all, file.path(out_dir, "pixel_deviation_all.csv"))
pixel_all <- read_csv("meta/pixel_explore_output/pixel_deviation_all.csv")

# Data-quality guard: a handful of HYPERNETS scans in thfr_2025.db carry the
# raw NetCDF fill-value sentinel (9.9692099683868690e36) in Hyp instead of NA
# Confirmed on 10 (matchup_id, sensor_Y) combinations across 2025-07-17,
# 2025-07-25, 2025-10-27 and 2025-11-03 (2026-08-27). Confirmed NOT present in
# the corresponding exported .csv files used by the production pipeline (e.g.
# RHOW_HYPERNETS_vs_AQUA/AQUA_20251103T142001_..._RHOW_R.csv has normal
# values). This is specific to the newer .db ingestion (2026-08-26), 
# not a manuscript stats bug, but must be filtered here before any mean-based 
# summary (a single 1e36 value swamps a mean). 
# Better to flag upstream at the soruce rather than fixing in place here.
n_before <- nrow(pixel_all)
pixel_all <- pixel_all |> dplyr::filter(abs(Hyp) < 10, abs(sat_value) < 10)
message("Dropped ", n_before - nrow(pixel_all), " of ", n_before,
        " rows with a corrupted (fill-value) Hyp/sat_value before summarising")

# Quadrant x sensor x wavelength summary
# quadrant_summary <- pixel_all |>
#   summarise(n = n(),
#             dev_mean = mean(deviation, na.rm = TRUE),
#             dev_median = median(deviation, na.rm = TRUE),
#             dev_iqr = IQR(deviation, na.rm = TRUE),
#             sat_mean = mean(sat_value, na.rm = TRUE),
#             .by = c(sensor_Y, wavelength, quadrant))
# write_csv(quadrant_summary, file.path(out_dir, "quadrant_summary.csv"))
quadrant_summary <- read_csv("meta/pixel_explore_output/quadrant_summary.csv")

# All 6 pairwise quadrant comparisons (NE-NW, NE-SE, NE-SW, NW-SE, NW-SW,
# SE-SW), per sensor x wavelength: average deviation within each quadrant per
# matchup first (so the test isn't inflated by unequal pixel counts per
# quadrant -- see the quadrant_summary pixel-count table, which is very
# unbalanced), then pair on matchup_id and run a paired Wilcoxon signed-rank
# test for every pair. Long format: one row per (sensor_Y, wavelength,
# quadrant_a, quadrant_b).
quadrant_wide <- pixel_all |>
  summarise(dev_mean = mean(deviation, na.rm = TRUE), .by = c(sensor_Y, wavelength, matchup_id, quadrant)) |>
  pivot_wider(names_from = quadrant, values_from = dev_mean)

quadrant_pair_list <- combn(c("NE", "NW", "SE", "SW"), 2, simplify = FALSE)

# quadrant_test_summary <- map_dfr(quadrant_pair_list, function(qp){
#   qa <- qp[1]; qb <- qp[2]
#   quadrant_wide |>
#     dplyr::rename(a = all_of(qa), b = all_of(qb)) |>
#     dplyr::select(sensor_Y, wavelength, matchup_id, a, b) |>
#     dplyr::filter(!is.na(a), !is.na(b)) |>
#     reframe({
#       if(n() >= 5){
#         wt <- suppressWarnings(wilcox.test(a, b, paired = TRUE))
#         tibble(n_matchups = n(), mean_a = mean(a), mean_b = mean(b),
#                diff_mean = mean(a - b), p_value = wt$p.value)
#       } else {
#         tibble(n_matchups = n(), mean_a = mean(a), mean_b = mean(b),
#                diff_mean = mean(a - b), p_value = NA_real_)
#       }
#     }, .by = c(sensor_Y, wavelength)) |>
#     mutate(quadrant_a = qa, quadrant_b = qb, .before = 1)
# })
# write_csv(quadrant_test_summary, file.path(out_dir, "quadrant_test_summary.csv"))
quadrant_test_summary <- read_csv("meta/pixel_explore_output/quadrant_test_summary.csv")

# Continuous bearing-vs-deviation check per sensor (pooled across wavelength,
# since the directional signal (if real) should show up regardless of band)
# bearing_cor <- pixel_all |>
#   reframe({
#     ct <- suppressWarnings(cor.test(bearing_deg, deviation, method = "spearman"))
#     tibble(n = n(), rho = unname(ct$estimate), p_value = ct$p.value)
#   }, .by = sensor_Y)
# write_csv(bearing_cor, file.path(out_dir, "bearing_correlation.csv"))
bearing_cor <- read_csv("meta/pixel_explore_output/bearing_correlation.csv")


# Figures ---------------------------------------------------------------------

pace_nm_breaks <- c(350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1050)
pace_nm_labels <- c("351-400", "401-450", "451-500", "501-550", "551-600", "601-650", "651-700", "701-750", "751-800", "801-900", "901-1050")

for(sY in unique(pixel_all$sensor_Y)){
  df_s <- pixel_all |> dplyr::filter(sensor_Y == sY)

  # PACE is hyperspectral (801 discrete bands), facetting on raw wavelength
  # would produce 801 near-empty panels, so bin into the same broad bands
  # colour_nm_func() already uses for PACE elsewhere in the codebase
  df_s <- df_s |>
    mutate(wavelength_facet = if(sY == "PACE"){
      as.character(cut(wavelength, breaks = pace_nm_breaks, labels = pace_nm_labels, include.lowest = TRUE))
    } else {
      as.character(wavelength)
    })

  # Boxplots per sector
  p_box <- df_s |>
    ggplot(aes(x = quadrant, y = deviation, fill = quadrant)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~wavelength_facet, scales = "free_y") +
    labs(title = paste0("THFR ", sY, " -- pixel deviation from HYPERNETS by quadrant"),
         x = "Quadrant (relative to HYPERNETS station)", y = "Sat pixel RHOW - Hyp RHOW") +
    theme_bw() +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, paste0("deviation_boxplot_", sY, ".png")), p_box, width = 10, height = 8)

  # Boxplots per exact pixel position (finer-grained than quadrant, above)
  pixel_pos_levels <- df_s |>
    distinct(pixel_pos) |>
    separate(pixel_pos, into = c("col", "row"), sep = ",", convert = TRUE, remove = FALSE) |>
    arrange(col, row) |>
    pull(pixel_pos)

  p_box_pixel <- df_s |>
    mutate(pixel_pos = factor(pixel_pos, levels = pixel_pos_levels)) |>
    ggplot(aes(x = pixel_pos, y = deviation, fill = pixel_pos)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~wavelength_facet, scales = "free_y") +
    labs(title = paste0("THFR ", sY, " -- pixel deviation from HYPERNETS by pixel position"),
         x = "Pixel position (col,row within 5x5 box)", y = "Sat pixel RHOW - Hyp RHOW") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, paste0("deviation_boxplot_pixelpos_", sY, ".png")), p_box_pixel, width = 12, height = 9)

  # Line plot of all degrees 0-360
  # NB: subsampled for the large sensors (PACE especially, 2.85M pooled rows)
  # Plotting every point isn't needed to see the pattern and is slow/heavy
  set.seed(1)
  df_plot <- if(nrow(df_s) > 100000) slice_sample(df_s, n = 100000) else df_s

  p_bearing <- df_plot |>
    ggplot(aes(x = bearing_deg, y = deviation)) +
    geom_point(alpha = 0.15, size = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(breaks = c(0, 90, 180, 270, 360), limits = c(0, 360)) +
    labs(title = paste0("THFR ", sY, " -- deviation vs. bearing from station (all wavelengths pooled)",
                         if(nrow(df_s) > 100000) paste0("\n(", nrow(df_plot), "-pt random subsample of ", nrow(df_s), ")") else ""),
         x = "Bearing from HYPERNETS station (deg, 0 = N, 90 = E)", y = "Sat pixel RHOW - Hyp RHOW") +
    theme_bw()
  ggsave(file.path(out_dir, paste0("bearing_scatter_", sY, ".png")), p_bearing, width = 8, height = 6)
}


# Visualise quadrant_test_summary ---------------------------------------------
# All 6 pairwise paired-Wilcoxon comparisons, one coloured line per pair,
# faceted by sensor. Filled points = significant (p < 0.05), open = not.

p_quadrant_test <- quadrant_test_summary |>
  mutate(pair_label = paste0(quadrant_a, " - ", quadrant_b),
         significant = if_else(!is.na(p_value) & p_value < 0.05, "p < 0.05", "n.s.")) |>
  ggplot(aes(x = wavelength, y = diff_mean, colour = pair_label)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point(aes(shape = significant), size = 1.6) +
  scale_shape_manual(values = c("p < 0.05" = 16, "n.s." = 1)) +
  facet_wrap(~sensor_Y, scales = "free") +
  labs(title = "THFR -- pairwise quadrant deviation tests, by sensor",
       subtitle = "diff_mean = mean(quadrant A deviation) - mean(quadrant B deviation), per matchup (paired Wilcoxon); open point = p >= 0.05",
       x = "Wavelength (nm)", y = "diff_mean  (A - B)", colour = "Pair (A - B)", shape = NULL) +
  theme_bw()
ggsave(file.path(out_dir, "quadrant_test_summary_plot.png"), p_quadrant_test, width = 12, height = 8)

