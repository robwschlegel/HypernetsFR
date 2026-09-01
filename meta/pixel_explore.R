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

# IGN Geoportail aerial orthophoto tile provider
# Esri.WorldImagery has a confirmed, persistent black tile-coverage
# gap over the western half of the lagoon, which cuts through the 
# oyster-bed field. IGN's national orthophoto layer has full gap-free 
# coverage of the whole lagoon, no API key required, so it's used 
# specifically for the oyster-bed section below where
# seeing the field's full extent matters.
ign_ortho <- create_provider(
  name = "IGN.Orthophotos",
  url = "https://data.geopf.fr/wmts?LAYER=ORTHOIMAGERY.ORTHOPHOTOS&EXCEPTIONS=text/xml&FORMAT=image/jpeg&SERVICE=WMTS&VERSION=1.0.0&REQUEST=GetTile&STYLE=normal&TILEMATRIXSET=PM&TILEMATRIX={z}&TILEROW={y}&TILECOL={x}",
  citation = "IGN-F/Geoportail"
)

# IGN's WMTS endpoint enforces a strict per-second rate limit. maptiles's
# internal downloader can trip it on any single get_tiles() call that needs
# enough tiles. forceDownload = TRUE matters here, not just the retry loop:
# get_tiles()'s default cachedir is tempdir(), which is stable for the life
# of an R session (unlike a fresh `Rscript` call, which always gets a new
# tempdir()) -- so in a long-running interactive session, a tile that failed
# or came back truncated on an earlier attempt can get cached in a corrupt
# state, and later calls may silently read that same bad cached tile and
# return an incomplete mosaic *without raising an R error at all*. Plain
# error-based retry can't catch that (nothing errors), which is exactly why
# a fresh `Rscript` run of this script can succeed while the same code
# re-run inside an already-used interactive session keeps reproducing the
# same cut-off tile. Forcing a fresh download every attempt avoids trusting
# that cache at all.
get_tiles_retry <- function(..., max_tries = 4, pause_sec = 3){
  for(i in seq_len(max_tries)){
    result <- tryCatch(get_tiles(..., forceDownload = TRUE), error = function(e) e)
    if(!inherits(result, "error")) return(result)
    message("get_tiles() attempt ", i, " failed (", conditionMessage(result), "), retrying...")
    Sys.sleep(pause_sec)
  }
  stop("get_tiles() failed after ", max_tries, " attempts")
}

# PACE is hyperspectral (801 discrete bands)
# Facetting a figure on raw wavelength produces ~801 near-empty panels, 
# so anywhere PACE needs faceting, its wavelengths get binned into these 
# same broad bands colour_nm_func() (code/0_functions.R) already uses for 
# PACE elsewhere in the codebase. Defined once here in Setup since multiple 
# sections below (both the Figures section and the grid-deviation pipeline) need it.
pace_nm_breaks <- c(350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1050)
pace_nm_labels <- c("351-400", "401-450", "451-500", "501-550", "551-600", 
                     "601-650", "651-700", "701-750", "751-800", "801-900", "901-1050")

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

# Map of the two coordinate sources against a real basemap tile.
db_centroid <- hyp_positions |> summarise(lon = mean(lon_Hyp), lat = mean(lat_Hyp))

pad <- 0.05
map_bbox <- st_bbox(c(xmin = min(c(hyp_positions$lon_Hyp, station_csv$lon)) - pad,
                       xmax = max(c(hyp_positions$lon_Hyp, station_csv$lon)) + pad,
                       ymin = min(c(hyp_positions$lat_Hyp, station_csv$lat)) - pad,
                       ymax = max(c(hyp_positions$lat_Hyp, station_csv$lat)) + pad),
                     crs = 4326) |> st_as_sfc() |> st_sf()
thau_tile <- get_tiles(map_bbox, provider = ign_ortho, crop = TRUE)

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
  labs(title = paste0("THFR station position: CSV vs. DB GPS (offset ~", round(mean(station_offset$dist_km), 2), " km)"),
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

  # Numeric col/row ordering (not alphabetical) drives the discrete legend
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
  pixel_pos_tile <- get_tiles(map_bbox3, provider = ign_ortho, crop = TRUE)

  p_pixel_all_pos <- ggplot() +
    geom_spatraster_rgb(data = pixel_pos_tile) +
    geom_point(data = pixel_positions_all, aes(x = pixel_lon, y = pixel_lat, colour = pixel_pos),
               size = 1.2, alpha = 0.5) +
    geom_point(data = station_pos, aes(x = lon_Hyp, y = lat_Hyp), colour = "red", shape = 4, size = 3, stroke = 1.5) +
    scale_colour_viridis_d(name = "pixel_pos", guide = guide_legend(ncol = 5, override.aes = list(size = 3, alpha = 1))) +
    coord_sf(crs = st_crs(pixel_pos_tile), expand = FALSE) +
    labs(title = paste0("THFR ", sY, ": all satellite pixel positions, coloured by grid cell"),
         x = NULL, y = NULL) +
    theme_bw()
  ggsave(file.path(out_dir, paste0("pixel_position_all_", sY, ".png")), p_pixel_all_pos, width = 9.5, height = 7)
}


# Hand-drawn "clean water" polygon, N/NE of station --------------------------
# clean_water_polygon / clean_water_sf are now defined in code/0_functions.R (single-sourced so
# this exploratory script and the THFR_poly pipeline site db_export_matchups_poly() /
# pixel_filter_clean_water() never drift apart), already available here via source() above.

pad4 <- 0.012
map_bbox4 <- st_bbox(c(xmin = db_centroid$lon - pad4, xmax = db_centroid$lon + pad4,
                        ymin = db_centroid$lat - pad4 * 0.6, ymax = db_centroid$lat + pad4),
                      crs = 4326) |> st_as_sfc() |> st_sf()
clean_water_tile <- get_tiles(map_bbox4, provider = ign_ortho, crop = TRUE, zoom = 17)

p_clean_water <- ggplot() +
  geom_spatraster_rgb(data = clean_water_tile) +
  geom_sf(data = clean_water_sf, fill = "cyan", alpha = 0.25, colour = "cyan", linewidth = 1.2) +
  geom_point(aes(x = db_centroid$lon, y = db_centroid$lat), colour = "red", shape = 4, size = 4, stroke = 2) +
  coord_sf(crs = st_crs(clean_water_tile), expand = FALSE) +
  labs(title = "THFR clean-water polygon N/NE of station", x = NULL, y = NULL) +
  theme_bw()
ggsave(file.path(out_dir, "clean_water_polygon_NE.png"), p_clean_water, width = 9, height = 9)


# Hand-drawn oyster-bed polygons, W/SW/NW of station -------------------------
# First-guess polygons tracing the visible oyster-farming rows: one large
# contiguous field W/SW/NW of the station, plus one small isolated cluster
# roughly NNE of the station. Traced by eye from ign_ortho at zoom 16-17 
# using 0.001-0.002 deg coordinate grid overlays for calibration.

main_field <- tribble(
  ~lon,   ~lat,
  3.6655, 43.4340,
  3.6645, 43.4355,
  3.6625, 43.4380,
  3.6600, 43.4400,
  3.6580, 43.4420,
  3.6550, 43.4445,
  3.6525, 43.4465,
  3.6400, 43.4440,
  3.6300, 43.4360,
  3.6180, 43.4290,
  3.6170, 43.4280,
  3.6180, 43.4230,
  3.6200, 43.4190,
  3.6240, 43.4165,
  3.6410, 43.4255,
  3.6530, 43.4275,
  3.6655, 43.4340  # closes the ring
)
cluster_nne <- tribble(
  ~lon,   ~lat,
  3.6660, 43.4390,
  3.6660, 43.4420,
  3.6690, 43.4420,
  3.6690, 43.4390,
  3.6660, 43.4390  # closes the ring
)

as_cluster_sf <- function(coords, label){
  st_sf(cluster = label, geometry = st_sfc(st_polygon(list(as.matrix(coords))), crs = 4326))
}
oyster_bed_sf <- bind_rows(
  as_cluster_sf(main_field, "main_field"),
  as_cluster_sf(cluster_nne, "cluster_NNE")
)

map_bbox5 <- st_bbox(c(xmin = 3.605, xmax = 3.71, ymin = 43.412, ymax = 43.448),
                      crs = 4326) |> st_as_sfc() |> st_sf()
oyster_bed_tile <- get_tiles_retry(map_bbox5, provider = ign_ortho, crop = TRUE, zoom = 16)

p_oyster_beds <- ggplot() +
  geom_spatraster_rgb(data = oyster_bed_tile) +
  geom_sf(data = oyster_bed_sf, aes(fill = cluster), alpha = 0.3, colour = "orange", linewidth = 1) +
  geom_point(aes(x = db_centroid$lon, y = db_centroid$lat), colour = "red", shape = 4, size = 4, stroke = 2) +
  coord_sf(crs = st_crs(oyster_bed_tile), expand = FALSE) +
  labs(title = "THFR oyster-bed polygons, W/SW/NW (+ 1 small NNE cluster)", x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "oyster_bed_polygons.png"), p_oyster_beds, width = 12, height = 7)


# Point-in-oyster-bed spatial test, S3A validation ----------------------------
# Core building block for the exclusion-logic plan: a point-in-polygon test
# against oyster_bed_sf (main_field + cluster_nne, both above), so any pixel's
# real lon/lat can be tagged in/out of the mask. 
# st_union() first so a pixel counts as "in" if it falls in EITHER
# polygon, without double-counting anything that might straddle both.

oyster_bed_union <- st_union(oyster_bed_sf)

pixel_in_oyster_bed <- function(lon, lat){
  pts <- st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326)
  lengths(st_intersects(pts, oyster_bed_union)) > 0
}

s3a_pixels <- db_matchup_pixels(db_path, "S3A") |>
  distinct(matchup_id, pixel_pos, pixel_lon, pixel_lat) |>
  mutate(in_oyster_bed = pixel_in_oyster_bed(pixel_lon, pixel_lat))

# Sanity check: fraction of matchups where each pixel_pos falls inside the
# mask, since S3A's grid is spatially stable, this should read close to 0
# or 1 per pixel_pos (not some in-between value), tracking a NW-SE gradient
# out from the station.
pixel_pos_summary <- s3a_pixels |>
  summarise(frac_in_bed = round(mean(in_oyster_bed), 2), n = n(), .by = pixel_pos) |>
  arrange(desc(frac_in_bed))
print(pixel_pos_summary)

# Visual validation: colour every individual pixel by the test result and
# overlay the mask outline, points should switch from cyan (out) to red
# (in) right at the orange polygon boundary, with no visible mismatches.
map_bbox6 <- st_bbox(c(xmin = min(s3a_pixels$pixel_lon) - 0.002, xmax = max(s3a_pixels$pixel_lon) + 0.002,
                        ymin = min(s3a_pixels$pixel_lat) - 0.002, ymax = max(s3a_pixels$pixel_lat) + 0.002),
                      crs = 4326) |> st_as_sfc() |> st_sf()
test_tile <- get_tiles(map_bbox6, provider = ign_ortho, crop = TRUE, zoom = 17)

p_test <- ggplot() +
  geom_spatraster_rgb(data = test_tile) +
  geom_sf(data = oyster_bed_sf, fill = NA, colour = "orange", linewidth = 1) +
  geom_point(data = s3a_pixels, aes(x = pixel_lon, y = pixel_lat, colour = in_oyster_bed), size = 1.5, alpha = 0.6) +
  geom_point(aes(x = db_centroid$lon, y = db_centroid$lat), colour = "yellow", shape = 4, size = 4, stroke = 2) +
  scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "deepskyblue"), name = "In oyster bed") +
  coord_sf(crs = st_crs(test_tile), expand = FALSE) +
  labs(title = "S3A pixel-in-oyster-bed test (validation)", x = NULL, y = NULL) +
  theme_bw()
ggsave(file.path(out_dir, "pixel_in_oyster_bed_test_S3A.png"), p_test, width = 9.5, height = 7.5)


# Grid-position deviation / oyster-bed contribution test ---------------------
# For every (matchup, wavelength): compute the mean/SD/CV of the 3x3 inner
# grid (|col|<=1 & |row|<=1, e.g. pixel_pos "1,1") and separately of the full
# 5x5 grid (e.g. pixel_pos "2,2", only in the 5x5), then rank every pixel by
# its deviation from the 5x5 (full-grid) mean and flag whether it falls inside
# the oyster-bed mask (pixel_in_oyster_bed(), above).
#
# Purpose: test whether the pixels contributing most to a matchup's variance
# are consistently the ones sitting over the oyster beds (supports the
# spatial-contamination hypothesis) or whether high-deviation pixels show up
# regardless of oyster-bed position (points instead to a sensor/waveband
# measurement issue in the Etang de Thau, not spatial contamination).
#
# Note on CV: at THFR, RHOW values are small and can sit close to zero, so
# CV% = 100*sd/mean is numerically unstable there (seen as high as 800%+ for
# SNPP below) -- it's reported because it was asked for, but abs_deviation
# (an absolute, not relative, measure) is the more robust ranking metric and
# is what drives `rank` and the in/out-of-bed comparison below.
#
# Implemented as one vectorised pipeline (summarise(.by=)/mutate(.by=), not
# group_modify()) since group_modify() over every (matchup_id, wavelength)
# combination -- SNPP alone has ~1880 matchups x 5 wavelengths -- took minutes;
# this runs in ~1 second for the same data.

pixel_pos_grid_class <- function(pixel_pos){
  cr <- str_split_fixed(pixel_pos, ",", 2)
  col <- as.integer(cr[, 1]); row <- as.integer(cr[, 2])
  if_else(abs(col) <= 1 & abs(row) <= 1, "3x3", "5x5_only")
}

# df_all_pixels: db_matchup_pixels() output for one sensor_Y (any number of
# matchups/wavelengths). sensor_col: the sat_value column name (== sensor_Y).
# Returns one row per pixel, with grid_class, in_oyster_bed, both grid sizes'
# mean/sd/cv, deviation from the 5x5 mean (deviation/abs_deviation/rank) and
# from the 3x3 mean (deviation_3x3/abs_deviation_3x3/rank_3x3 -- NA for
# 5x5_only rows, since they aren't members of the 3x3 grid), each rank
# computed independently within its own (matchup_id, wavelength[, grid_class])
# group.
pixel_deviation_pipeline <- function(df_all_pixels, sensor_col){
  # Data-quality guard (same issue/fix as the Figures section further down):
  # a handful of HYPERNETS scans in thfr_2025.db carry the raw NetCDF
  # fill-value sentinel (9.9692099683868690e36) in Hyp instead of NA, and
  # occasionally in the satellite value too. A single such row swamps any
  # mean-based summary (deviation_hyp especially, since it's a direct Hyp -
  # sat difference) -- must be dropped before any of the stats below.
  df_all_pixels <- df_all_pixels |> dplyr::filter(abs(Hyp) < 10, abs(.data[[sensor_col]]) < 10)

  df_base <- df_all_pixels |>
    mutate(grid_class = pixel_pos_grid_class(pixel_pos),
           in_oyster_bed = pixel_in_oyster_bed(pixel_lon, pixel_lat))

  stats_5x5 <- df_base |>
    summarise(grid5x5_mean = mean(.data[[sensor_col]], na.rm = TRUE),
              grid5x5_sd = sd(.data[[sensor_col]], na.rm = TRUE),
              .by = c(matchup_id, wavelength)) |>
    mutate(grid5x5_cv_pct = round(100 * grid5x5_sd / grid5x5_mean, 1))

  stats_3x3 <- df_base |> dplyr::filter(grid_class == "3x3") |>
    summarise(grid3x3_mean = mean(.data[[sensor_col]], na.rm = TRUE),
              grid3x3_sd = sd(.data[[sensor_col]], na.rm = TRUE),
              .by = c(matchup_id, wavelength)) |>
    mutate(grid3x3_cv_pct = round(100 * grid3x3_sd / grid3x3_mean, 1))

  df_base |>
    left_join(stats_5x5, by = c("matchup_id", "wavelength")) |>
    left_join(stats_3x3, by = c("matchup_id", "wavelength")) |>
    mutate(deviation = .data[[sensor_col]] - grid5x5_mean,
           abs_deviation = abs(deviation),
           deviation_3x3 = .data[[sensor_col]] - grid3x3_mean,
           abs_deviation_3x3 = abs(deviation_3x3),
           # In situ minus satellite -- how far this pixel's own value sits
           # from the paired HYPERNETS measurement, i.e. actual matchup
           # quality, rather than internal grid-mean consistency. Hyp is the
           # same single in-situ value for every pixel in a matchup, so this
           # doesn't depend on grid size -- only which pixel population the
           # rank/comparison is computed over does (hence _3x3 variant below,
           # mirroring rank_3x3).
           deviation_hyp = Hyp - .data[[sensor_col]],
           abs_deviation_hyp = abs(deviation_hyp)) |>
    arrange(matchup_id, wavelength, desc(abs_deviation)) |>
    mutate(rank = row_number(), .by = c(matchup_id, wavelength)) |>
    arrange(matchup_id, wavelength, grid_class, desc(abs_deviation_3x3)) |>
    mutate(rank_3x3 = if_else(grid_class == "3x3", row_number(), NA_integer_),
           .by = c(matchup_id, wavelength, grid_class)) |>
    arrange(matchup_id, wavelength, desc(abs_deviation_hyp)) |>
    mutate(rank_hyp = row_number(), .by = c(matchup_id, wavelength)) |>
    arrange(matchup_id, wavelength, grid_class, desc(abs_deviation_hyp)) |>
    mutate(rank_hyp_3x3 = if_else(grid_class == "3x3", row_number(), NA_integer_),
           .by = c(matchup_id, wavelength, grid_class))
}


# Manual per-day test (EDIT ME) -----------------------------------------------
# Pick one sensor + date, see the full ranked pixel table for every wavelength.

sensor_Y_dev <- "SNPP"               # EDIT ME
date_dev <- NA                       # EDIT ME (e.g. as.Date("2025-11-30")) or NA for the first available day

df_dev_all <- db_matchup_pixels(db_path, sensor_Y_dev)
matchup_id_dev <- if(is.na(date_dev)){
  df_dev_all$matchup_id[1]
} else {
  df_dev_all |> dplyr::filter(match_date == date_dev) |> pull(matchup_id) |> head(1)
}

day_deviation <- df_dev_all |>
  dplyr::filter(matchup_id == matchup_id_dev) |>
  pixel_deviation_pipeline(sensor_Y_dev)

print(day_deviation |>
        dplyr::select(wavelength, rank, pixel_pos, grid_class, in_oyster_bed,
                       all_of(sensor_Y_dev), deviation, grid3x3_cv_pct, grid5x5_cv_pct) |>
        arrange(wavelength, rank))

# Spatial map for the same day: every pixel's real position on a real
# (ign_ortho) basemap, oyster-bed mask drawn on top, coloured by deviation
# from HYPERNETS (deviation_hyp = Hyp - sat pixel value, from
# pixel_deviation_pipeline() above), faceted by waveband. PACE gets binned
# into pace_nm_labels (its ~800 native bands would make an unreadable facet
# grid otherwise); every other sensor facets on native wavelength.
day_deviation_facet <- day_deviation |>
  mutate(wavelength_facet = if(sensor_Y_dev == "PACE"){
    as.character(cut(wavelength, breaks = pace_nm_breaks, labels = pace_nm_labels, include.lowest = TRUE))
  } else {
    as.character(wavelength)
  })

pad_day <- 0.005
map_bbox_day <- st_bbox(c(xmin = min(day_deviation_facet$pixel_lon) - pad_day,
                           xmax = max(day_deviation_facet$pixel_lon) + pad_day,
                           ymin = min(day_deviation_facet$pixel_lat) - pad_day,
                           ymax = max(day_deviation_facet$pixel_lat) + pad_day),
                         crs = 4326) |> st_as_sfc() |> st_sf()
day_tile <- get_tiles(map_bbox_day, provider = ign_ortho, crop = TRUE, zoom = 17)

p_day_map <- ggplot() +
  geom_spatraster_rgb(data = day_tile) +
  geom_sf(data = oyster_bed_sf, fill = "cyan", alpha = 0.15, colour = "orange", linewidth = 1) +
  geom_point(data = day_deviation_facet, aes(x = pixel_lon, y = pixel_lat, colour = deviation_hyp), size = 2.5) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Hyp - Sat") +
  facet_wrap(~wavelength_facet) +
  coord_sf(crs = st_crs(day_tile), expand = FALSE) +
  labs(title = paste0(sensor_Y_dev, " -- pixel deviation from HYPERNETS, ", day_deviation_facet$match_date[1]),
       x = NULL, y = NULL) +
  theme_bw()
ggsave(file.path(out_dir, paste0("daily_check_map_", sensor_Y_dev, ".png")), p_day_map, width = 12, height = 10)


# Full-sensor pipeline: every daily matchup, all sensors ---------------------
# Same computation as the SNPP test case, looped over every sensor family
# present for THFR. Writes one CSV + two comparison figures per sensor, plus
# one combined summary/test table across all sensors so the in-bed vs.
# out-of-bed pattern can be compared side by side -- SNPP alone (see above)
# showed no meaningful difference, so the interesting question is whether
# that holds for every sensor or whether some (e.g. PACE, already shown to
# have the most scattered/least stable pixel grid) look different.
#
# PACE is hyperspectral: the underlying grid mean/sd/deviation is still
# computed per exact native wavelength (correct -- binning first would blend
# RHOW across a ~50 nm band into one number), but its comparison figures
# facet on pace_nm_breaks/pace_nm_labels bins instead of ~800 individual
# panels, same convention as the pixel_pos boxplots in the Figures section.

oyster_bed_summary_all <- tibble()
oyster_bed_test_all <- tibble()
oyster_bed_summary_all_3x3 <- tibble()
oyster_bed_test_all_3x3 <- tibble()
oyster_bed_summary_all_hyp_3x3 <- tibble()
oyster_bed_test_all_hyp_3x3 <- tibble()

for(sY in sensor_Y_thfr){
  t0 <- Sys.time()
  sensor_deviation <- tryCatch({
    db_matchup_pixels(db_path, sY) |> pixel_deviation_pipeline(sY)
  }, warning = function(w){ message(conditionMessage(w)); tibble() })
  if(nrow(sensor_deviation) == 0) next
  message(sY, " deviation pipeline: ", round(difftime(Sys.time(), t0, units = "secs"), 1), " sec, ",
          nrow(sensor_deviation), " rows")
  write_csv(sensor_deviation, file.path(out_dir, paste0("pixel_deviation_all_", sY, ".csv")))

  sensor_deviation <- sensor_deviation |>
    mutate(wavelength_facet = if(sY == "PACE"){
      as.character(cut(wavelength, breaks = pace_nm_breaks, labels = pace_nm_labels, include.lowest = TRUE))
    } else {
      as.character(wavelength)
    })

  p_dev <- sensor_deviation |>
    ggplot(aes(x = in_oyster_bed, y = abs_deviation, fill = in_oyster_bed)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~wavelength_facet, scales = "free_y") +
    labs(title = paste0(sY, " -- pixel deviation from 5x5 grid mean, in vs. out of oyster bed"),
         x = "In oyster bed", y = "|Sat pixel RHOW - grid mean|", fill = "In oyster bed") +
    theme_bw()
  ggsave(file.path(out_dir, paste0("deviation_vs_oyster_bed_", sY, ".png")), p_dev, width = 10, height = 8)

  p_rank <- sensor_deviation |>
    ggplot(aes(x = in_oyster_bed, y = rank, fill = in_oyster_bed)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~wavelength_facet) +
    labs(title = paste0(sY, " -- deviation rank (1 = biggest outlier), in vs. out of oyster bed"),
         x = "In oyster bed", y = "Rank within 5x5 grid (1 = furthest from mean)", fill = "In oyster bed") +
    theme_bw()
  ggsave(file.path(out_dir, paste0("deviation_rank_vs_oyster_bed_", sY, ".png")), p_rank, width = 10, height = 8)

  # Quantitative summary: mean |deviation| and mean rank in vs. out of the
  # oyster bed, per wavelength, plus a Wilcoxon rank-sum test (unpaired -- in-
  # and out-of-bed pixels aren't a 1:1 pairing) for whether the in/out
  # difference is significant. Guarded against wavelengths where every pixel
  # falls on only one side of the mask (can't run a two-sample test then).
  sensor_summary <- sensor_deviation |>
    summarise(n = n(), mean_abs_deviation = mean(abs_deviation, na.rm = TRUE),
              mean_rank = round(mean(rank, na.rm = TRUE), 2), .by = c(wavelength, in_oyster_bed)) |>
    mutate(sensor_Y = sY, .before = 1) |>
    arrange(wavelength, in_oyster_bed)
  oyster_bed_summary_all <- bind_rows(oyster_bed_summary_all, sensor_summary)

  sensor_test <- sensor_deviation |>
    reframe({
      if(n_distinct(in_oyster_bed) == 2 && n() >= 10){
        wt <- suppressWarnings(wilcox.test(abs_deviation ~ in_oyster_bed))
        tibble(p_value = wt$p.value)
      } else {
        tibble(p_value = NA_real_)
      }
    }, .by = wavelength) |>
    mutate(sensor_Y = sY, .before = 1)
  oyster_bed_test_all <- bind_rows(oyster_bed_test_all, sensor_test)

  # 3x3-grid equivalent of the whole block above -- same structure, but
  # restricted to the 9 inner (grid_class == "3x3") pixels and using the
  # deviation_3x3/rank_3x3 columns (deviation from grid3x3_mean) instead of
  # the 5x5 ones.
  sensor_deviation_3x3 <- sensor_deviation |> dplyr::filter(grid_class == "3x3")

  p_dev_3x3 <- sensor_deviation_3x3 |>
    ggplot(aes(x = in_oyster_bed, y = abs_deviation_3x3, fill = in_oyster_bed)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~wavelength_facet, scales = "free_y") +
    labs(title = paste0(sY, " -- pixel deviation from 3x3 grid mean, in vs. out of oyster bed"),
         x = "In oyster bed", y = "|Sat pixel RHOW - grid mean|", fill = "In oyster bed") +
    theme_bw()
  ggsave(file.path(out_dir, paste0("deviation_vs_oyster_bed_3x3_", sY, ".png")), p_dev_3x3, width = 10, height = 8)

  p_rank_3x3 <- sensor_deviation_3x3 |>
    ggplot(aes(x = in_oyster_bed, y = rank_3x3, fill = in_oyster_bed)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~wavelength_facet) +
    labs(title = paste0(sY, " -- deviation rank (1 = biggest outlier), in vs. out of oyster bed"),
         x = "In oyster bed", y = "Rank within 3x3 grid (1 = furthest from mean)", fill = "In oyster bed") +
    theme_bw()
  ggsave(file.path(out_dir, paste0("deviation_rank_vs_oyster_bed_3x3_", sY, ".png")), p_rank_3x3, width = 10, height = 8)

  sensor_summary_3x3 <- sensor_deviation_3x3 |>
    summarise(n = n(), mean_abs_deviation = mean(abs_deviation_3x3, na.rm = TRUE),
              mean_rank = round(mean(rank_3x3, na.rm = TRUE), 2), .by = c(wavelength, in_oyster_bed)) |>
    mutate(sensor_Y = sY, .before = 1) |>
    arrange(wavelength, in_oyster_bed)
  oyster_bed_summary_all_3x3 <- bind_rows(oyster_bed_summary_all_3x3, sensor_summary_3x3)

  sensor_test_3x3 <- sensor_deviation_3x3 |>
    reframe({
      if(n_distinct(in_oyster_bed) == 2 && n() >= 10){
        wt <- suppressWarnings(wilcox.test(abs_deviation_3x3 ~ in_oyster_bed))
        tibble(p_value = wt$p.value)
      } else {
        tibble(p_value = NA_real_)
      }
    }, .by = wavelength) |>
    mutate(sensor_Y = sY, .before = 1)
  oyster_bed_test_all_3x3 <- bind_rows(oyster_bed_test_all_3x3, sensor_test_3x3)

  # In-situ-minus-satellite version of the 3x3 test: same structure as
  # sensor_deviation_3x3 above, but using abs_deviation_hyp/rank_hyp_3x3 (how
  # far each pixel sits from the paired HYPERNETS value) instead of
  # abs_deviation_3x3/rank_3x3 (how far each pixel sits from its own grid's
  # mean). This is the more direct test of the manuscript question -- actual
  # matchup quality, not just internal grid consistency.
  p_dev_hyp_3x3 <- sensor_deviation_3x3 |>
    ggplot(aes(x = in_oyster_bed, y = abs_deviation_hyp, fill = in_oyster_bed)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~wavelength_facet, scales = "free_y") +
    labs(title = paste0(sY, " -- pixel deviation from HYPERNETS (3x3 grid), in vs. out of oyster bed"),
         x = "In oyster bed", y = "|Hyp RHOW - sat pixel RHOW|", fill = "In oyster bed") +
    theme_bw()
  ggsave(file.path(out_dir, paste0("deviation_vs_oyster_bed_hyp_3x3_", sY, ".png")), p_dev_hyp_3x3, width = 10, height = 8)

  p_rank_hyp_3x3 <- sensor_deviation_3x3 |>
    ggplot(aes(x = in_oyster_bed, y = rank_hyp_3x3, fill = in_oyster_bed)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~wavelength_facet) +
    labs(title = paste0(sY, " -- deviation-from-HYPERNETS rank (1 = worst, 3x3 grid), in vs. out of oyster bed"),
         x = "In oyster bed", y = "Rank within 3x3 grid (1 = furthest from Hyp)", fill = "In oyster bed") +
    theme_bw()
  ggsave(file.path(out_dir, paste0("deviation_rank_vs_oyster_bed_hyp_3x3_", sY, ".png")), p_rank_hyp_3x3, width = 10, height = 8)

  sensor_summary_hyp_3x3 <- sensor_deviation_3x3 |>
    summarise(n = n(), mean_abs_deviation_hyp = mean(abs_deviation_hyp, na.rm = TRUE),
              mean_rank_hyp = round(mean(rank_hyp_3x3, na.rm = TRUE), 2), .by = c(wavelength, in_oyster_bed)) |>
    mutate(sensor_Y = sY, .before = 1) |>
    arrange(wavelength, in_oyster_bed)
  oyster_bed_summary_all_hyp_3x3 <- bind_rows(oyster_bed_summary_all_hyp_3x3, sensor_summary_hyp_3x3)

  sensor_test_hyp_3x3 <- sensor_deviation_3x3 |>
    reframe({
      if(n_distinct(in_oyster_bed) == 2 && n() >= 10){
        wt <- suppressWarnings(wilcox.test(abs_deviation_hyp ~ in_oyster_bed))
        tibble(p_value = wt$p.value)
      } else {
        tibble(p_value = NA_real_)
      }
    }, .by = wavelength) |>
    mutate(sensor_Y = sY, .before = 1)
  oyster_bed_test_all_hyp_3x3 <- bind_rows(oyster_bed_test_all_hyp_3x3, sensor_test_hyp_3x3)
}

write_csv(oyster_bed_summary_all, file.path(out_dir, "oyster_bed_summary_all_sensors.csv"))
write_csv(oyster_bed_test_all, file.path(out_dir, "oyster_bed_test_all_sensors.csv"))
write_csv(oyster_bed_summary_all_3x3, file.path(out_dir, "oyster_bed_summary_all_sensors_3x3.csv"))
write_csv(oyster_bed_test_all_3x3, file.path(out_dir, "oyster_bed_test_all_sensors_3x3.csv"))
write_csv(oyster_bed_summary_all_hyp_3x3, file.path(out_dir, "oyster_bed_summary_all_sensors_hyp_3x3.csv"))
write_csv(oyster_bed_test_all_hyp_3x3, file.path(out_dir, "oyster_bed_test_all_sensors_hyp_3x3.csv"))
print(oyster_bed_summary_all, n = 100)
print(oyster_bed_test_all, n = 100)
print(oyster_bed_summary_all_3x3, n = 100)
print(oyster_bed_test_all_3x3, n = 100)
print(oyster_bed_summary_all_hyp_3x3, n = 100)
print(oyster_bed_test_all_hyp_3x3, n = 100)


# Single-matchup inspection --------------------------------------------------
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
# pace_nm_breaks / pace_nm_labels defined once in Setup, above.

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

