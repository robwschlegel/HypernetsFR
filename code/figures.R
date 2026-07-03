# code/figures.R
# It does what it says on the tin


# Setup -------------------------------------------------------------------

source("code/functions.R")

# For tailor diagrams
# library(openair)

# Luna package is used to access large spatial data products
# install.packages('luna', repos = 'https://rspatial.r-universe.dev')
# library(terra)
# library(luna)


# TODO: Rather change the colour palette so that the bands for each satellite are shown in 
# the legend for what exactly they are. Keeping the current spectrum range colours open
# for PACE.

# Global scatterplots -----------------------------------------------------

# TODO: Add the stacked scatterplots by sensors
# Get the code from HypernetsTara
global_scatterplot("S3A")
global_scatterplot("S3B")
global_scatterplot("S3")

# Satellite
# global_scatterplot_stack("MODIS") # To add
# global_scatterplot_stack("VIIRS")
# global_scatterplot_stack("OLCI")
# global_scatterplot_stack("OCI")


# Figure 1 ---------------------------------------------------------------

## Station data ------------------------------------------------------------

# Get GPS coordinates
file_list_GPS <- list.files("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/Tara_GPS", pattern = "*.csv", full.names = TRUE)
base_GPS <- plyr::ldply(file_list_GPS, read_csv, .parallel = TRUE, skip = 2, 
                        col_names = c("time",	"datetime", "latitude",	"longitude", "altitude", "gps_qual", "num_sats",
                                      "horizontal_dil",	"true_course", "true_track", "mag_track",	"spd_over_grnd_kmph"))
station_GPS <- base_GPS |> 
  dplyr::select(time:longitude) |> 
  filter(!is.nan(latitude), !is.nan(longitude))
station_GPS_unique <- station_GPS |> 
  dplyr::select(longitude, latitude) |> 
  mutate(longitude = round(longitude, 2),
         latitude = round(latitude, 2)) |> 
  distinct()

# Get raw Trios coordinates
station_Trios <- map_dfr(dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/Trios_processed_data/", full.names = TRUE),
                         read_csv, col_select = 1:4, skip = 31, col_names = c("date", "time", "lat", "lon"))
station_Trios_unique <- station_Trios |> 
  dplyr::rename(latitude = lat, longitude = lon) |> 
  mutate(date = as.Date(as.character(date), format = "%Y%m%d")) |> 
  filter(date >= "2024-08-06", date <= "2024-08-21") |>
  dplyr::select(date, longitude, latitude) |> 
  mutate(longitude = round(longitude, 2),
         latitude = round(latitude, 2)) |> 
  distinct()

# Get raw HyperPRO coords
# station_Pro <- map_dfr(dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/hyperpro_processed_data", full.names = TRUE), load_HyperPRO_coords)
# station_Pro_unique <- station_Pro |> 
#   dplyr::rename(latitude = lat, longitude = lon) |> 
#   dplyr::select(longitude, latitude) |> 
#   mutate(longitude = round(longitude, 2),
#          latitude = round(latitude, 2)) |> 
#   distinct()

# Get HYPERNETS raw coordinates
# NB: Only contains four distinct points
station_Hyp <- plyr::ldply(dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/Hypernets_processed_data", 
                               recursive = TRUE, pattern = "\\.nc", full.names = TRUE), load_HYPERNETS_coords, .parallel = TRUE)
station_Hyp_unique <- station_Hyp |> 
  dplyr::rename(latitude = site_latitude, longitude = site_longitude) |> 
  dplyr::select(longitude, latitude) |> 
  # mutate(longitude = round(longitude, 2),
         # latitude = round(latitude, 2)) |> 
  distinct()

# Get file lists for all in situ measurements
file_list_Hyp <- list.files(dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260203", 
                                pattern = "HYPERNETS", full.names = TRUE), pattern = "*.csv", full.names = TRUE)
file_list_Hyp <- file_list_Hyp[!grepl("HYPERNETS_vs_TRIOS_vs_HYPERPRO", file_list_Hyp)]
file_list_Trios <- list.files(dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260203", 
                                pattern = "TRIOS", full.names = TRUE), pattern = "*.csv", full.names = TRUE)
file_list_Trios <- file_list_Trios[!grepl("HYPERNETS_vs_TRIOS_vs_HYPERPRO", file_list_Trios)]
file_list_Pro <- list.files(dir("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260203", 
                                pattern = "HYPERPRO", full.names = TRUE), pattern = "*.csv", full.names = TRUE)
file_list_Pro <- file_list_Pro[!grepl("HYPERNETS_vs_TRIOS_vs_HYPERPRO", file_list_Pro)]

# Load all matchups made with in situ data
# NB: This is heavy, rather just load the output below
# base_Hyp <- plyr::ldply(file_list_Hyp, load_matchup_mean, .parallel = TRUE)
# base_Trios <- plyr::ldply(file_list_Trios, load_matchup_mean, .parallel = TRUE)
# base_Pro <- plyr::ldply(file_list_Pro, load_matchup_mean, .parallel = TRUE)

# Filter out satellite coordinates and reduce to unique measurements
# station_in_situ <- bind_rows(base_Hyp, base_Trios, base_Pro) |> 
#   dplyr::select(sensor:longitude) |> 
#   filter(sensor %in% c("Hyp", "HYPERPRO", "TRIOS")) |> 
#   distinct() |> 
#   mutate(date = parse_date(day, format = "%Y%m%d", locale = locale(tz = "UTC")),
#          dateTime = as.POSIXct(paste0(day," ",time), format = "%Y%m%d %H%M%S", locale = locale(tz = "UTC")), .before = "latitude")
# write_csv(station_in_situ, "meta/station_in_situ.csv")
station_in_situ <- read_csv("meta/station_in_situ.csv")

# Get count of samples per ~1 km lon/lat
station_count <- station_in_situ |> 
  mutate(longitude = round(longitude, 2),
         latitude = round(latitude, 2)) |> 
  summarise(count_n = n(), .by = c("longitude", "latitude"))


## MODIS data --------------------------------------------------------------

# NB: Most of the code here only needed to be run once

# Load username and password
# earth_up <- read_csv("~/pCloudDrive/Documents/info/earthdata_pswd.csv")

# Lists all products that are currently searchable
# https://lpdaac.usgs.gov/documents/925/MOD09_User_Guide_V61.pdf
# MODIS_prod <- luna::getProducts()
# productInfo("MODISA_L3m_RRS")

# MODIS/Aqua Surface Reflectance Daily L2G Global 250m SIN Grid V061
# productInfo("MYD09GQ")

# Niveau 1
# productInfo("MYD01")

# MODIS/Aqua Surface Reflectance 8-Day L3 Global 250m SIN Grid V006
# productInfo("MYD09Q1")

# MODIS/Aqua Surface Reflectance 8-Day L3 Global 500m SIN Grid V006
# moproductInfo("MYD09A1")

# MODIS/Terra Land Water Mask Derived from MODIS and SRTM L3 Global 250m SIN Grid V061
# https://lpdaac.usgs.gov/documents/1915/MOD44W_User_Guide_ATBD_V61.pdf
# productInfo("MOD44W")

# Coords for bbox
# coords <- matrix(c(
#   3, 32,  # Bottom-left corner
#   27, 32, # Bottom-right corner
#   27, 45, # Top-right corner
#   3, 45,  # Top-left corner
#   3, 32   # Close the polygon (same as first point)
# ), ncol = 2, byrow = TRUE)

# Create a SpatVector object
# bbox <- vect(coords, crs = "EPSG:4326", type = "polygons")

# Print the object to verify
# print(bbox)
# plot(bbox)

# Chosen start and end dates for downloading
# dl_start <- "2024-08-09"; dl_end <- "2024-08-16"

# Download data
# MODIS_dl("MYD09A1", dl_start, bbox, earth_up$usrname, earth_up$psswrd)

# Set file pathways
# mask_files <- list.files(path = "data/MODIS", pattern = "MOD", full.names = TRUE)
# rast_files <- list.files(path = "data/MODIS", pattern = "MYD", full.names = TRUE)

# Run on all of them
# MODIS_mask <- MODIS_proc(mask_files, bbox = bbox, water_mask = TRUE)
# writeRaster(MODIS_mask, "data/MODIS/study_area_mask.tif", overwrite = TRUE)
MODIS_mask <- rast("data/MODIS/study_area_mask.tif")
plot(MODIS_mask)

# Prep one day of MODIS data
# MODIS_rast <- MODIS_proc(rast_files, bbox = bbox)
# writeRaster(MODIS_rast, "data/MODIS/study_area_rast.tif", overwrite = TRUE)
MODIS_rast <- rast("data/MODIS/study_area_rast.tif")
plot(MODIS_rast)

# Projected the 250 m mask to the same grid as the 500 m raster data
MODIS_mask_proj <- project(MODIS_mask, MODIS_rast)

# Mask the raster data
MODIS_water <- mask(MODIS_rast, MODIS_mask_proj)
plot(MODIS_water)

# Convert to data.frame for easy plotting
MODIS_water_df <- as.data.frame(MODIS_water, xy = TRUE, na.rm = TRUE) |> 
  mutate(sur_refl_b03 = case_when(sur_refl_b03 > 0.1 ~ 0.1, 
                                  sur_refl_b03 < 0.0 ~ 0.0,
                                  TRUE ~ sur_refl_b03))


## Map of FR Stations -----------------------------------------------------

# Map
pl_map <- ggplot(data = station_in_situ, aes(x = longitude, y = latitude)) +
  annotation_borders(fill = "grey80") +
  geom_tile(data = MODIS_water_df, aes(x = x, y = y, fill = sur_refl_b03)) +
  # geom_path(data = station_GPS_unique) +
  geom_point(data = station_GPS_unique, size = 0.1) +
  # Show Trios raw collection points
  # geom_path(data = station_Trios_unique, colour = "red") +
  # geom_point(data = station_Trios_unique, colour = "black") +
  # The sampling stations
  geom_point(colour = "black", size = 5.5) +
  geom_point(aes(colour = as.factor(date)), size = 5) +
  geom_point(data = filter(station_in_situ, sensor == "HYPERPRO"), 
             colour = "maroon", size = 7, shape = 5, stroke = 2) +
  scale_fill_viridis_c(option = "D", breaks = c(0.03, 0.06, 0.1), labels = c("0.03", "0.06", ">0.10")) +
  scale_colour_discrete(palette = "Set3") +
  # scale_fill_viridis_c() +
  # guides(fill = guide_colorbar(barwidth = 20, barheight = 2)) +
  guides(fill = guide_colourbar(theme = theme(legend.key.width  = unit(25, "lines"),
                                              legend.key.height = unit(1.5, "lines")))) +
  guides(color = guide_legend(nrow = 2,
                              override.aes = list(shape = 21, color = "black", size = 5,
                                                  fill = c(RColorBrewer::brewer.pal(n = 11, name = "Set3"))))) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)", 
       colour = "Measurement<br>date", 
       fill = "<i>ρ<sub>w</sub></i>  (Band 3; 459-479 nm)") +
  coord_quickmap(xlim = c(7, 25), ylim = c(35, 42)) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top",
        legend.box = "vertical",
        # legend.key.width = unit(2, "cm"),
        legend.title = element_markdown(size = 26, halign = 0.5),
        legend.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 22))
# pl_map
ggsave("figures/fig_1.png", pl_map, height = 12, width = 18)


# Figure 3 ----------------------------------------------------------------

# Create a hyperspectral plot that shows all in situ sensors vs PACE and one other multispectral satellite
# Add variance bars for multispectral data
# Plus the photos from HYPERNETS

# Using the HYPERPRO 2024-08-12 10:53:00 measurement as a reference as this has a PACE and VIIRS_N matchup
pro_pace2 <- read_delim("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260504/RHOW_HYPERPRO_vs_PACE_V2/HYPERPRO_vs_PACE_V2_vs_20240812T105300_RHOW.csv", delim = ";", col_types = "ccccnnic")
pro_pace3 <- read_delim("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260504/RHOW_HYPERPRO_vs_PACE_V30/HYPERPRO_vs_PACE_V30_vs_20240812T105300_RHOW.csv", delim = ";", col_types = "ccccnnic")
pro_pace31 <- read_delim("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260504/RHOW_HYPERPRO_vs_PACE_V31/HYPERPRO_vs_PACE_V31_vs_20240812T105300_RHOW.csv", delim = ";", col_types = "ccccnnic")
pro_viirsn <- read_delim("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260504/RHOW_HYPERPRO_vs_VIIRS_N/HYPERPRO_vs_VIIRS_N_vs_20240812T105300_RHOW.csv", delim = ";", col_types = "ccccnnic")
pro_hyp <- read_delim("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260504/RHOW_HYPERNETS_vs_HYPERPRO/HYPERNETS_vs_HYPERPRO_vs_20240812T104500_RHOW.csv", delim = ";", col_types = "ccccnnic")
pro_tri <- read_delim("~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260504/RHOW_TRIOS_vs_HYPERPRO/TRIOS_vs_HYPERPRO_vs_20240812T104536_RHOW.csv", delim = ";", col_types = "ccccnnic")

# Load the RGB images taken by HYPERNETS at the closest available timing
hyp_LD <- "~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/Hypernets_processed_data/12/SEQ20240812T104553/image/HYPERNETS_W_MAFR_IMG_20240812T1045_20240912T1030_006_140_90_v2.0.jpg"
hyp_LU <- "~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/Hypernets_processed_data/12/SEQ20240812T104553/image/HYPERNETS_W_MAFR_IMG_20240812T1045_20240912T1030_009_40_90_v2.0.jpg"
hyp_ED <- "~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/Hypernets_processed_data/12/SEQ20240812T104553/image/HYPERNETS_W_MAFR_IMG_20240812T1045_20240912T1030_015_180_90_v2.0.jpg"
hyp_SUN <- "~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/Hypernets_processed_data/12/SEQ20240812T104553/image/HYPERNETS_W_MAFR_IMG_20240812T1045_20240912T1030_016_0_0_v2.0.jpg"

# Combine wide for reference
pro_all <- rbind(pro_pace2, pro_pace3, pro_pace31, pro_viirsn, pro_hyp, pro_tri) |> 
  `colnames<-`(c("sensor", colnames(pro_pace2)[2:329])) |> 
  filter(!sensor %in% c("HYPERPRO 4", "HYPERPRO 5", "HYPERPRO 6", "HYPERPRO 7")) |> 
  distinct() |> mutate(sensor = gsub(" 1$| 2$| 3$| 4$| 5$| 6$| 7$| 8$| 9$", "", sensor))

# Melt and prep for plotting
pro_all_long <- pro_all |> 
  dplyr::select(-(day:type)) |> 
  pivot_longer(`380`:`700`, names_to = "wavelength") |> 
  filter(!is.na(value)) |> 
  mutate(wavelength = as.numeric(wavelength),
         data_type = case_when(data_type == "weighted" ~ "rhow", TRUE ~ data_type),
         sensor = case_when(sensor == "HYPERPRO" ~ "HyperPRO",
                            sensor == "PACE_V2" ~ "PACE v2.0",
                            sensor == "PACE_V30" ~ "PACE v3.0",
                            sensor == "PACE_V31" ~ "PACE v3.1",
                            sensor == "VIIRS_N" ~ "VIIRS SNPP",
                            sensor == "Hyp" ~ "HYPERNETS",
                            sensor == "Hyp_nosc" ~ "HYPERNETS (nosc)",
                            sensor == "TRIOS" ~ "So-Rad")) |> 
  pivot_wider(names_from = data_type, values_from = value) |> 
  filter(wavelength >= 380, wavelength <= 700) |> 
  filter(sensor != "HYPERNETS (nosc)") |>  # Not interesting for this match-up
  mutate(sensor = factor(sensor, levels = c("HyperPRO", "So-Rad", "HYPERNETS", "PACE v2.0", "PACE v3.0", "PACE v3.1", "VIIRS SNPP")))

# The hyperspectral plot
pl_spect <- ggplot(data = pro_all_long, aes(x = wavelength, y = rhow, colour = sensor, fill = sensor)) +
  # geom_ribbon(data = filter(pro_all_long, sensor != "VIIRS SNPP"),
  # aes(ymin = std_min, ymax = std_max), alpha = 0.1, colour = NA) +
  geom_line(data = filter(pro_all_long, sensor != "VIIRS SNPP"), linewidth = 2, alpha = 0.7) +
  geom_point(data = filter(pro_all_long, sensor == "VIIRS SNPP"), size = 3, show.legend = FALSE) +
  geom_errorbar(data = filter(pro_all_long, sensor == "VIIRS SNPP"), 
                aes(ymin = std_min, ymax = std_max), width = 5) +
  # scale_colour_brewer(palette = "Dark2") +
  scale_colour_manual(values = c("darkred", "darkorange", "goldenrod", "skyblue", "dodgerblue", "royalblue", "limegreen")) +
  # scale_x_continuous(expand = FALSE) +
  guides(colour = guide_legend(override.aes = list(alpha = 1.0, linewidth = 3))) +
  labs(x = "Wavelength (nm)",
       y = "<i>ρ<sub>w</sub></i>",
       colour = "Sensor") +
  theme_minimal() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.75),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_markdown(size = 12),
        axis.text = element_text(size = 10))

# Prep each jpg
pl_ld <- ggplot() +
  geom_image(aes(x = 0, y = 0), image = hyp_LD, size = 1) +
  labs(title = "<i>L<sub>d</sub></i>") + theme_void() +
  theme(plot.title = element_markdown(hjust = 0.5, vjust = 0))
pl_lu <- ggplot() +
  geom_image(aes(x = 0, y = 0), image = hyp_LU, size = 1) +
  labs(title = "<i>L<sub>u</sub></i>") + theme_void() +
  theme(plot.title = element_markdown(hjust = 0.5, vjust = 0))
pl_ed <- ggplot() +
  geom_image(aes(x = 0, y = 0), image = hyp_ED, size = 1) +
  labs(title = "<i>E<sub>d</sub></i>") + theme_void() +
  theme(plot.title = element_markdown(hjust = 0.5, vjust = 0))
pl_sun <- ggplot() +
  geom_image(aes(x = 0, y = 0), image = hyp_SUN, size = 1) +
  labs(title = "<i>Sun</i>") + theme_void() +
  theme(plot.title = element_markdown(hjust = 0.5, vjust = 0))

# Combine photos into one row
fig_3_bottom <- ggpubr::ggarrange(pl_ld, pl_lu, pl_ed, pl_sun, nrow = 1, ncol = 4)
fig_3 <- ggpubr::ggarrange(pl_spect, fig_3_bottom, ncol = 1, nrow = 2, heights = c(2, 1)) +
  ggpubr::bgcolor("white") + ggpubr::border("white", size = 2)
ggsave("figures/fig_3.png", fig_3, width = 12, height = 9)


# Figure 4 ----------------------------------------------------------------

# Matchup scatterplots for everything shown in Figure 3, with statistic panels in each pane;
fig_4_a <- plot_matchup_single_nm(pro_all_long, "HyperPRO", "So-Rad")
fig_4_b <- plot_matchup_single_nm(pro_all_long, "HyperPRO", "HYPERNETS")
fig_4_c <- plot_matchup_single_nm(pro_all_long, "HyperPRO", "VIIRS SNPP")
fig_4_d <- plot_matchup_single_nm(pro_all_long, "HyperPRO", "PACE v2.0")
fig_4_e <- plot_matchup_single_nm(pro_all_long, "HyperPRO", "PACE v3.0")
fig_4_f <- plot_matchup_single_nm(pro_all_long, "HyperPRO", "PACE v3.1")

# Combine and save
fig_4 <- ggpubr::ggarrange(fig_4_a, fig_4_b, fig_4_c, fig_4_d, fig_4_e, fig_4_f, 
                           labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
                           ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom") +
  ggpubr::bgcolor("white") + ggpubr::border("white", size = 2)
ggsave("figures/fig_4.png", fig_4, width = 12, height = 9)


# Figure 11 ---------------------------------------------------------------

# Tailor diagram showing comparison of all satellites per waveband

# Load all matchup data at once into one dataframe and filter for wavebands of interests

## Load single matchup results
df_matchups_single <- map_dfr(dir(path = "output", pattern = "matchup_stats_RHOW", full.names = TRUE), 
                              read_csv, show_col_types = FALSE) |> 
  filter(sensor_X %in% c("Hyp", "TRIOS", "HYPERPRO"),
        !(sensor_Y %in% c("Hyp", "TRIOS", "HYPERPRO")))

## Load global matchuups
df_matchups_global <- read_csv("output/global_stats_all.csv", show_col_types = FALSE) |> 
  filter(sensor_X %in% c("HYPERNETS", "TRIOS", "HYPERPRO"),
        !(sensor_Y %in% c("HYPERNETS", "TRIOS", "HYPERPRO"))) |> 
  mutate(sensor_Y = case_when(sensor_Y == "AQUA" ~ "Aqua",
                              sensor_Y == "S3A" ~ "S3A",
                              sensor_Y == "S3B" ~ "S3B",
                              sensor_Y == "S3_all" ~ "S3 all",
                              sensor_Y == "PACE_V2" ~ "PACE v2.0",
                              sensor_Y == "PACE_V30" ~ "PACE v3.0",
                              sensor_Y == "PACE_V31" ~ "PACE v3.1",
                              sensor_Y == "VIIRS_N" ~ "SNPP",
                              sensor_Y == "VIIRS_J1" ~ "JPSS1",
                              sensor_Y == "VIIRS_J2" ~ "JPSS2"),
        sensor_X = case_when(sensor_X == "HYPERPRO" ~ "HyperPRO",
                             sensor_X == "TRIOS" ~ "So-Rad",
                             sensor_X == "HYPERNETS" ~ "HYPERNETS"))

## List all files used in matchups
### NB: Careful with exact indexing
match_base_details <- df_matchups_single |> dplyr::select(file_name) |> distinct()

## List all RHOW matchup satellite files and filter to those that passed QC
file_list_sat <- dir(path = "~/pCloudDrive/Documents/OMTAB/HYPERNETS/Tara/tara_matchups_results_20260504", 
                     pattern = "RHOW", full.names = TRUE, recursive = TRUE)
file_list_sat <- file_list_sat[grepl("AQUA|PACE|S3A|S3B|VIIRS", file_list_sat)]
file_list_sat <- file_list_sat[basename(file_list_sat) %in% match_base_details$file_name]

## Split into the three sensors for ease of management below
file_list_sat_hyperpro <- file_list_sat[grepl("RHOW_HYPERPRO", file_list_sat)]
file_list_sat_trios <- file_list_sat[grepl("RHOW_TRIOS", file_list_sat)]
file_list_sat_hypernets <- file_list_sat[grepl("RHOW_HYPERNETS", file_list_sat)]

## Load all needed files
df_sat_hyperpro <- map_dfr(file_list_sat_hyperpro, load_matchup_long) |> 
  pivot_longer(AQUA:VIIRS_N, values_to = "value_sat", names_to = "sensor_sat") |> 
  mutate(sensor_is = "HYPERPRO") |> dplyr::rename(value_is = HYPERPRO)
df_sat_trios <- map_dfr(file_list_sat_trios, load_matchup_long) |> 
  pivot_longer(AQUA:VIIRS_N, values_to = "value_sat", names_to = "sensor_sat") |> 
  mutate(sensor_is = "TRIOS") |> dplyr::rename(value_is = TRIOS)
df_sat_hypernets <- map_dfr(file_list_sat_hypernets, load_matchup_long) |> 
  pivot_longer(AQUA:VIIRS_N, values_to = "value_sat", names_to = "sensor_sat") |> 
  mutate(sensor_is = "Hyp") |> dplyr::rename(value_is = Hyp)

## Combine into one dataframe
df_sat_all <- bind_rows(df_sat_hyperpro, df_sat_trios, df_sat_hypernets) |>
  mutate(value_sat = case_when(sensor_sat %in% c("PACE_V2", "PACE_V30", "PACE_V31") & 
                                  !(wavelength %in% c(380, 400, 412, 443, 490, 510, 560, 620, 673, 700)) ~ NA , TRUE ~ value_sat)) |>
  filter(!is.na(value_sat)) |>
  mutate(sensor_sat = case_when(sensor_sat == "AQUA" ~ "Aqua MODIS",
                                sensor_sat == "S3A" ~ "Sentinel-3A OLCI",
                                sensor_sat == "S3B" ~ "Sentinel-3B OLCI",
                                sensor_sat == "PACE_V2" ~ "PACE v2.0",
                                sensor_sat == "PACE_V30" ~ "PACE v3.0",
                                sensor_sat == "PACE_V31" ~ "PACE v3.1",
                                sensor_sat == "VIIRS_N" ~ "VIIRS SNPP",
                                sensor_sat == "VIIRS_J1" ~ "VIIRS J1",
                                sensor_sat == "VIIRS_J2" ~ "VIIRS J2"),
         sensor_is = case_when(sensor_is == "HYPERPRO" ~ "HyperPRO",
                               sensor_is == "TRIOS" ~ "So-Rad",
                               sensor_is == "Hyp" ~ "HYPERNETS"),
         sensor_group = paste0(sensor_is , " vs ", sensor_sat, " at ", wavelength))

# NB: This doesn't work well for everything in one whack
# Rather need to subset the data further
# TaylorDiagram(
#   mydata       = df_sat_all[df_sat_all$sensor_is == "HyperPRO",],
#   obs          = "value_is",
#   mod          = "value_sat",
#   group        = "sensor_group",
#   main         = "Taylor Diagram — All Wavebands",
#   normalise    = TRUE,        # normalise so bands are comparable
#   cols         = colour_nm[1:7],
#   pch          = 19,
#   cex          = 1.1,
#   key.title    = "Wavelength",
#   annotate     = "RMSE"
# )

# Another option is to create correlation matrices
# Or potentially heatmaps of statistics for matchups

# Prep data for matrix plots
df_matchups_global_pretty <- df_matchups_global |> 
  filter(Wavelength_nm >= 400, Wavelength_nm <= 600) |> 
  mutate(sensor_sat = case_when(sensor_Y %in% c("Aqua") ~ "MODIS",
                                sensor_Y %in% c("S3A", "S3B", "S3 all") ~ "OLCI",
                                sensor_Y %in% c("PACE v2.0", "PACE v3.0", "PACE v3.1") ~ "OCI",
                                sensor_Y %in% c("SNPP", "JPSS1", "JPSS2") ~ "VIIRS"),
         wavelength_clean = case_when(sensor_sat == "VIIRS" & Wavelength_nm %in% c(410, 411) ~ "410/411", 
                                      sensor_sat == "VIIRS" & Wavelength_nm %in% c(443, 445) ~ "443/445",
                                      sensor_sat == "VIIRS" & Wavelength_nm %in% c(486, 489) ~ "486/489",
                                      sensor_sat == "VIIRS" & Wavelength_nm %in% c(551, 556) ~ "551/556",  
                                      TRUE ~ as.character(Wavelength_nm))) |> 
  mutate(sensor_Y = factor(sensor_Y, levels = c("Aqua", "S3 all", "S3B", "S3A", 
                                                "PACE v3.1", "PACE v3.0", "PACE v2.0", 
                                                "JPSS2", "JPSS1", "SNPP")),
         sensor_sat = factor(sensor_sat, levels = c("MODIS", "OLCI", "OCI", "VIIRS")),
         sensor_X = factor(sensor_X, levels = c("HyperPRO", "So-Rad", "HYPERNETS")))

# Matrix plot
plot_matrix_error <- function(df, val_range){
  # Create column for all wavelengths
  df_all <- df |> 
    group_by(sensor_Y, sensor_X, sensor_sat) |> 
    summarise(Error = mean(Error, na.rm = TRUE), .groups = "drop") |> 
    mutate(wavelength_clean = "All")
  df <- bind_rows(df, df_all) #|> 
    # mutate(wavelength_clean = factor(wavelength_clean)#, 
      # levels = c("410/411", "443/445", "486/489", "551/556", "400", "412", "443", "490", "510", "560", "All")))
  # Round data to given range
  df_round <- df |> 
    mutate(Error = case_when(Error > val_range[2] ~ val_range[2],
                             TRUE ~ Error))
  ggplot(data = df_round, aes(x = wavelength_clean, y = sensor_Y)) +
    geom_tile(aes(fill = Error), colour = "black") +
    geom_label(data = df, aes(label = sprintf("%.1f", round(Error, 1))), size = 3) +
    labs(x = "Waveband (nm)", y = NULL, fill = "Error (%)") +
    facet_grid(sensor_sat~sensor_X, scales = "free") +
    # scale_y_reverse() +
    scale_fill_viridis_c(limits = val_range, breaks = c(10, 20, 30), labels = c("10", "20", "30")) +
    # theme_minimal() +
    coord_cartesian(expand = FALSE) +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          # legend.position = "bottom",
          axis.title.x = element_markdown(size = 12),
          axis.title.y = element_markdown(size = 12),
          axis.text = element_text(size = 10))
}

# Select the sensors to plot and create a list of plots by sensor
sensors <- unique(df_matchups_global_pretty$sensor_sat)
plots_by_sensor <- purrr::set_names(
  purrr::map(sensors, function(s) {
    # val_range <- range(df_matchups_global_pretty$Error, na.rm = TRUE)
    val_range <- c(0, 40)
    df_sub <- df_matchups_global_pretty |> dplyr::filter(sensor_sat == s)
    plot_matrix_error(df_sub, val_range)
  }),
  sensors
)
# plots_by_sensor

fig_11 <- plots_by_sensor$MODIS / plots_by_sensor$VIIRS / plots_by_sensor$OLCI / plots_by_sensor$OCI + 
  plot_annotation(tag_levels = "a", tag_suffix = ")") +
  patchwork::plot_layout(guides = "collect", axis_titles = "collect", heights = c(0.35, 1, 1, 1))
ggsave("figures/fig_11.png", fig_11, width = 12, height = 9)


# Fig S1 ------------------------------------------------------------------

# Load one slice of PACE v3.1 data at 413 nm
## Visualise all five days to pick the best coverage
# map_PACE(read_csv("data/PACE_OCI.20240809T105059.L2.OC_AOP.V3_1_rrs_413.csv"))
# map_PACE(read_csv("data/PACE_OCI.20240812T105611.L2.OC_AOP.V3_1_rrs_413.csv"))
# map_PACE(read_csv("data/PACE_OCI.20240813T113041.L2.OC_AOP.V3_1_rrs_413.csv")) # Winner
# map_PACE(read_csv("data/PACE_OCI.20240814T120512.L2.OC_AOP.V3_1_rrs_413.csv"))
# map_PACE(read_csv("data/PACE_OCI.20240815T110624.L2.OC_AOP.V3_1_rrs_413.csv"))
PACE_swath <- read_csv("data/PACE/PACE_OCI.20240813T113041.L2.OC_AOP.V3_1_rrs_413.csv") |> 
  filter(!is.na(Rrs))

# Load the full spectra for one point
v_all_spectra <- read_delim("data/csv_pour_générer_spectres_pace_V20_V30_V31/PACE_20240809T0900.csv", delim = ";")
colnames(v_all_spectra)[1] <- "version"
v_all_spectra_long <- v_all_spectra |> 
  pivot_longer(`356`:`718`, names_to = "nm") |> 
  mutate(nm = as.integer(nm))

# Load data extracted via Python script
v2_413 <- read_csv("data/PACE/PACE_OCI.20240809T105059.L2.OC_AOP.V2_0.NRT_rrs_413.csv") |> 
  filter(!is.na(Rrs)) |> 
  mutate(longitude = plyr::round_any(longitude, 0.02),
         latitude = plyr::round_any(latitude, 0.02)) |>
  summarise(Rrs = mean(Rrs, na.rm = TRUE), .by = c("longitude", "latitude")) |> 
  dplyr::rename(v2_Rrs = Rrs)
v3_413 <- read_csv("data/PACE/PACE_OCI.20240809T105059.L2.OC_AOP.V3_0_rrs_413.csv") |> 
  filter(!is.na(Rrs)) |> 
  mutate(longitude = plyr::round_any(longitude, 0.02),
         latitude = plyr::round_any(latitude, 0.02)) |>
  summarise(Rrs = mean(Rrs, na.rm = TRUE), .by = c("longitude", "latitude")) |> 
  dplyr::rename(v3_Rrs = Rrs)
v31_413 <- read_csv("data/PACE/PACE_OCI.20240809T105059.L2.OC_AOP.V3_1_rrs_413.csv") |> 
  filter(!is.na(Rrs)) |> 
  mutate(longitude = plyr::round_any(longitude, 0.02),
         latitude = plyr::round_any(latitude, 0.02)) |>
  summarise(Rrs = mean(Rrs, na.rm = TRUE), .by = c("longitude", "latitude")) |> 
  dplyr::rename(v31_Rrs = Rrs)

# Determine the cuts for the discrete groups
# comp_cut <- c(-200, -150, -100, -50, -25, -5, 5, 25, 50, 100, 150, 200)
comp_cut <- c(-200, -100, -50, -25, -5, 5, 25, 200)

# Combine and compare
vall_413 <- left_join(v2_413, v3_413, by = join_by(latitude, longitude)) |> 
  left_join(v31_413, by = join_by(latitude, longitude)) |> 
  mutate(v2_v3 = round(((v2_Rrs / v3_Rrs)*100)-100),
         v2_v31 = round(((v2_Rrs / v31_Rrs)*100)-100),
         v3_v31 = round(((v3_Rrs / v31_Rrs)*100)-100)) |> 
  mutate(v2_v3 = ifelse(v2_v3 > 200, 200, v2_v3),
         v2_v31 = ifelse(v2_v31 > 200, 200, v2_v31),
         v3_v31 = ifelse(v3_v31 > 200, 200, v3_v31),
         v2_v3 = ifelse(v2_v3 < -200, -200, v2_v3),
         v2_v31 = ifelse(v2_v31 < -200, -200, v2_v31),
         v3_v31 = ifelse(v3_v31 < -200, -200, v3_v31)) |>
  # mutate(v2_v3 = cut(v2_v3, seq(-200, 200, 50)),
  #        v2_v31 = cut(v2_v31, seq(-200, 200, 50)),
  #        v3_v31 = cut(v3_v31, seq(-200, 200, 50)))
  mutate(v2_v3_cut = cut(v2_v3, comp_cut, include.lowest = TRUE),
         v2_v31_cut = cut(v2_v31, comp_cut, include.lowest = TRUE),
         v3_v31_cut = cut(v3_v31, comp_cut, include.lowest = TRUE))

# Pivot longer for plotting
v_comp_long <- vall_413 |> 
  dplyr::select(longitude, latitude, v2_v3:v3_v31) |>
  pivot_longer(v2_v3:v3_v31, names_to = "ver") |>
  mutate(ver = factor(ver, levels = c("v2_v3", "v2_v31", "v3_v31"),
                      labels = c("v2 / v3", "v2 / v3.1", "v3 / v3.1"))) |>
  filter(!is.na(value))

# Another version as discrete cuts
v_comp_long_cut <- vall_413 |> 
  dplyr::select(longitude, latitude, v2_v3_cut:v3_v31_cut) |>
  pivot_longer(v2_v3_cut:v3_v31_cut, names_to = "ver") |>
  mutate(ver = factor(ver, levels = c("v2_v3_cut", "v2_v31_cut", "v3_v31_cut"),
                      labels = c("v2.0 / v3.0", "v2.0 / v3.1", "v3.0 / v3.1"))) |> 
  filter(!is.na(value))

# Get map data for plotting
map_df <- map_data("world")

# Pre-determine discrete colours for plotting
pixel_fill_d <- RColorBrewer::brewer.pal(8, "RdBu")
pixel_fill_d <- c(pixel_fill_d[1:4], "white", pixel_fill_d[6], pixel_fill_d[8])

# Plot map differences
pl_top <- ggplot(data = v_comp_long_cut) +
  geom_polygon(data = map_df, fill = "grey80", #colour = "black",
               aes(x = long, y = lat, group = group)) +
  geom_tile(aes(x = longitude, y = latitude, fill = value)) +
  annotate(geom = "point", x = v_all_spectra$longitude[1], y = v_all_spectra$latitude[1]) +
  # geom_contour(aes(x = longitude, y = latitude, z = value),
  #              breaks = c(1.0), color = "black") +
  scale_fill_manual(values = pixel_fill_d) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)", fill = "Difference (%)") +
  facet_wrap(~ver) +
  coord_quickmap(xlim = c(min(v_comp_long$longitude), max(v_comp_long$longitude)),
                 ylim = c(min(v_comp_long$latitude), max(v_comp_long$latitude))) +
  theme(panel.border = element_rect(colour = "black", fill = "grey90"),
        legend.position = "top")
# pl_top

# Plot percent difference as non-map
pl_left <- v_comp_long_cut |> 
  summarise(cut_n = n(), .by = c("ver", "value")) |> 
  # complete(ver, value) |> 
  ggplot() +
  geom_col(aes(x = value, y = cut_n, fill = ver), 
           position = "dodge", colour = "black") +
  scale_y_continuous(expand = c(0, 2000), 
                     breaks = c(0, 50000, 100000, 150000, 200000, 250000),
                     labels = c("0", "50K", "100K", "150K", "200K", "250K")) +
  scale_fill_brewer(palette = "Accent") +
  # scale_fill_viridis_d(option = "A") + # yuck
  labs(x = "Difference (%)", y = "Pixel count (n)", fill = "Comparison") +
  theme(panel.border = element_rect(colour = "black", fill = "grey90"),
        legend.position = "bottom")
# pl_left

# Plot spectra differences
pl_right <- ggplot(v_all_spectra_long) +
  geom_path(aes(x = nm, y = value, colour = version),
            linewidth = 2, alpha = 0.8) +
  geom_vline(xintercept = 413) +
  labs(x = "Wavelength (nm)", y = "Remote sensing reflectance (Rrs)") +
  scale_color_brewer("Version", palette = "YlGnBu") +
  # scale_colour_viridis_d(option = "B") + # yuck
  scale_y_continuous(expand = c(-0.1, 0.005)) +
  theme(panel.border = element_rect(colour = "black", fill = "grey90"),
        legend.position = "bottom")
# pl_right

# Put it all together
pl_all <- (pl_top / (pl_left + pl_right)) + 
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ')')
ggsave("figures/fig_S1.png", pl_all, height = 9, width = 14)


# Fig S2 ------------------------------------------------------------------

# Uniquely for this plot
library(patchwork)

# Barplot of all PACE wavelength matchups
global_stats_OCI <- read_csv("output/global_stats_RHOW_OCI.csv")

# Prep for plotting
global_stats_OCI_pretty <- global_stats_OCI |> 
  filter(sensor_X %in% c("HYPERNETS", "TRIOS", "HYPERPRO")) |> 
  mutate(sensor_X = factor(sensor_X,
                           levels = c("HYPERPRO", "TRIOS", "HYPERNETS"),
                           labels = c("HyperPRO", "So-Rad", "HYPERNETS")),
         sensor_Y = factor(sensor_Y,
                           levels = c("PACE_V2", "PACE_V30", "PACE_V31"),
                           labels = c("PACE v2.0", "PACE v3.0", "PACE v3.1"))) |> 
  dplyr::rename(`Wavelength (nm)` = Wavelength_nm)

# Create the barplot
## Error
pl_Error_OCI <- ggplot(data = global_stats_OCI_pretty, aes(x = `Wavelength (nm)`, y = Error)) +
  geom_col(aes(fill = sensor_Y), position = "dodge", show.legend = FALSE) +
  labs(y = "Error  (%)") +
  facet_grid(sensor_Y~sensor_X) +
  # theme_minimal() +
  theme(panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        axis.text = element_text(size = 10))
# pl_Error_OCI

## Bias
pl_Bias_OCI <- ggplot(data = global_stats_OCI_pretty, aes(x = `Wavelength (nm)`, y = Bias)) +
  geom_col(aes(fill = sensor_Y), position = "dodge", show.legend = FALSE) +
  labs(y = "Bias  (%)") +
  facet_grid(sensor_Y~sensor_X) +
  # theme_minimal() +
  theme(panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        axis.text = element_text(size = 10))
# pl_Bias_OCI

# Steek'em
fig_S2 <- pl_Error_OCI / pl_Bias_OCI + plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave("figures/fig_S2.png", fig_S2, width = 9, height = 12)

