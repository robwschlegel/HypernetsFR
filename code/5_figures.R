# code/5_figures.R
# It does what it says on the tin


# Setup -------------------------------------------------------------------

source("code/0_functions.R")

# For tailor diagrams
# library(openair)

# Luna package is used to access large spatial data products
# install.packages('luna', repos = 'https://rspatial.r-universe.dev')
# library(terra)
# library(luna)


# Global scatterplots -----------------------------------------------------

# Satellite
global_scatterplot_stack("MODIS")
global_scatterplot_stack("VIIRS")
global_scatterplot_stack("OLCI")
global_scatterplot_stack("OCI")

# Stop here if running from terminal
if (!interactive()) quit(save = "no", status = 0)


# Figure 1 ---------------------------------------------------------------

station_in_situ <- read_csv("meta/station_in_situ.csv")

map_france <- map_data("world", region = "France")

pl_map <- ggplot() +
  geom_polygon(data = map_france, aes(x = long, y = lat, group = group),
               fill = "grey80", colour = "grey40", linewidth = 0.3) +
  geom_point(data = station_in_situ, aes(x = lon, y = lat, colour = site),
             size = 4) +
  geom_label(data = station_in_situ, aes(x = lon, y = lat, label = site),
             nudge_y = 0.6, size = 4) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Longitude (°E)", y = "Latitude (°N)", colour = "Site") +
  coord_quickmap(xlim = c(-5, 10), ylim = c(41, 52)) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "right",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
# pl_map
ggsave("figures/fig_1.png", pl_map, height = 10, width = 8)


# Figure 3 ----------------------------------------------------------------
# HYPERNETS hyperspectral spectrum vs satellite band-equivalent Rhow for one
# representative matchup date at MAFR. Dates with OCI + OLCI + VIIRS coverage:
# 20240605, 20240619, 20240811, 20241011, 20250711

match_site <- "MAFR"
match_date <- "20240811"

# Load one matchup file for site / sat_name / date_str; return long data frame
# with columns: sensor, wavelength, rhow, std_min, std_max.
# HYPERNETS rows contain the full spectrum; satellite rows contain band-only values.
load_fig3_file <- function(site, sat_name, date_str) {
  folder <- file_path_build(site, sat_name)
  files  <- list.files(folder, pattern = date_str, full.names = TRUE)
  if (length(files) == 0) return(NULL)

  # load_matchup_mean handles all data_type quirks across sensor families
  # (PACE "rhow", JPSS1 "rhow weighted", replicate-scan fallback, etc.)
  df_mean <- tryCatch(load_matchup_mean(files[1]), error = function(e) NULL)
  if (is.null(df_mean) || nrow(df_mean) == 0) return(NULL)

  sat_label <- df_mean |> filter(sensor != "Hyp") |> pull(sensor) |> first()

  # load_matchup_mean has already removed data_type/type/pixel_pos/variability_centered;
  # pivot everything that isn't a metadata column
  meta_cols <- c("sensor", "day", "time", "latitude", "longitude")
  df_long <- df_mean |>
    pivot_longer(-any_of(meta_cols), names_to = "wavelength", values_to = "rhow") |>
    filter(!is.na(rhow)) |>
    mutate(wavelength = as.numeric(wavelength))

  hyp <- df_long |>
    filter(sensor == "Hyp") |>
    mutate(sensor = "HYPERNETS", std_min = NA_real_, std_max = NA_real_)

  # Reload raw file only to extract pixel-box std_min / std_max for satellite bands
  suppressMessages(df_raw <- read_delim(files[1], delim = ";", col_types = "ccccnnic"))
  colnames(df_raw)[1] <- "sensor"
  df_raw <- df_raw |> mutate(sensor = gsub(" [0-9]+$", "", sensor))

  non_wave_cols <- c("sensor", "day", "time", "latitude", "longitude",
                     "radiometer_id", "type", "data_type", "pixel_pos", "variability_centered")
  sat_lo <- df_raw |> filter(sensor == sat_label, data_type == "std_min") |> slice(1) |>
    dplyr::select(-any_of(non_wave_cols)) |>
    pivot_longer(everything(), names_to = "wavelength", values_to = "std_min") |>
    filter(!is.na(std_min)) |>
    mutate(wavelength = as.numeric(wavelength))
  sat_hi <- df_raw |> filter(sensor == sat_label, data_type == "std_max") |> slice(1) |>
    dplyr::select(-any_of(non_wave_cols)) |>
    pivot_longer(everything(), names_to = "wavelength", values_to = "std_max") |>
    filter(!is.na(std_max)) |>
    mutate(wavelength = as.numeric(wavelength))

  sat <- df_long |>
    filter(sensor == sat_label) |>
    left_join(sat_lo |> dplyr::select(wavelength, std_min), by = "wavelength") |>
    left_join(sat_hi |> dplyr::select(wavelength, std_max), by = "wavelength")

  bind_rows(
    dplyr::select(hyp, sensor, wavelength, rhow, std_min, std_max),
    dplyr::select(sat, sensor, wavelength, rhow, std_min, std_max)
  )
}

# Try each sensor_Y in preference order; PACE first so its denser HYPERNETS
# spectrum is used as the canonical in-situ line
sat_names_fig3 <- c("PACE", "S3A", "S3B", "SNPP", "JPSS1", "JPSS2", "AQUA")
fig3_list <- Filter(Negate(is.null),
                    lapply(sat_names_fig3, function(s) load_fig3_file(match_site, s, match_date)))

# HYPERNETS from the first available file; satellite bands from all files
hyp_spectrum <- fig3_list[[1]] |> filter(sensor == "HYPERNETS")
sat_bands <- map_dfr(fig3_list, ~ filter(.x, sensor != "HYPERNETS")) |>
  mutate(sensor = case_when(
    sensor == "PACE"  ~ "PACE OCI",
    sensor == "S3A"   ~ "S3A OLCI",
    sensor == "S3B"   ~ "S3B OLCI",
    sensor == "SNPP"  ~ "VIIRS SNPP",
    sensor == "JPSS1" ~ "VIIRS JPSS-1",
    sensor == "JPSS2" ~ "VIIRS JPSS-2",
    sensor == "AQUA"  ~ "MODIS Aqua",
    TRUE ~ sensor
  ))

sat_levels <- sort(unique(sat_bands$sensor))
sensor_colours <- c(
  "HYPERNETS" = "black",
  setNames(RColorBrewer::brewer.pal(max(3, length(sat_levels)), "Dark2")[seq_len(length(sat_levels))],
           sat_levels)
)

pl_fig3 <- ggplot(data = hyp_spectrum, aes(x = wavelength, y = rhow, colour = sensor)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(data = sat_bands, aes(colour = sensor), size = 3) +
  geom_errorbar(data = filter(sat_bands, !is.na(std_min)),
                aes(ymin = std_min, ymax = std_max, colour = sensor), width = 5) +
  scale_colour_manual(values = sensor_colours) +
  coord_cartesian(xlim = c(380, 900)) +
  labs(x = "Wavelength (nm)",
       y = "<i>ρ<sub>w</sub></i>",
       colour = "Sensor",
       caption = paste0(match_site, "  ·  ", match_date)) +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_markdown(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.background = element_rect(colour = "black", fill = "white"),
        plot.caption = element_text(size = 10, hjust = 0))
# pl_fig3
ggsave("figures/fig_3.png", pl_fig3, width = 10, height = 6)


# Figure 4 ----------------------------------------------------------------
# HYPERNETS vs satellite scatterplots for the single representative matchup
# date defined in match_date above. One panel per satellite found on that date.
# Uses load_matchup_long() which internally calls load_matchup_mean(), so all
# sensor-family data_type quirks (PACE, JPSS1, etc.) are handled correctly.

load_fig4_file <- function(site, sat_name, date_str) {
  folder <- file_path_build(site, sat_name)
  files  <- list.files(folder, pattern = date_str, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  tryCatch(load_matchup_long(files[1]), error = function(e) NULL)
}

fig4_list <- Filter(Negate(is.null),
                    lapply(sat_names_fig3, function(s) load_fig4_file(match_site, s, match_date)))

plot_fig4_panel <- function(df) {
  sat_col    <- setdiff(names(df), c("file_name", "wavelength", "Hyp"))[1]
  colour_pal <- colour_nm_func(sat_col)

  # PACE has continuous wavelengths; group into the same bands colour_nm_func uses.
  # All other sensors have discrete bands that match the palette names directly.
  if (sat_col == "PACE") {
    breaks <- c(350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1050)
    df <- df |> mutate(wl_col = cut(wavelength, breaks = breaks,
                                    labels = names(colour_pal), include.lowest = TRUE))
  } else {
    df <- df |> mutate(wl_col = factor(as.character(wavelength), levels = names(colour_pal)))
  }

  max_val  <- max(c(df$Hyp, df[[sat_col]]), na.rm = TRUE) * 1.05
  df_stats <- base_stats(df$Hyp, df[[sat_col]])

  ggplot(df, aes(x = Hyp, y = .data[[sat_col]], colour = wl_col)) +
    geom_abline(slope = 1, intercept = 0, colour = "black") +
    geom_point(size = 3, alpha = 0.8) +
    annotate("text", x = 0, y = max_val, hjust = 0, vjust = 1, size = 3.5,
             label = paste0("n = ", df_stats$n,
                            "\nS = ", sprintf("%.2f", df_stats$Slope_II),
                            "\nBias = ", sprintf("%.1f", df_stats$Bias_50), "%",
                            "\nError = ", sprintf("%.1f", df_stats$Error_50), "%")) +
    scale_colour_manual(values = colour_pal, drop = FALSE) +
    coord_fixed(xlim = c(0, max_val), ylim = c(0, max_val)) +
    labs(x = "HYPERNETS *ρ*<sub>w</sub>",
         y = paste0(sat_col, " *ρ*<sub>w</sub>"),
         colour = "Wavelength (nm)") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          axis.title.x = element_markdown(size = 11),
          axis.title.y = element_markdown(size = 11),
          axis.text = element_text(size = 9),
          legend.position = "bottom",
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))
}

fig4_panels <- lapply(fig4_list, plot_fig4_panel)
fig_4 <- wrap_plots(fig4_panels, ncol = 2) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave("figures/fig_4.png", fig_4, width = 12,
       height = ceiling(length(fig4_panels) / 2) * 6)


# Figure 11 ---------------------------------------------------------------

## Load global matchups
df_matchups_global <- read_csv("output/global_stats_all.csv", show_col_types = FALSE) |>
  filter(sensor_X == "HYPERNETS", sensor_Y != "HYPERNETS") |>
  rename(Error = Error_50) |>
  mutate(sensor_Y = case_when(sensor_Y == "AQUA"   ~ "Aqua",
                              sensor_Y == "S3A"    ~ "S3A",
                              sensor_Y == "S3B"    ~ "S3B",
                              sensor_Y == "S3_all" ~ "S3 all",
                              sensor_Y == "PACE"   ~ "PACE",
                              sensor_Y == "SNPP"   ~ "SNPP",
                              sensor_Y == "JPSS1"  ~ "JPSS1",
                              sensor_Y == "JPSS2"  ~ "JPSS2"))

# Prep data for matrix plots
df_matchups_global_pretty <- df_matchups_global |>
  filter(wavelength >= 400, wavelength <= 600) |>
  mutate(sensor_sat = case_when(sensor_Y %in% c("Aqua")                   ~ "MODIS",
                                sensor_Y %in% c("S3A", "S3B", "S3 all")   ~ "OLCI",
                                sensor_Y %in% c("PACE")                    ~ "OCI",
                                sensor_Y %in% c("SNPP", "JPSS1", "JPSS2") ~ "VIIRS"),
         wavelength_clean = case_when(sensor_sat == "VIIRS" & wavelength %in% c(410, 411) ~ "410/411",
                                      sensor_sat == "VIIRS" & wavelength %in% c(443, 445) ~ "443/445",
                                      sensor_sat == "VIIRS" & wavelength %in% c(486, 489) ~ "486/489",
                                      sensor_sat == "VIIRS" & wavelength %in% c(551, 556) ~ "551/556",
                                      TRUE ~ as.character(wavelength))) |>
  mutate(sensor_Y   = factor(sensor_Y,   levels = c("Aqua", "S3 all", "S3B", "S3A",
                                                     "PACE", "JPSS2", "JPSS1", "SNPP")),
         sensor_sat = factor(sensor_sat, levels = c("MODIS", "OLCI", "OCI", "VIIRS")))

# Matrix plot
plot_matrix_error <- function(df, val_range) {
  df_all <- df |>
    group_by(sensor_Y, sensor_sat) |>
    summarise(Error = mean(Error, na.rm = TRUE), .groups = "drop") |>
    mutate(wavelength_clean = "All")
  df <- bind_rows(df, df_all)
  df_round <- df |>
    mutate(Error = case_when(Error > val_range[2] ~ val_range[2], TRUE ~ Error))
  ggplot(data = df_round, aes(x = wavelength_clean, y = sensor_Y)) +
    geom_tile(aes(fill = Error), colour = "black") +
    geom_label(data = df, aes(label = sprintf("%.1f", round(Error, 1))), size = 3) +
    labs(x = "Waveband (nm)", y = NULL, fill = "Error (%)") +
    facet_grid(sensor_sat~., scales = "free") +
    scale_fill_viridis_c(limits = val_range, breaks = c(10, 20, 30), labels = c("10", "20", "30")) +
    coord_cartesian(expand = FALSE) +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_text(size = 14),
          legend.text  = element_text(size = 12),
          axis.title.x = element_markdown(size = 12),
          axis.title.y = element_markdown(size = 12),
          axis.text    = element_text(size = 10))
}

# Build one panel per sensor family and stack
sensors <- c("MODIS", "VIIRS", "OLCI", "OCI")
plots_by_sensor <- purrr::set_names(
  purrr::map(sensors, function(s) {
    val_range <- c(0, 40)
    df_sub <- df_matchups_global_pretty |> dplyr::filter(sensor_sat == s)
    plot_matrix_error(df_sub, val_range)
  }),
  sensors
)

fig_11 <- plots_by_sensor$MODIS / plots_by_sensor$VIIRS / plots_by_sensor$OLCI / plots_by_sensor$OCI +
  plot_annotation(tag_levels = "a", tag_suffix = ")") +
  patchwork::plot_layout(guides = "collect", axis_titles = "collect", heights = c(0.35, 1, 1, 1))
ggsave("figures/fig_11.png", fig_11, width = 12, height = 9)


# Fig S1 ------------------------------------------------------------------

# Barplot of PACE OCI Error and Bias across all wavelengths (HYPERNETS vs PACE at MAFR)
global_stats_OCI <- read_csv("output/global_stats_RHOW_OCI.csv", show_col_types = FALSE) |>
  filter(sensor_X == "HYPERNETS", sensor_Y == "PACE")

theme_S1 <- theme(panel.border  = element_rect(fill = NA, color = "black"),
                  axis.title.x  = element_markdown(size = 12),
                  axis.title.y  = element_markdown(size = 12),
                  axis.text     = element_text(size = 10))

pl_Error_OCI <- ggplot(global_stats_OCI, aes(x = wavelength, y = Error_50)) +
  geom_col(fill = "steelblue") +
  labs(x = "Wavelength (nm)", y = "Error (%)") +
  theme_S1

pl_Bias_OCI <- ggplot(global_stats_OCI, aes(x = wavelength, y = Bias_50)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(x = "Wavelength (nm)", y = "Bias (%)") +
  theme_S1

fig_S1 <- pl_Error_OCI / pl_Bias_OCI + plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave("figures/fig_S1.png", fig_S1, width = 9, height = 8)

