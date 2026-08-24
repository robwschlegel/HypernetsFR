# code/4_matchups_global.R
# Aggregate per-file match-ups into global (per-wavelength) statistics, once
# outlier screening (2_outliers.R) and the matchup-protocol sensitivity checks
# (3_sensitivity.R) have both been run.


# Setup -------------------------------------------------------------------

source("code/0_functions.R")

# NB: this script assumes code/1_matchups_single.R, code/2_outliers.R, and
# code/3_sensitivity.R have already been run and their outputs already exist
# on disk), since global_stats() reads both.


# Global statistics --------------------------------------------------------

process_sensor("MODIS", "global", daily_average = TRUE)
process_sensor("VIIRS", "global", daily_average = TRUE)
process_sensor("OLCI", "global", daily_average = TRUE)
process_sensor("OCI", "global", daily_average = TRUE)


# Combine all outputs ----------------------------------------------------

# Load outlier reports
outliers_sat <- read_csv("meta/satellite_outliers.csv", show_col_types = FALSE) |> distinct()

# Decide which wavelengths to compare with wavebands
wave_length_bands <- c(380, 400, 412, 443, 490, 510, 560, 620, 673, 700, 885, 1050)

# Load PACE data and subset by wavelengths of interest (see wave_length_bands above)
global_stats_wavelengths <- read_csv("output/global_stats_RHOW_OCI.csv", show_col_types = FALSE) |>
  filter(wavelength %in% wave_length_bands)

# Load all global stats
# NB: Careful with the exact indexing of files here
global_stats_all <- map_dfr(dir("output", pattern = "global", full.names = TRUE)[c(2, 4:5)],
                                read_csv, show_col_types = FALSE) |>
  bind_rows(global_stats_wavelengths)
write_csv(global_stats_all, file = "output/global_stats_all.csv")

# Stop here if running from terminal
if (!interactive()) quit(save = "no", status = 0)


# Investigate results ----------------------------------------------------

# Get matchups counts, outliers, etc.
global_count_var_name <- global_stats_all |>
  filter(sensor_X %in% c("HYPERNETS")) |>
  filter(var_name == "RHOW") |>
  dplyr::select(sensor_X, sensor_Y, var_name, n_clean, n_no_out) |>
  group_by(var_name, sensor_X, sensor_Y) |>
  filter(n_clean == max(n_clean, na.rm = TRUE)) |>
  ungroup() |>
  distinct()

# Quantify which sensors matched most closely to which satellites
global_match_mean <- global_stats_all |>
  filter(sensor_X %in% c("HYPERNETS")) |>
  filter(var_name == "RHOW") |>
  filter(wavelength >= 400, wavelength <= 600) |>
  summarise(Slope = mean(Slope_II, na.rm = TRUE),
            Bias_mean = mean(Bias_50, na.rm = TRUE),
            Bias_abs = mean(abs(Bias_50), na.rm = TRUE),
            Error = mean(Error_50, na.rm = TRUE), .by = c("site_name", "var_name", "sensor_X", "sensor_Y"))
global_match_mean_red <- global_stats_all |>
  filter(sensor_X %in% c("HYPERNETS", "TRIOS", "HYPERPRO")) |>
  filter(var_name == "RHOW") |>
  filter(wavelength > 600) |>
  summarise(Slope = mean(Slope_II, na.rm = TRUE),
            MRD_mean = mean(MRD_50, na.rm = TRUE),
            MRD_abs = mean(abs(MRD_50), na.rm = TRUE),
            MARD = mean(MARD_50, na.rm = TRUE), .by = c("site_name", "var_name", "sensor_X", "sensor_Y"))

# Mean just for in situ systems
global_match_is_mean <- global_stats_all |>
  filter(sensor_X %in% c("HYPERNETS")) |>
  filter(var_name == "RHOW") |>
  filter(wavelength >= 400, wavelength <= 600) |>
  # filter(!(sensor_Y %in% c("HYPERNETS", "TRIOS", "HYPERPRO"))) |>
  summarise(Slope = mean(Slope_II, na.rm = TRUE),
            Bias_mean = mean(Bias_50, na.rm = TRUE),
            Bias_abs = mean(abs(Bias_50), na.rm = TRUE),
            Error = mean(Error_50, na.rm = TRUE), .by = c("site_name", "var_name", "sensor_X"))
global_match_is_mean_red <- global_stats_all |>
  filter(sensor_X %in% c("HYPERNETS")) |>
  filter(var_name == "RHOW") |>
  filter(wavelength > 600) |>
  # filter(sensor_Y != "S3B") |>
  # filter(!(sensor_Y %in% c("HYPERNETS", "TRIOS", "HYPERPRO"))) |>
  summarise(Slope = mean(Slope_II, na.rm = TRUE),
            MRD_mean = mean(MRD_50, na.rm = TRUE),
            MRD_abs = mean(abs(MRD_50), na.rm = TRUE),
            MARD = mean(MARD_50, na.rm = TRUE), .by = c("site_name", "var_name", "sensor_X"))

# Global mean matchups from the perspective of the satellites
global_match_sat_mean <- global_stats_all |>
  filter(sensor_X %in% c("HYPERNETS")) |>
  filter(var_name == "RHOW") |>
  filter(wavelength >= 400, wavelength <= 600) |>
  # filter(!(sensor_Y %in% c("HYPERNETS", "TRIOS", "HYPERPRO"))) |>
  summarise(Slope = mean(Slope_II, na.rm = TRUE),
            Bias_mean = mean(Bias_50, na.rm = TRUE),
            Bias_abs = mean(abs(Bias_50), na.rm = TRUE),
            Error = mean(Error_50, na.rm = TRUE), .by = c("sensor_Y"))
global_match_sat_mean_red <- global_stats_all |>
  filter(sensor_X %in% c("HYPERNETS")) |>
  filter(var_name == "RHOW") |>
  filter(wavelength > 600) |>
  # filter(!(sensor_Y %in% c("HYPERNETS", "TRIOS", "HYPERPRO"))) |>
  summarise(Slope = mean(Slope_II, na.rm = TRUE),
            MRD_mean = mean(MRD_50, na.rm = TRUE),
            MRD_abs = mean(abs(MRD_50), na.rm = TRUE),
            MARD = mean(MARD_50, na.rm = TRUE), .by = c("sensor_Y"))

# Count of number of wavebands with negative or positive biases
global_match_bias_sign <- global_stats_all |>
  filter(sensor_X %in% c("HYPERNETS")) |>
  # filter(var_name == "RHOW") |>
  mutate(Bias_positive = case_when(Bias_50 > 0 ~ 1, TRUE ~ 0),
         Bias_negative = case_when(Bias_50 < 0 ~ 1, TRUE ~ 0)) |>
  summarise(Bias_positive = sum(Bias_positive),
            Bias_negative = sum(Bias_negative), .by = c("site_name", "var_name", "sensor_X", "sensor_Y"))

# Summarise per in situ platform against satellites
global_match_bias_sign_remote <- global_match_bias_sign |>
  filter(!(sensor_Y %in% c("HYPERNETS", "TRIOS", "HYPERPRO"))) |>
  filter(var_name == "RHOW") |>
  summarise(Bias_positive = sum(Bias_positive),
            Bias_negative = sum(Bias_negative), .by = c("sensor_X"))

# Average bias and error values per waveband
global_waveband_mean <- global_stats_all |>
  filter(sensor_X %in% c("HYPERNETS")) |>
  # filter(!(sensor_Y %in% c("HYPERNETS", "TRIOS", "HYPERPRO"))) |>
  filter(var_name == "RHOW") |>
  summarise(Slope = mean(Slope_II, na.rm = TRUE),
            Bias_mean = mean(Bias_50, na.rm = TRUE),
            Bias_abs = mean(abs(Bias_50), na.rm = TRUE),
            Error = mean(Error_50, na.rm = TRUE), .by = c("wavelength"))
global_waveband_mean_red <- global_stats_all |>
  filter(sensor_X %in% c("HYPERNETS")) |>
  # filter(!(sensor_Y %in% c("HYPERNETS", "TRIOS", "HYPERPRO"))) |>
  # filter(sensor_Y != "S3B") |>
  filter(var_name == "RHOW") |>
  filter(wavelength > 600) |>
  summarise(Slope = mean(Slope_II, na.rm = TRUE),
            MRD_mean = mean(MRD_50, na.rm = TRUE),
            MRD_abs = mean(abs(MRD_50), na.rm = TRUE),
            MARD = mean(MARD_50, na.rm = TRUE), .by = c("wavelength"))

# Plot the global mean matchups for HYPERNETS vs. satellites, per site and sensor
global_stats_all |>
  filter(sensor_X %in% c("HYPERNETS")) |>
  ggplot(aes(x = sensor_Y, y = Error_50)) +
  geom_boxplot(aes(fill = sensor_Y)) +
  labs(x = "In situ platform", y = "Error (%)", fill = "Satellite sensor") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = "black"))


# Check individual matchups -----------------------------------------------
# Ad hoc diagnostic template for inspecting one matchup file in detail.
# Change site_name, sat_name, and file_name below to target any file of interest.

site_name <- "MAFR"
sat_name  <- "S3A"
file_name <- "S3A_20240531T101658_vs_HYPERNETS_20240531T100000_RHOW_C.csv"

# Load file
match_1 <- load_matchup_long(file.path(file_path_build(site_name, sat_name), file_name))

# Create vectors from filtered columns
(x_vec <- match_1[["Hyp"]])
(y_vec <- match_1[[sat_name]])

# Calculate stats one-by-one
message("Slope : ", round(coef(lm(y_vec ~ x_vec))[2], 4))
message("RMSE : ", round(sqrt(mean((y_vec - x_vec)^2, na.rm = TRUE)), 4))
message("MSA : ", round(mean(abs(y_vec - x_vec), na.rm = TRUE), 4))
message("MAPE (%) : ", round(mean(abs((y_vec - x_vec) / x_vec), na.rm = TRUE) * 100, 2))

# Calculate Bias and Error (Pahlevan's method)
log_ratio <- log10(y_vec / x_vec)
bias_pahlevan <- median(log_ratio, na.rm = TRUE)
bias_pahlevan_final <- sign(bias_pahlevan) * (10^abs(bias_pahlevan) - 1)
message("Bias (%) : ", round(bias_pahlevan_final * 100, 2))

error_pahlevan <- median(abs(log_ratio), na.rm = TRUE)
error_pahlevan_final <- 10^error_pahlevan - 1
message("Error (%) : ", round(error_pahlevan_final * 100, 2))

# Plot data
match_1 |>
  ggplot(aes(x = Hyp, y = .data[[sat_name]])) +
  geom_point(aes(colour = wavelength))


# Check individual global values ------------------------------------------
# Ad hoc diagnostic template for inspecting per-wavelength stats across all
# QC-passed files for a given satellite. Change the three variables below to
# target any site/satellite combination.

site_name <- "MAFR"
sat_name  <- "S3A"
sensor_Z  <- "OLCI"   # sensor family: OLCI / MODIS / VIIRS / OCI

# Build file list filtered to files that passed the spatiotemporal QC gate
folder_path <- file_path_build(site_name, sat_name)
file_list <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)
match_base_details <- read_csv(paste0("output/matchup_stats_RHOW_", sensor_Z, ".csv"),
                               show_col_types = FALSE) |>
  dplyr::select(file_name) |> distinct()
file_list_clean <- file_list[basename(file_list) %in% match_base_details$file_name]

# Load data
match_base <- map_dfr(file_list_clean, load_matchup_long)
# print(unique(match_base$wavelength))

# Choose wavelengths for this satellite
(W_nm <- W_nm_out(sat_name))

# Get data.frame for matchup based on the wavelength of choice
(matchup_filt <- filter(match_base, wavelength == W_nm[10]))

# Create vectors from filtered columns (x = HYPERNETS, y = satellite)
(x_vec <- matchup_filt[["Hyp"]])
(y_vec <- matchup_filt[[sat_name]])

# Calculate stats one-by-one
message("Slope : ", round(coef(lm(y_vec ~ x_vec))[2], 4))
message("RMSE : ", round(sqrt(mean((y_vec - x_vec)^2, na.rm = TRUE)), 4))
message("MSA : ", round(mean(abs(y_vec - x_vec), na.rm = TRUE), 4))
message("MAPE (%) : ", round(mean(abs((y_vec - x_vec) / x_vec), na.rm = TRUE) * 100, 2))

# Calculate Bias and Error (Pahlevan's method)
log_ratio <- log10(y_vec / x_vec)
bias_pahlevan <- median(log_ratio, na.rm = TRUE)
bias_pahlevan_final <- sign(bias_pahlevan) * (10^abs(bias_pahlevan) - 1)
message("Bias (%) : ", round(bias_pahlevan_final * 100, 2))

error_pahlevan <- median(abs(log_ratio), na.rm = TRUE)
error_pahlevan_final <- 10^error_pahlevan - 1
message("Error (%) : ", round(error_pahlevan_final * 100, 2))

# Plot data
match_base |> filter(wavelength == W_nm[8]) |>
  ggplot(aes(x = Hyp, y = .data[[sat_name]])) +
  geom_point(aes(colour = wavelength))

