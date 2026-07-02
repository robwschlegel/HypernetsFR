# code/matchups.R
# Get the stats for all matchups and visualise results


# Setup -------------------------------------------------------------------

source("code/functions.R")


# Individual matchup stats ------------------------------------------------

# process_sensor("MODIS") # To add
# process_sensor("VIIRS")
process_sensor("OLCI")
# process_sensor("OCI")

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


# Outliers ---------------------------------------------------------------

# Run this to populate the necessary files used in the global matchups
# NB: The outlier identification process is more manually focussed, hence the separate script
source("code/outliers.R")


# Global statistics --------------------------------------------------------

# NB: Run code/outliers.R before the global stats in order to filter outliers

# process_sensor("MODIS", "global") # To add
# process_sensor("VIIRS", "global")
process_sensor("OLCI", "global")
# process_sensor("OCI", "global")

# Load outlier reports
outliers_sat <- read_csv("meta/satellite_outliers.csv", show_col_types = FALSE) |> distinct()

# Decide which wavelengths to compare with wavebands
wave_length_bands <- c(380, 400, 412, 443, 490, 510, 560, 620, 673, 700, 885, 1050)

# TODO: Add back in once PACE data are available
# Load PACE data separately to filter specific wavebands
# NB: Careful with the exact indexing of files here
# global_stats_wavelengths <- map_dfr(dir("output", pattern = "global", full.names = TRUE)[c(4)], read_csv, show_col_types = FALSE) |> 
#   filter(Wavelength_nm %in% wave_length_bands)

# Load all global stats
# NB: Careful with the exact indexing of files here
global_stats_all <- map_dfr(dir("output", pattern = "global", full.names = TRUE)[c(2)], 
                                read_csv, show_col_types = FALSE) #|> 
  # bind_rows(global_stats_wavelengths)
write_csv(global_stats_all, file = "output/global_stats_all.csv")

# Visualise difference between linear-space and log-space slopes
# global_stats_all |> 
#   filter(sensor_X %in% c("HYPERNETS")) |> 
#   ggplot() +
#   geom_histogram(aes(x = Slope), colour = "green", alpha = 0.3, binwidth = 0.1) +
#   geom_histogram(aes(x = Slope_log), colour = "red", alpha = 0.3, binwidth = 0.1) +
#   facet_grid(sensor_X ~ sensor_Y)

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

# Plot the global mean matchups per in situ sensor
global_match_mean |> 
  filter(!(sensor_Y %in% c("HYPERNETS"))) |> 
  ggplot(aes(x = sensor_X, y = Error)) +
  geom_boxplot(aes(fill = sensor_X))


# Check individual matchups -----------------------------------------------

# Files of interest
# HYPERNETS_vs_TRIOS_vs_20240816T081000_LW.csv # Negative values cause massive MAPE (%)

# Load file
match_1 <- load_matchup_long(paste0(file_path_build("LW", "HYPERNETS", "TRIOS"), "HYPERNETS_vs_TRIOS_vs_20240816T081000_LW.csv"))

# Create vectors from filtered columns
(x_vec <- match_1[["Hyp"]])
(y_vec <- match_1[["TRIOS"]])

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
match_base |> 
  ggplot(aes(x = Hyp, y = TRIOS)) +
  geom_point(aes(colour = wavelength))


# Check individual global values ------------------------------------------

# Choose acordingly
folder_path <- file_path_build("RHOW", "HYPERNETS", "HYPERPRO")

# List all files in directory
file_list <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)
match_base_details <- read_csv("~/HypernetsTara/output/matchup_stats_RHOW_in_situ.csv") |>
  dplyr::select(file_name) |> distinct()
file_list_clean <- file_list[basename(file_list) %in% match_base_details$file_name]

# Load data
match_base <- map_dfr(file_list_clean, load_matchup_long)
# print(unique(match_base$wavelength))

# Auto-generate chosen wavelengths
(W_nm <- W_nm_out("HYPERNETS"))
# W_nm <- c(380, 400, 412, 443, 490, 510, 560, 620, 673, 700)

# Get data.frame for matchup based on the wavelength of choice
(matchup_filt <- filter(match_base, wavelength == W_nm[10]))

# Create vectors from filtered columns
(x_vec <- matchup_filt[["HYPERPRO"]])
(y_vec <- matchup_filt[["Hyp"]])

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
  ggplot(aes(x = HYPERPRO, y = Hyp)) +
  geom_point(aes(colour = wavelength))

