# code/3_sensitivity.R
# Sensitivity analyses used to justify the match-up protocol's time-window and
# distance-window choices (see manuscript/roadmap.md, "site-specific matchup
# criteria" item, and site_diff_time_limit()/daily_average_matchups() in
# code/0_functions.R).
#
# NB: unlike Doxaran et al. 2024 (who used a fixed 3x3-pixel-box/+-30 min window at
# the clear-ish Berre lagoon and a nearest-pixel/+-15 min window at turbid Gironde),
# this pipeline uses ONE distance ceiling (dist_limit = 10 km; raised from 5 km on
# 2026-07-14 to accommodate a known THFR PACE pixel-extraction offset -- see
# manuscript/upstream-data-bugs.md Bug 9) for every site, but a SITE-SPECIFIC time
# window (site_diff_time_limit() in
# code/0_functions.R: 15 min at MAFR, 30 min at THFR). This script empirically checks
# both choices, and is also where the eventual before/after comparison of
# daily_average_matchups() (see code/0_functions.R) should live once that first-draft
# function has been validated.
#
# This script is intended to be run AFTER 1_matchups_single.R (needs
# output/matchup_stats_RHOW_*.csv and output/matchup_noQC_stats_RHOW_*.csv) and is
# read by no other script, so it can safely be skipped in a quick production run once
# its conclusions have been folded into the Methods text -- but re-run it whenever
# THFR data change materially, since its whole purpose is to justify the window
# choices with real data.


# Setup -------------------------------------------------------------------

source("code/0_functions.R")


# Distance sanity check ------------------------------------------------------

# Load all individual matchup stats across sensor families, both QC-passed and not,
# so the full range of observed distances (not just the already-filtered subset) is visible
matchup_all_noQC <- map_dfr(dir("output", pattern = "matchup_stats_RHOW_|matchup_noQC_stats_RHOW_", full.names = TRUE),
                            read_csv, show_col_types = FALSE)

# Confirm nearest-pixel distances relative to the 10 km ceiling.
# NB: dist_limit was raised from 5 km to 10 km on 2026-07-14 to accommodate a systematic
# pixel-extraction offset in THFR PACE data (Bug 9 in manuscript/upstream-data-bugs.md).
# For MAFR and all sensors except THFR PACE, distances remain well under 1 km in practice.
dist_summary <- matchup_all_noQC |>
  filter(sensor_X == "Hyp") |>
  summarise(dist_min = min(dist, na.rm = TRUE),
            dist_median = median(dist, na.rm = TRUE),
            dist_p95 = quantile(dist, 0.95, na.rm = TRUE),
            dist_max = max(dist, na.rm = TRUE),
            n_over_1km = sum(dist > 1, na.rm = TRUE),
            n_over_5km  = sum(dist > 5,  na.rm = TRUE),
            n_over_10km = sum(dist > 10, na.rm = TRUE),
            n = dplyr::n(),
            .by = c("site_name", "sensor_Y"))
print(dist_summary)

# Visualise the distance distribution per site/sensor, with the 5 km ceiling marked
pl_dist <- matchup_all_noQC |>
  filter(sensor_X == "Hyp") |>
  ggplot(aes(x = dist)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = 10, colour = "red", linetype = "dashed") +
  geom_vline(xintercept = 5,  colour = "orange", linetype = "dotted") +
  labs(x = "Distance between HYPERNETS station and nearest satellite pixel (km)",
       y = "Count",
       title = "Distance check (10 km ceiling in red; original 5 km threshold in orange)") +
  facet_grid(site_name ~ sensor_Y, scales = "free_y") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = "black"))
ggsave("figures/sensitivity_distance_check.png", pl_dist, width = 12, height = 4)


# Time-window sensitivity -----------------------------------------------------

# For each site, look at how Error_50/Bias_50 vary with diff_time, to justify the
# site-specific time-window choice (site_diff_time_limit(): MAFR = 15 min, THFR = 30 min).
# NB: follows the same approach as the Tara "in review" paper's time-window check (which found
# no significant trend there); here a site-dependent trend is plausible given MAFR's much
# faster tidal turbidity dynamics relative to THFR's more stable lagoon water -- that is
# exactly the asymmetry the two different site-specific limits are meant to capture.
pl_time_sensitivity <- matchup_all_noQC |>
  filter(sensor_X == "Hyp") |>
  filter(Error_50 < 1000) |> # TODO: Loop back and look at these high error values
  ggplot(aes(x = diff_time, y = Error_50)) +
  geom_point(aes(colour = dist), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_vline(data = data.frame(site_name = c("MAFR", "THFR"), time_limit = c(15, 30)),
             aes(xintercept = time_limit), colour = "red", linetype = "dashed") +
  scale_colour_viridis_c() +
  labs(x = "Time difference between HYPERNETS scan and satellite overpass (minutes)",
       y = "Error (%)", colour = "Distance\n(km)",
       title = "Time-window sensitivity per site (site-specific limit shown in red)") +
  facet_grid(site_name ~ sensor_Y, scales = "free") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = "black"))
ggsave("figures/sensitivity_time_window.png", pl_time_sensitivity, width = 12, height = 8)

# Formal test: does Error_50 trend significantly with diff_time, per site and sensor?
# (mirrors the check performed in the Tara "in review" paper's Discussion)
time_trend_test <- matchup_all_noQC |>
  filter(sensor_X == "Hyp") |>
  nest_by(site_name, sensor_Y) |>
  mutate(model = list(lm(Error_50 ~ diff_time, data = data)),
         slope = coef(model)[["diff_time"]],
         p_value = summary(model)$coefficients["diff_time", "Pr(>|t|)"]) |>
  dplyr::select(site_name, sensor_Y, slope, p_value)
print(time_trend_test)


# Daily-averaging sensitivity (site-specific time window) -----------------------

# TODO: once daily_average_matchups() (code/0_functions.R) has been validated against real
# multi-day MAFR/THFR data, compare global_stats(..., daily_average = FALSE) against
# global_stats(..., daily_average = TRUE) here, per site and sensor family, to quantify how
# much the day-averaging step changes the headline Error/Bias values, and to check whether it
# meaningfully addresses the "not all matchups are independent of one another" caveat raised in
# the Tara "in review" paper's Conclusion. Not yet run -- needs THFR data, and confirmation
# that daily_average_matchups()'s file-name date-parsing is robust across all four sensor
# families' naming conventions (see CLAUDE.md's note on satellite-family naming inconsistencies).
#
compare_daily_avg <- function(site_name, sensor_Y){
  no_avg <- global_stats(site_name, sensor_Y, daily_average = FALSE) |> mutate(daily_average = FALSE)
  avg    <- global_stats(site_name, sensor_Y, daily_average = TRUE)  |> mutate(daily_average = TRUE)
  bind_rows(no_avg, avg)
}
daily_avg_comparison <- purrr::pmap_dfr(sensor_grid("OLCI"), compare_daily_avg)
pl_daily_avg <- daily_avg_comparison |>
  filter(sensor_X == "HYPERNETS") |>
  ggplot(aes(x = wavelength, y = Error_50, colour = daily_average, size = n_w_nm_clean)) +
  geom_point() +
  facet_grid(site_name ~ sensor_Y, scales = "free_y") +
  labs(title = "Effect of daily averaging on Error (%) per wavelength", size = "Unique \n data points")
ggsave("figures/sensitivity_daily_average.png", pl_daily_avg, width = 12, height = 8)

