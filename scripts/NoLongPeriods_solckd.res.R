library(dplyr)
library(ggplot2)
library(broom)

# ==============================================================================
# 0. SETUP & DATA PREPARATION
# ==============================================================================
# Assuming 'Ocean_solar_anom' contains:
#   - dt.mnth: monthly time index (e.g., 1850.0, 1850.083, ...)
#   - Anomaly: global ocean SST anomaly
#   - resids_solar: residuals from your solar-phase locked model

# Ensure time step in months for harmonic frequency calculations
dt_start <- min(Ocean_solar_anom$dt.mnth)
Ocean_solar_anom <- Ocean_solar_anom %>%
  mutate(time_months = (dt.mnth - dt_start) * 12)

# ==============================================================================
# 1. SPECTRAL DIAGNOSTICS OF RESIDUALS
# ==============================================================================
# A. Standard Periodogram (Check if low frequencies slope upward or stay flat)
spec_pgram <- spec.pgram(Ocean_solar_anom$resids_solar, plot = FALSE, log = "yes")
df_pgram <- data.frame(
  period_years = (1 / spec_pgram$freq) / 12,
  spec = spec_pgram$spec
)

ggplot(df_pgram %>% filter(period_years <= 176), aes(x = period_years, y = spec)) +
  geom_line(color = "steelblue") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Residual Periodogram (Log-Log)",
       x = "Period (Years)", y = "Spectral Density") +
  theme_minimal()

# B. Maximum Entropy / AR Spectrum (Superior resolution for short records)
spec_ar_fit <- spec.ar(Ocean_solar_anom$resids_solar, plot = FALSE, log = "yes")
df_ar <- data.frame(
  period_years = (1 / spec_ar_fit$freq) / 12,
  spec = spec_ar_fit$spec
)

# Plot AR spectrum focusing on long periods
ggplot(df_ar %>% filter(period_years >= 10), aes(x = period_years, y = spec)) +
  geom_line(color = "darkred", size = 1) +
  scale_x_log10() +
  labs(title = "AR (Maximum Entropy) Spectrum of Residuals",
       subtitle = "Flatness at high periods confirms complete low-frequency capture",
       x = "Period (Years)", y = "Power") +
  theme_minimal()

# ==============================================================================
# 2. FIT KNOWN ASTROPHYSICAL HARMONICS
# ==============================================================================
# Define periods in years converted to months
T_gleissberg <- 88 * 12     # ~88 years
T_suess      <- 208 * 12    # ~208 years
T_eddy       <- 1000 * 12   # ~1000 years (Millennial scale)

# Construct harmonic predictors
Ocean_solar_anom <- Ocean_solar_anom %>%
  mutate(
    sin_gleiss = sin(2 * pi * time_months / T_gleissberg),
    cos_gleiss = cos(2 * pi * time_months / T_gleissberg),
    sin_suess  = sin(2 * pi * time_months / T_suess),
    cos_suess  = cos(2 * pi * time_months / T_suess),
    sin_eddy   = sin(2 * pi * time_months / T_eddy),
    cos_eddy   = cos(2 * pi * time_months / T_eddy)
  )

# Fit harmonic model to the residual series or total low-frequency baseline
fit_astro_resids <- lm(
  resids_solar ~ sin_gleiss + cos_gleiss +
    sin_suess  + cos_suess  +
    sin_eddy   + cos_eddy,
  data = Ocean_solar_anom
)

# View statistical significance of candidate long periods
summary(fit_astro_resids)

# Calculate variance explained by long periods in the residuals
var_total_resids <- var(Ocean_solar_anom$resids_solar)
var_astro_fitted <- var(fitted(fit_astro_resids))
leakage_ratio    <- var_astro_fitted / var_total_resids

cat(sprintf("Low-frequency variance leakage in residuals: %.2f%%\n", leakage_ratio * 100))

# ==============================================================================
# 3. BACKCASTING TO THE MAUNDER MINIMUM (~1650 - 2026)
# ==============================================================================
# Fit the harmonic terms directly on the LOW-FREQUENCY BASELINE (Trend)
fit_astro_baseline <- lm(
  Anomaly ~ sin_gleiss + cos_gleiss +
    sin_suess  + cos_suess  +
    sin_eddy   + cos_eddy,
  data = Ocean_solar_anom
)

# Create an extended timeline back to 1645
df_backcast <- data.frame(
  dt.mnth = seq(from = 1645, to = 2026, by = 1/12)
) %>%
  mutate(
    time_months = (dt.mnth - dt_start) * 12,
    sin_gleiss  = sin(2 * pi * time_months / T_gleissberg),
    cos_gleiss  = cos(2 * pi * time_months / T_gleissberg),
    sin_suess   = sin(2 * pi * time_months / T_suess),
    cos_suess   = cos(2 * pi * time_months / T_suess),
    sin_eddy    = sin(2 * pi * time_months / T_eddy),
    cos_eddy    = cos(2 * pi * time_months / T_eddy)
  )

# Generate backcasted baseline values
df_backcast$baseline_backcast <- predict(fit_astro_baseline, newdata = df_backcast)

# ==============================================================================
# 4. PLOT THE EXTENDED BASELINE
# ==============================================================================
ggplot() +
  # Historical backcast (1645-2026)
  geom_line(data = df_backcast, aes(x = dt.mnth, y = baseline_backcast),
            color = "firebrick", size = 1) +
  # Observed data window overlay (1850-2026)
  geom_line(data = Ocean_solar_anom, aes(x = dt.mnth, y = Anomaly),
            alpha = 0.35, color = "gray20") +
  # Highlight Little Ice Age / Maunder Minimum window
  annotate("rect", xmin = 1645, xmax = 1715, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "blue") +
  annotate("text", x = 1680, y = max(df_backcast$baseline_backcast) * 0.8,
           label = "Maunder Minimum", fontface = "italic") +
  labs(
    title = "Astrophysical Harmonic Reconstruction & Backcast",
    subtitle = "Extending long-period baseline (Gleissberg, Suess, Eddy) prior to 1850",
    x = "Year", y = "SST Anomaly (°C)"
  ) +
  theme_minimal()
