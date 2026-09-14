## =============================================================================
##  Long-period "astrophysical harmonic" fit and backcast   -- SET ASIDE
##  Status: kept as a starting point, NOT a working method.  Do not revive it
##          on the instrumental record alone; see the conditions at the end.
##  Reviewed 2026-09-14 against ERSSTv6 (1850.00-2026.33, N = 2117)
## -----------------------------------------------------------------------------
##  What it does: regresses the SST anomaly, and separately the solar-locked
##  residual, on sine/cosine pairs at named solar periods -- Hale 20 yr,
##  Gleissberg 88 yr, Suess 208 yr, Eddy 1000 yr -- and extrapolates the fit
##  back to 1600 to compare with the Maunder Minimum.
##
##  Why it was set aside.  Three checks, each independently disqualifying.
##
##  (1) There is no power at those periods to fit.  Tested against an AR(1)
##      red-noise null (the residual has lag-1 autocorrelation 0.77, so the
##      rise towards long periods in the log-log periodogram is red noise, not
##      peaks), the significant bins between 3 and 200 yr are:
##          19.6 yr  x5.95 ] 
##          14.7 yr  x4.86 ]  all > 99 %
##           9.3 yr  x7.44 ]
##        3.2-6.1 yr x3.4-4.6  > 95 %
##      10 of 58 bins clear 95 % where 2.9 are expected by chance, so the
##      structure is real -- but it is all at 20 yr and below.  At 88 yr and at
##      176/208 yr the power is ~0.00 times the red-noise expectation.  Hale and
##      Schwabe are present; Gleissberg, Suess and Eddy are not.
##
##  (2) Over a 176-yr window these are not harmonics.  Regressing each basis
##      vector on {1, t, t^2}:
##          sin_eddy, cos_eddy   R^2 = 0.9998
##          cos_suess            R^2 = 0.926      sin_suess  R^2 = 0.847
##          sin_gleiss           R^2 = 0.151      cos_hale   R^2 = 0.004
##      A 1000-yr sinusoid sampled over 176 yr IS a quadratic.  The design
##      matrix is singular (condition number infinite), so fit_astro_baseline
##      fits a parabola and labels it Eddy + Suess; the backcast is that
##      parabola extrapolated 250 yr.  Hence the minimum at 1731 -- outside the
##      Maunder window -- and the monotone ramp before it.  The p-values from
##      summary() are meaningless under this collinearity.
##
##  (3) The backcast is not stable to the analysis window.  Same model, same
##      four periods, different fitting windows; mean of the backcast over
##      1645-1715:
##          1850-2026   -0.61 K      1870-2026   -1.33 K
##          1850-2010   -2.99 K      1890-2026  +24.25 K
##          1850-2000   +3.26 K
##          1850-1990  +24.95 K
##      Spread 28 K.  The plausible-looking -0.5 K of the full-record fit is an
##      accident of the record ending in 2026.  This is exactly the
##      record-length dependence the solar-locked construction was adopted to
##      avoid (see Draft_Report_files/Paper_SolarLocked_Decomposition.qmd,
##      sections 3.1 and sec-window).
##
##  Searching for a combination of long periods that reproduces the Maunder
##  Minimum cannot fail, and that is the objection: the target is known in
##  advance, the basis functions are indistinguishable from polynomials over
##  the window, and the region being fitted contains no data.  A hit would
##  carry no information about the ocean.
##
##  (4) The long-wavelength term is a usable BASIS but not an identified
##      COMPONENT.  Using a period far longer than the record to carry the
##      secular baseline is a legitimate parametrisation -- but its amplitude
##      is then not a measurement of anything, because it is collinear with
##      any other trend term.  Refitting the anomaly with an explicit
##      industrial term (the Richards curve of scripts/Industrial_growthcurve.R)
##      beside the same four harmonics:
##
##                        harmonics only    + Richards term
##          R^2              0.8539            0.8554
##          Hale             0.027 K           0.026 K
##          Gleissberg       0.065 K           0.090 K
##          Suess            0.201 K           0.671 K
##          Eddy             0.640 K           9.074 K
##          backcast
##          1645-1715       -0.61  K          +5.18  K
##
##      The Eddy amplitude grows fourteenfold and the fit improves by 0.0015
##      in R^2: the trend is being split between two near-duplicate regressors
##      and the split is arbitrary.  Hale is untouched, being genuinely
##      orthogonal to the trend at 20 yr.
##
##      This is what specifically invalidates the BACKCAST, as distinct from
##      the fit.  A backcast interprets the long term as a component -- "this
##      much of the rise was a slow natural cycle, therefore before 1850 it was
##      down there".  The regression can fit the rise but cannot attribute it,
##      so the 17th-century value is set by an arbitrary split.
##
##      DESIGN RULE for any future version:  fit the long-period terms and the
##      industrial term JOINTLY, on a record long enough that the long period
##      is realised several times, and check that the long-period amplitude is
##      stable when the industrial term is added or removed.  If it moves by a
##      factor of ten, as above, the decomposition is not identified however
##      good the fit looks.  Run this test before the window-stability test of
##      (3); it is cheaper and it fails faster.
##
##  What would make this script legitimate.  Fit it on data that CONTAINS the
##  period in question, and use the instrumental era as the test rather than
##  the training set:
##    - a proxy temperature reconstruction covering 1600 onwards (PAGES 2k,
##      or an equivalent ocean-only compilation);
##    - and/or the telescopic sunspot record from 1610, which spans the
##      Maunder Minimum directly and gives ~2.4 Gleissberg cycles instead of
##      1.9 -- still marginal, but the periods would at least be estimated
##      from a window that contains them;
##    - fit over the proxy era, then validate by predicting 1850-2026 and
##      comparing with ERSSTv6.  Reversing that order is what fails above.
##  Under that design the Little Ice Age becomes a test of the model instead
##  of a target for it, and the window-stability check in (3) should be run
##  again before any result is believed.
##
##  Note: scripts/NoLongPeriods_solckd.res.R is a near-duplicate of this file
##  (54 diff lines, mostly the residual column name).  Its commit message
##  reports a leakage ratio of 4.81e-9 and concludes the solar filter isolates
##  the long-period energy completely; a recomputation on ERSSTv6 gave ~5 % for
##  the same quantity, so that figure should be re-derived before it is cited.
## =============================================================================

library(dplyr)
library(ggplot2)
library(broom)

# ==============================================================================
# 0. SETUP & DATA PREPARATION
# ==============================================================================
# Assuming 'Ocean_solar_anom' contains:
#   - dt.mnth: monthly time index (e.g., 1850.0, 1850.083, ...)
#   - Anomaly: global ocean SST anomaly
#   - Solar_Retained_Residuals: residuals from solar-phase locked model
Ocean_solar_anom<-readRDS("data/Ocean_solar_anom.rds")
# Ensure time step in months for harmonic frequency calculations
dt_start <- min(Ocean_solar_anom$dt.mnth)
Ocean_solar_anom <- Ocean_solar_anom %>%
  mutate(time_months = (dt.mnth - dt_start) * 12)

# ==============================================================================
# 1. SPECTRAL DIAGNOSTICS OF RESIDUALS
# ==============================================================================
# A. Standard Periodogram (Check if low frequencies slope upward or stay flat)
spec_pgram <- spec.pgram(Ocean_solar_anom$Solar_Retained_Residuals , plot = FALSE, log = "yes")
df_pgram <- data.frame(
  period_years = (1 / spec_pgram$freq) / 12,
  spec = spec_pgram$spec
)

ggplot(df_pgram %>% filter(period_years <= 176), aes(x = period_years, y = spec)) +
  geom_point(color = "steelblue",shape=1,size= 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Residual Periodogram (Log-Log)",
       x = "Period (Years)", y = "Spectral Density") +
  theme_minimal()

# B. Maximum Entropy / AR Spectrum (Superior resolution for short records)
spec_ar_fit <- spec.ar(Ocean_solar_anom$Solar_Retained_Residuals , plot = FALSE, log = "yes")
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
T_hale       <- 20 *12      # ~ 20 years
T_gleissberg <- 88 * 12     # ~88 years
T_suess      <- 208 * 12    # ~208 years
T_eddy       <- 1000 * 12   # ~1000 years (Millennial scale)

# Construct harmonic predictors
Ocean_solar_anom <- Ocean_solar_anom %>%
  mutate(sin_hale = sin(2 * pi * time_months / T_hale),
         cos_hale = cos(2 * pi * time_months / T_hale),


    sin_gleiss = sin(2 * pi * time_months / T_gleissberg),
    cos_gleiss = cos(2 * pi * time_months / T_gleissberg),
    sin_suess  = sin(2 * pi * time_months / T_suess),
    cos_suess  = cos(2 * pi * time_months / T_suess),
    sin_eddy   = sin(2 * pi * time_months / T_eddy),
    cos_eddy   = cos(2 * pi * time_months / T_eddy)
  )

# Fit harmonic model to the residual series or total low-frequency baseline
fit_astro_resids <- lm(
  Solar_Retained_Residuals  ~ sin_gleiss + cos_gleiss +
    sin_suess  + cos_suess  +
    sin_eddy   + cos_eddy+ sin_hale+cos_hale,
  data = Ocean_solar_anom
)

# View statistical significance of candidate long periods
summary(fit_astro_resids)

# Calculate variance explained by long periods in the residuals
var_total_resids <- var(Ocean_solar_anom$Solar_Retained_Residuals )
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
    sin_eddy   + cos_eddy+sin_hale+cos_hale,
  data = Ocean_solar_anom
)

# Create an extended timeline back to 1645
df_backcast <- data.frame(
  dt.mnth = seq(from = 1600, to = 2026, by = 1/12)
) %>%
  mutate(
    time_months = (dt.mnth - dt_start) * 12,
    sin_hale    = sin(2*pi*time_months / T_hale),
    cos_hale    = cos(2*pi*time_months/T_hale),
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
