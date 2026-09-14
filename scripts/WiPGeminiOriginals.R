library(itsmr)

# 1. Residuals containing absorbed long-period dynamics
resids_solar <- Ocean_solar_anom$Anomaly - Ocean_solar_anom$trd.solar
N <- length(resids_solar) # ~2116 months

# 2. Define exact sequence of periods (from 2116 months down to 105.8 months)
# N / (1:20) gives periods: 2116, 1058, 705.3, ..., 105.8 months
periods <- N / (1:20)

# 3. Pass non-integer harmonic ratios directly to hr()
# Since h = N / period, passing N / periods evaluates to (1:20) exactly,
# or you can pass custom exact solar periods directly in months:
# e.g., custom_periods <- c(2116, 1020, 264, 132) # Gleissberg, Hale, Schwabe
fit_hr <- hr(resids_solar, N / periods)

# 4. Combine reconstructed absorbed periods with your smooth secular trend
Ocean_solar_anom$sol_recstr_hr <- fit_hr
Ocean_solar_anom$trend_with_absorbed_solar <- Ocean_solar_anom$secular_trd + fit_hr
#222222222
# Exact periods in months:
# Schwabe (~11 yrs = 132 mo), Hale (~22 yrs = 264 mo), Gleissberg (~85 yrs = 1020 mo)
solar_periods_months <- c(1020, 540, 264, 132, 105.8)

# Calculate fractional harmonic numbers h = N / T
h_custom <- N / solar_periods_months

# Fit exact physical solar harmonics
fit_hr_exact <- hr(resids_solar, h_custom)
#3333333
#Step-by-Step Rssa Workflow
library(Rssa)

# 1. Prepare the target series (e.g., your hr() reconstructed wave or filtered residuals)
target_series <- Ocean_solar_anom$sol_recstr_hr # or resids_solar

# 2. Choose window length L
# L must be large enough to resolve your longest expected frequency (e.g., at least
# 2-3 times the period of interest, typically 240 to 400 months for decadal/multi-decadal bands)
N <- length(target_series)
L_window <- 360 # 30 years window

# Run SSA decomposition
s_ssa <- ssa(target_series, L = L_window, kind = "toeplitz-ssa")

# 3. Inspect eigenvalues (scree plot) to find dominant variance and coupled pairs
plot(s_ssa, type = "values", main = "SSA Eigenvalues (Screeplot)")

# 4. Check the w-correlation matrix for groups 1 to 20
# Blocks of high off-diagonal w-correlation indicate coupled pairs representing oscillations
plot(s_ssa, type = "wcor", groups = 1:20)

# 5. Periodogram of individual Elementary Components (ECs)
# This directly reveals the frequency of each reconstructed component
plot(s_ssa, type = "values",main = "SSA Eigenvalues (Screeplot)"
# or periodogram via spec.pgram on individual reconstructions
#44444444
#Interpreting the SSA Results for Solar Cycles
# Reconstruct specific groups identified from the w-correlation and eigenvalue plots
# (Replace c(1,2) with the actual paired indices found in your SSA analysis)
rec <- reconstruct(s_ssa, groups = list(Cycle_11yr = c(1, 2), Cycle_multidecadal = c(3, 4)))

# Plot the individual reconstructed components
plot(rec)

# Check the exact frequency/period of the reconstructed components using spec.pgram
spec.pgram(rec$Cycle_11yr, main = "Periodogram of SSA Reconstructed Schwabe Component")
#555555555
library(Rssa)
library(tidyverse)

# 1. Prepare series and run SSA
target_series <- Ocean_solar_anom$sol_recstr_hr # or resids_solar
L_window <- 360 # 30-year window (adjust based on series length)

s_ssa <- ssa(target_series, L = L_window, kind = "toeplitz-ssa")

# 2. Inspect w-correlation to group pairs (optional, adjust group numbers based on your scree plot)
# Reconstruct top component pairs (e.g., pairs 1-2, 3-4, 5-6)
rec <- reconstruct(s_ssa, groups = list(
  Comp_1_2 = c(1, 2),
  Comp_3_4 = c(3, 4),
  Comp_5_6 = c(5, 6)
))

# 3. Build a tidy data frame for ggplot
df_ssa <- Ocean_solar_anom %>%
  select(dt.mnth) %>%
  mutate(
    `Harmonic Input` = target_series,
    `Mode 1 (Pair 1-2)` = as.numeric(rec$Comp_1_2),
    `Mode 2 (Pair 3-4)` = as.numeric(rec$Comp_3_4),
    `Mode 3 (Pair 5-6)` = as.numeric(rec$Comp_5_6)
  ) %>%
  pivot_longer(
    cols = -dt.mnth,
    names_to = "Component",
    values_to = "Amplitude"
  )

# 4. High-resolution ggplot with custom dimensions
ggplot(df_ssa, aes(x = dt.mnth, y = Amplitude, color = Component)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ Component, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 14) + # Larger font scaling for MacBook screens
  labs(
    title = "SSA Decomposition of Absorbed Solar Dynamics",
    x = "Date",
    y = "Anomaly / Amplitude"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 16)
  )
#6666666
# Tips for plotts
ggsave("SSA_Components_Plot.png", width = 12, height = 8, dpi = 300)
# Or export to vector PDF for crisp zooming:
ggsave("SSA_Components_Plot.pdf", width = 12, height = 8)
#777777777777
#Find Eigenvalue Pairs & Periods
library(Rssa)
library(tidyverse)

target_series <- Ocean_solar_anom$sol_recstr_hr
L_window <- 360 # 30-year window

s_ssa <- ssa(target_series, L = L_window, kind = "toeplitz-ssa")

# Extract top 20 elementary components
n_comp <- 20
recs_all <- reconstruct(s_ssa, groups = as.list(1:n_comp))

# Find dominant period (in months) for each component using spec.pgram
comp_summary <- tibble(Component = 1:n_comp) %>%
  mutate(
    Variance_Pct = s_ssa$sigma[1:n_comp]^2 / sum(s_ssa$sigma^2) * 100,
    Dominant_Period_Months = map_dbl(1:n_comp, function(i) {
      s <- as.numeric(recs_all[[i]])
      spec <- spec.pgram(s, plot = FALSE)
      # Convert dominant frequency to period (in months)
      dom_freq <- spec$freq[which.max(spec$spec)]
      return(1 / dom_freq)
    })
  )

# Print the summary table to your R console
print(comp_summary, n = 20)
#Large, Scalable $W$-Correlation Matrix in ggplot2
# Extract W-correlation matrix for the first 20 components
w_mat <- wcor(s_ssa, groups = 1:20)

# Convert matrix to a long data frame for ggplot
df_wcor <- expand.grid(F1 = 1:20, F2 = 1:20) %>%
  mutate(
    Correlation = as.vector(w_mat[1:20, 1:20]),
    F1 = factor(F1),
    F2 = factor(F2, levels = rev(1:20)) # Reverse so 1 is at top
  )

# Plot large heatmap
ggplot(df_wcor, aes(x = F1, y = F2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), size = 3) + # Print numerical correlation values
  scale_fill_gradient(low = "white", high = "firebrick") +
  theme_minimal(base_size = 14) +
  labs(
    title = "SSA W-Correlation Matrix (Groups 1 to 20)",
    x = "Component Number",
    y = "Component Number"
  ) +
  coord_fixed()
#888888888
library(Rssa)
library(tidyverse)

# 1. Name your identified groups based on their periodic bands
solar_groups <- list(
  Gleissberg_Secular1 = c(1, 3),
  Gleissberg_Secular2 = c(2, 4),
  Hale_22yr          = c(5, 6),
  Decadal_SubHale    = c(7, 8),
  Schwabe_11yr       = c(9, 11, 12),
  SubDecadal_5.5yr   = c(13, 14))
)

# 2. Reconstruct components using Rssa
rec_solar <- reconstruct(s_ssa, groups = solar_groups)

# 3. Build a clean, tidy data frame for ggplot
df_rec <- Ocean_solar_anom %>%
  select(dt.mnth) %>%
  mutate(
    `1. Gleissberg/Secular A` = as.numeric(rec_solar$Gleissberg_Secular1),
    `2. Gleissberg/Secular B` = as.numeric(rec_solar$Gleissberg_Secular2),
    `3. Hale (~22 yr)`        = as.numeric(rec_solar$Hale_22yr),
    `4. Decadal Mode`         = as.numeric(rec_solar$Decadal_SubHale),
    `5. Schwabe (~11 yr)`     = as.numeric(rec_solar$Schwabe_11yr),
    `6. Sub-decadal (~5.5 yr)`= as.numeric(rec_solar$SubDecadal_5.5yr)
  ) %>%
  pivot_longer(
    cols = -dt.mnth,
    names_to = "Mode",
    values_to = "Amplitude"
  )

# 4. Large, stacked ggplot with facets
ggplot(df_rec, aes(x = dt.mnth, y = Amplitude, color = Mode)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Mode, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "SSA Decomposed Solar Cycles (Absorbed Signal)",
    x = "Date",
    y = "Temperature Anomaly Contribution (°C)"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 16)
  )
#999999999
#What to Look for in the Facet Plot

#Facet 5 (Schwabe_11yr): You should see clean, quasi-11-year oscillations whose peak heights rise and fall over time (reflecting the modulation carried by component 12).

#Facet 3 (Hale_22yr): Should show smooth ~22-year waves (half the frequency of the Schwabe mode).

#Summing the Reconstruction: If you sum these reconstructed modes together with your original smooth secular_trd, you will have fully recovered the absorbed solar cycles while keeping your zero-crossing spline framework completely intact!

#  What are the exact periods (in months) reported by your summary list for the c(9,11,12) triplet?