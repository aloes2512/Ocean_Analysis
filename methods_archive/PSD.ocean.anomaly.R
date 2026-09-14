library(Rssa)
library(tidyverse)

Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
df<-Ocean_solar_anomaly%>%dplyr::select(dt.mnth,Anomaly,trd.solar)
# 1. Compare Power Spectral Density (PSD) of Trend vs Original Residual
spec_trd <- spec.pgram(df$trd.solar, plot = FALSE)
spec_res <- spec.pgram(df$Anomaly - df$trd.solar, plot = FALSE)

df_psd <- data.frame(
  freq = spec_trd[["freq"]],
  period = 1/as.numeric(spec_trd[["freq"]]),
  PSD_Trend = spec_trd[["spec"]],
  PSD_Residual = spec_res[["spec"]]
)

# Plot spectra to identify overlapping harmonic frequencies
subset(df_psd,period< 100)%>%ggplot(aes(x = period)) +
  geom_line(aes(y = PSD_Trend, color = "Solar-Locked Trend"), size = 0.9) +
  geom_line(aes(y = PSD_Residual, color = "Original Residual"), size = 0.9, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Spectral Distribution: Solar-Locked Trend vs Residual",
    x = "Period (Years, log scale)",
    y = "Power Spectral Density",
    color = "Series"
  ) +
  theme_minimal()

# 2. Joint SSA Decomposition on raw Anomaly to group all distributed harmonics
# Choose L to cover multi-decadal cycles (~50-60 years of annual data)
L_window <- floor(nrow(df) / 3)
s_joint <- ssa(df$Anomaly, L = L_window)

# 3. Compute W-Correlation Matrix to identify harmonic pair groups
w_mat <- wcor(s_joint, groups = 1:20)
plot(w_mat, main = "W-Correlation Matrix (Identify Harmonic Pairs)")
w_df=tibble(w_mat)
library(Rssa)
library(tidyverse)

# 1. Calculate W-correlation matrix (e.g., for components 1 to 20)
w <- wcor(s, groups = 1:20)

# 2. Extract matrix & convert to tidy format
w_mat <- as.matrix(w)

w_df <- w_mat %>%
  tibble() %>%
  rownames_to_column(var = "F1") %>%
  pivot_longer(-F1, names_to = "F2", values_to = "correlation") %>%
  mutate(
    # Clean up column names (removes "F" prefix if present) and set ordering
    F1 = factor(gsub("F", "", F1), levels = 1:ncol(w_mat)),
    F2 = factor(gsub("F", "", F2), levels = ncol(w_mat):1) # Reversed for top-left matrix layout
  )
w_df<-w_df%>%dplyr::select(-F2)
colnames(w_df)
# 3. Plot high-resolution heat map
ggplot(w_df, aes(x = F1, y = F2, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = "grey10", limits = c(0, 1)) +
  coord_fixed() + # Forces square matrix tiles
  labs(
    title = "W-Correlation Matrix",
    x = "Eigenvector Index",
    y = "Eigenvector Index",
    fill = "W-Cor"
  ) +
  theme_minimal(base_size = 14) + # Scaled text size for MacBook Retina displays
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black")
  )

# 4. Extract Scree Plot / Singular Values to check component significance
plot(s_joint, type = "values", main = "Eigenvalue Spectrum")

# 5. Group components based on W-cor inspection:
# Example grouping based on spectral paired structures:
# Group 'trend': Lead singular components (1 or 1-2)
# Group 'harmonics': Identified pairs (e.g., 2:3 for ~60yr, 4:5 for ~22yr, 6:7 for ~11yr)
rec <- reconstruct(s_joint, groups = list(
  secular_trend = 1,
  all_harmonics = c(2, 3, 4, 5, 6, 7) # Adjust according to wcor pairing
))

# 6. Isolated Non-Harmonic Residual (Anomaly - Secular Trend - All Harmonics)
df <- df %>%
  mutate(
    ssa_secular = rec$secular_trend,
    ssa_harmonics = rec$all_harmonics,
    non_harmonic_residual = Anomaly - (ssa_secular + ssa_harmonics)
  )

# Plot isolated non-harmonic residual showing 1910 and 1950 dips cleanly
ggplot(df, aes(x = dt.mnth, y = non_harmonic_residual)) +
  geom_line(color = "#D55E00", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title = "True Non-Harmonic Residual Profile (Stripped of Distributed Harmonics)",
    y = "Anomaly Residual (°C)",
    x = "Year"
  ) +
  theme_minimal()
