library(tidyverse)
library(Rssa)
library(gsignal)
sp_nms <- c("Year", "mth", "dec.dt", "SP_N", "SP_sd", "N.obs", "prov_Fin")

sn_month <- read.table(
  file = file.path(Sys.getenv("HOME"), "Downloads", "SN_m_tot_V2.0.csv"),
  header = FALSE,
  sep = ";",          # The file is using semicolons as delimiters here
  col.names = sp_nms,
  comment.char = "#",
  fill = TRUE,
  stringsAsFactors = FALSE
)
# Clean up any trailing empty columns or semicolons if necessary
sn_month <- sn_month %>%
  dplyr::select(1,2,4)%>%mutate(dt.mnth=Year+(mth-1)/12)
#browseURL("http://www.lajpe.org/dec11/LAJPE_576_Ambelu_Tebabal_prreprint_corr.pdf")
a= 1365.66 # W/m^2
b= 0.0136  # W/(SP*m^2)
c=-0.0000076 # W/(SP*m)^2
#TSI= a+ b*SP+c*SP^2
S_power_ext=sn_month%>%mutate(TSI=a+b*SP_N+c*SP_N^2,TSI=1000*TSI/a)%>%
  dplyr::select(dt.mnth,TSI)%>%
  mutate(S.irr=TSI-mean(TSI,na.rm=T))
S_power_variation<-S_power_ext%>%mutate(SI=TSI-mean(TSI,na.rm=T))

# 1. Load the full SILSO sunspot series from 1749 onward
# (Assuming your source script or file brings in monthly sunspot numbers from 1749)
#source("scripts/Sun_power_full_history.R") # Or your extended loading script

# Ensure your extended solar power data frame spans from 1749 to present
# e.g., S_power_full has columns: dt.mnth, TSI, SI (solar peak power)

# 2. Convert to a time series object covering the full historical span
solar_full_ts <- ts(S_power_variation$SI, start = c(1749, 1), frequency = 12)

# 3. Apply SSA with an appropriately large embedding window
# (e.g., L = 1200 months / 100 years, to fully capture Gleissberg and Hale cycles)
s_ssa_full <- ssa(solar_full_ts, L = 1200)

# 4. Reconstruct with correct harmonic pairing based on the extended spectrum
# Component 1-2: Gleissberg cycle / Long-term envelope (~86-90 yr)
# Component 3-6: Hale cycle (~22 yr)
# Component 7+: Schwabe fundamental (~11 yr)
solar_recon_full <- reconstruct(s_ssa_full, groups = list(
  Trend = c(1, 2),
  Hale       = c(3, 4, 5, 6),
  Schwabe    = c(7:14)
))

# 5. Extract Hilbert analytic envelope for the Gleissberg component over the full timeline
Trend_vec <- as.numeric(solar_recon_full$Trend)
analytic_sig_full <- gsignal::hilbert(Trend_vec, length(Trend_vec))

# 6. Build the extended plotting data frame
plot_df_full <- S_power_variation %>%
  mutate(
    Trend = Trend_vec,
    Hale_Component = as.numeric(solar_recon_full$Hale),
    Schwabe_Component = as.numeric(solar_recon_full$Schwabe),
    Envelope_Upper = Mod(analytic_sig_full),
    Envelope_Lower = -Mod(analytic_sig_full)
  )

# 7. Plot the extended solar power variation from 1749 to present
ggplot(plot_df_full, aes(x = dt.mnth)) +
  geom_line(aes(y = SI, color = "Raw Series (1749+)"), alpha = 0.25, linewidth = 0.3) +
  geom_line(aes(y = Trend, color = "Trend Wave "), linewidth = 0.7) +
  geom_line(aes(y = Envelope_Upper, color = "Trend Envelope"), linetype = "dashed", linewidth = 0.5) +
  geom_line(aes(y = Envelope_Lower, color = "Trend Envelope"), linetype = "dashed", linewidth = 0.5) +
  geom_line(aes(y = Hale_Component, color = "Hale Cycle (~22 yr)"), linewidth = 0.7)+
  scale_color_manual(values =c(
                     "Raw Series (1749+)" = "grey50",
                     "Trend Wave " = "firebrick",
                     "Trend Envelope" = "black",
                     "Hale Cycle (~22 yr)" = "blue3"
                     )) +
  labs(
    y = expression("Solar Power Variation [" * W/m^2 * "]"),
    title = "SSA Harmonic Decomposition of Reconstructed Solar Irradiance",
    subtitle = "Quadratic Sunspot-to-TSI Transformation (1749–Present):\n Schwabe Carrier & Multi-Decadal Modulation",
    color = "Solar Component"
  ) +
  theme_minimal()
