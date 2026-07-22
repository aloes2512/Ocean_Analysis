# Install and load the Rssa package
library(tidyverse)
ssa_url="https://cran.r-project.org/web/packages/Rssa/index.html"
#browseURL(ssa_url)
"Rssa vignette Golyandina"
"asl/rssa"

#browseURL("https://github.com/cran/Rssa/blob/master/man/AustralianWine.Rd")

library(Rssa)
# data  formatting
df=readRDS("data/NOAA.ocean.anomalies.list.rds") #2116 obs
Ocean_anoma<-df$data
library(itsmr)
M=c("season",12,"season",6)
Ocean_anomaly=Ocean_anoma%>%  mutate(anom=Resid(anoma.mean,M),
                                     trd3=trend(anom,3),
                                     res3=Resid(anom,3))
solar_mnthly=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))
Ocean_solar_anom=Ocean_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,"Anomaly"=anoma.mean,"Solar_Variation"=SI,trd3)
# anchor timeseries solar power
library(splines)
# find zero crossings of SI
find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}
nodes <- find_nodes(solar_mnthly$SI)

zero_crossing_times <- solar_mnthly$dt.mnth[find_nodes(solar_mnthly$SI)]

trend_fit <- lm(Anomaly ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                data = Ocean_solar_anom)
Ocean_solar_anom$trend.solar<-predict(trend_fit)
Ocean_solar_anom=Ocean_solar_anom%>%
  mutate(trd5=trend(Anomaly,5),difference=trend.solar-trd5)
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.solar,col="trd.sol"))+
  geom_line(aes(y=difference,col="differnc"))+
  geom_line(aes(y=trd5,col="trd5"))+
  labs(x="",title="Poly5.trend,Sol.lckd.trend")
#difference=trend.solar-trd5
differnc=Ocean_solar_anom$difference
N=length(differnc) # 2116
# 1. Run the SSA decomposition
# L is the window length. A good rule of thumb is L ≈ N/4 to N/2.
# With N = 2116, L = 500 is a robust choice to capture lower-frequency cycles.
s_ocean <- ssa(differnc, L = 500) # maybe 529
str(s_ocean)
# 2. Plot the "Scree Plot" (eigenvalues) to see the variance distribution
# You will likely see a steep drop-off after the first few components.
# Pure oscillations appear as "steps" (pairs of similar eigenvalues).
plot(s_ocean, type = "values")

# 3. Plot the pairing of the eigenvectors (W-correlation matrix)
# This helps you visually group which components belong to the same harmonic.
W <- wcor(s_ocean,groups = 1:15)
class(W) # "wcor.matrix"
dim(W) # 50 50

#alternative plot
# Convert matrix to a long data frame for ggplot
df_wcor <- expand.grid(F1 = 1:15, F2 = 1:15)
df_wcor$Wcor <- as.vector(W)
ggplot(df_wcor, aes(x = F1, y = F2, fill = Wcor)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred") +
  scale_x_continuous(breaks = 1:15, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:15, expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "W-Correlation Matrix",
    x = "Component",
    y = "Component",
    fill = "|W-cor|"
  ) +
  theme_minimal(base_size = 10)
# Estimate frequencies for individual components/pairs
# method = "esprit" is very precise for discrete periodic signals
pe_1.2 <- parestimate(s_ocean, groups = list(Schwabe = c(1, 2)), method = "esprit")
long.period=691.629/12 # 57.63575
pe_34 <-parestimate(s_ocean, groups = list(Schwabe = c(3, 4)), method = "esprit")
pe_34 # 24 years
# 1. Force to standard matrix, convert to table, then to a long data frame.
# This instantly creates 3 clean columns: Var1, Var2, and Freq (the correlation values).
w_tidy <- as.data.frame(as.table(as.matrix(W)))

# 2. Rename the columns and clean up the numeric indices
colnames(w_tidy) <- c("Comp1", "Comp2", "W_Corr")

w_tidy <- w_tidy %>%
  mutate(
    # Extract just the numbers from the component names (e.g., "F1" -> 1)
    Comp1 = parse_number(as.character(Comp1)),
    Comp2 = parse_number(as.character(Comp2)),
    Abs_Corr = abs(W_Corr)
  )

# 3. Create the ggplot
w_plot <- ggplot(w_tidy, aes(x = Comp1, y = Comp2, fill = Abs_Corr)) +
  geom_tile(color = "gray95", linewidth = 0.1) +
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(1, 20, by = 5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(1, 20, by = 5), expand = c(0, 0)) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "W-Correlation Matrix (First 20 Components)",
    subtitle = "Look for dark blue 2x2 squares along the diagonal indicating paired cycles",
    x = "Component Index",
    y = "Component Index",
    fill = "|W-Corr|"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    # This rotates the horizontal axis text by 45 degrees and aligns them perfectly:
    axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

# Display the plot in RStudio
print(w_plot)

# Save a high-resolution, zoomable version for your MacBook screen
ggsave("figs/wcor_matrix_large.png", plot = w_plot, width = 11, height = 10, dpi = 300)
# 4. Reconstruct the clean, periodic signal using only the dominant components
# For example, if the first 4 components capture your harmonics:
reconstructed <- reconstruct(s_ocean, groups = list(Harmonics = 1:20))

# Plot the original difference against your clean SSA reconstruction
plot(differnc, type = "l", col = "red", main = "SSA Reconstruction of Harmonics")
lines(reconstructed$Harmonics, col = "blue", lwd = 2)
#======
# Estimate the exact, continuous frequencies and periods using ESPRIT
pe <- parestimate(s_ocean,
                  groups = list(Cycle1 = c(1, 2),
                                Cycle2 = c(3, 4),
                                Cycle3 = c(5, 6),
                                Cycle4 = c(7, 8),
                                Cycle5 = c(9, 10)),
                  method = "esprit")

# Print the results
print(pe)
307/12


# Check the dominant frequency/period of the first 10 components
# Corrected loop using base R's built-in spectral tools
for (i in 1:10) {
  # 1. Reconstruct the single component (returns a list)
  rec_single <- reconstruct(s_ocean, groups = list(v = i))

  # 2. Extract the actual numeric vector from the reconstructed list
  # Rssa reconstruct returns a list of data frames; we grab the column corresponding to 'v'
  series_data <- rec_single$v

  # 3. Compute the periodogram using base R (fast, robust, and built-in)
  # fast = FALSE prevents padding that can shift the peak slightly
  spec_obj <- spec.pgram(series_data, plot = FALSE, fast = FALSE)

  # 4. Find the frequency with the highest spectral power
  dom_freq <- spec_obj$freq[which.max(spec_obj$spec)]

  # 5. Calculate periods (Frequency is cycles per observation, which is monthly)
  period_months <- 1 / dom_freq
  period_years <- period_months / 12

  cat(sprintf("Component %d: Period = %.2f months (%.2f years)\n",
              i, period_months, period_years))
}

