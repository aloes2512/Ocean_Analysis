library(tidyverse)
library(Rssa)
Ocean_solar_anomaly=readRDS("data/Ocean_solar_anomaly.rds")
lcked.trend<-Ocean_solar_anomaly$trend.solar
time<-Ocean_solar_anomaly$dt.mnth
# 1. Fit LOESS to extract non-harmonic solar trend
fit_loess <- loess(lcked.trend ~ time, span = 0.25, degree = 2)
my_dat<-Ocean_solar_anomaly%>%dplyr::select(dt.mnth,trend.solar)
my_dat<-my_dat%>%mutate(trend_non_harmonic = predict(fit_loess))
my_dat%>%ggplot(aes(x=dt.mnth))+geom_line(aes(y=trend_non_harmonic))
# 2. Extract stationary residuals for harmonic analysis
my_dat<-my_dat%>%mutate(trd.har=trend.solar-trend_non_harmonic)
library(itsmr)
N<-NROW(Ocean_solar_anomaly)
Ocean_solar_anomaly<-Ocean_solar_anomaly%>%mutate(anoma.res=Anomaly-trend.solar,
                                                  Resd.filt=hr(anoma.res,N/1:20))
wind.1950=time[which.max(Ocean_solar_anomaly$Resd.filt)] #1942.25
resd1950=subset(Ocean_solar_anomaly,dt.mnth>1950)%>%pull(Resd.filt)
end.1950<-1950+which.max(resd1950)/12
window.1950=c(wind.1950,end.1950)
Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+geom_line(aes(y=Resd.filt))
my_dat%>%ggplot((aes(x=dt.mnth)))+
  geom_line(aes(y=trd.har),linetype = 3)+
  geom_line(data=Ocean_solar_anomaly,aes(y=Resd.filt),linetype = 2)+
  geom_vline(xintercept = wind.1950,col=2)+
  geom_vline(xintercept = end.1950,col=2)+
  labs(title="Harmonic Part of trend.solar")
# add the two harmonics
my_dat$Resd.filt<-Ocean_solar_anomaly$Resd.filt
my_dat%>%mutate(tot.Resd=Resd.filt+trd.har)%>%
  ggplot(aes(x=dt.mnth))+geom_line(aes(y=tot.Resd))+
  geom_vline(xintercept = wind.1950,col=2)+
  geom_vline(xintercept = end.1950,col=2)+
  geom_vline(xintercept = 1999,col=3)+
  geom_vline(xintercept = 2004.917,col=3)+
  labs(x="",title = "Long Periods Solar Ocean Anomaly",
       subtitle = "Low-Pass periods > 9 yrs")
library(gsignal)
resds=my_dat%>%subset(dt.mnth>1970)%>%pull(Resd.filt)
pk.lc=1970+findpeaks(resds,DoubleSided = T)$loc /12
findpeaks(resds,DoubleSided = T)$height
pk.lc[6:9]
y_detrended<-my_dat$trd.har
# 3. Clean SSA on residuals (no trend-interference)
s <- ssa(y_detrended, L = 1058) # adjust L to your primary harmonic window
plot(s, type = "paired")        # look for clean circular harmonic pairs

# 4. Reconstruct pure harmonics
# Reconstruct isolated physical solar modes
rec <- reconstruct(s, groups = list(
  Schwabe_11yr = c(1, 2),        # Clean circular pair
  Hale_Overtone = c(3, 4),       # Lissajous harmonic pair
  Amplitude_Mod = c(5, 6)        # Concentric spiral pair
))

# Plot the reconstructed physical components
plot(rec)
library(dplyr)
library(tidyr)
library(ggplot2)

# Convert rec list components to a data frame with date axis
# Inspect the structure and names of the list
names(rec)
str(rec)

# Extract individual components directly as numeric vectors
comp1 <- as.numeric(rec$Schwabe_11yr)
comp2 <- as.numeric(rec$Hale_Overtone)
comp3 <- as.numeric(rec$Amplitude_Mod)
df_rec <- tibble(dt.mnth=my_dat$dt.mnth,
                 Schwabe=comp1,
                 Hale_Overtone=comp2,
                 ampl.mod=comp3) %>%
  pivot_longer(
    cols = -dt.mnth,
    names_to = "Component",
    values_to = "Value"
  )

# Plot overlaid series
ggplot(df_rec, aes(x = dt.mnth, y = Value, color = Component)) +
  geom_line(linewidth = 0.5) +
  theme_minimal() +
  labs(
    title = "SSA Reconstructed Components",
    x = "Date",
    y = "Amplitude",
    color = "Component"
  )
# plot components separated
df_rec <- tibble(dt.mnth=my_dat$dt.mnth,
                 Schwabe=comp1,
                 Hale_Overtone=comp2,
                 ampl.mod=comp3)
df_rec%>%  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Hale_Overtone,col="hale"))+
  geom_line(data=my_dat,aes(y=trd.har,col="har.trd"))

library(seewave) # or practical implementation via FFT

# Calculate the instantaneous analytic amplitude envelope
schwabe_vec <- as.numeric(rec$Schwabe)
env <- abs(seewave::rms(schwabe_vec)) # or hilbert envelope

# Alternatively via basic R built-ins:
h <- fft(schwabe_vec)
h[2:(length(h)/2)] <- 2 * h[2:(length(h)/2)]
h[(length(h)/2 + 1):length(h)] <- 0
analytic_signal <- fft(h, inverse = TRUE) / length(h)
my_dat$schwabe_amplitude <- Mod(analytic_signal)

# Plotting the continuous amplitude envelope curve
ggplot(my_dat, aes(x = dt.mnth)) +
  geom_line(aes(y = trd.har), color = "grey70") +
  geom_line(aes(y = as.numeric(rec$Schwabe)), color = "dodgerblue", linewidth = 0.8) +
  geom_line(aes(y = schwabe_amplitude), color = "firebrick", linetype = "dashed", linewidth = 1) +
  theme_minimal() +
  labs(title = "Schwabe Cycle",subtitle="with Instantaneous Amplitude Envelope", x = "", y = "Amplitude")

#----------
# Calculate average peak-to-peak distance in years
get_peak_distance <- function(vec) {
  v <- as.numeric(vec)
  # Find indices where slope turns from positive to negative
  peaks <- which(diff(sign(diff(v))) == -2) + 1
  mean_diff_months <- mean(diff(peaks))
  return(mean_diff_months / 12) # convert to years
}

sapply(rec, get_peak_distance)
#------
# Helper to test peak distance for any eigenvector pair
check_pair_period <- function(ssa_obj, pair_idx) {
  # Reconstruct just that specific pair
  rec_pair <- reconstruct(ssa_obj, groups = list(p = pair_idx))
  v <- as.numeric(rec_pair$p)
  peaks <- which(diff(sign(diff(v))) == -2) + 1
  period_yrs <- mean(diff(peaks)) / 12
  return(period_yrs)
}

# Scan pairs 1-2, 3-4, 5-6, 7-8, 9-10, 11-12
sapply(list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
       function(pair) check_pair_period(s, pair))
#===========
library(ggplot2)

# 1. Extract raw numeric matrix (first 20 components)
w_obj <- wcor(s, groups = 1:20)
w_mat <- as.matrix(w_obj)

# 2. Build grid directly from matrix dimensions
n_comp <- nrow(w_mat) #20
df_wcor <- expand.grid(Comp1 = 1:n_comp, Comp2 = 1:n_comp)
df_wcor$W_cor <- as.vector(w_mat)

# 3. Plot clean matrix with legible axis labels
ggplot(df_wcor, aes(x = factor(Comp1), y = factor(Comp2), fill = W_cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1)) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "W-Correlation Matrix (Components 1-20)",
    x = "Eigenvector Index",
    y = "Eigenvector Index",
    fill = "W-Cor"
  ) +
  theme(
    axis.text = element_text(size = 9, face = "bold"),
    panel.grid = element_blank()
  )
# tripplet reconstruction because of variable length of intervalls
# Reconstruct the true quasi-periodic Schwabe mode
rec_schwabe <- reconstruct(s, groups = list(Schwabe_Quasi = c(9, 10, 11)))
sapply(rec_schwabe, get_peak_distance)# 15.12
rec_schwabe2 <- reconstruct(s,groups = list(Unknown=c(12,13,14)))
sapply(rec_schwabe2, get_peak_distance) # 7.54
# from the cor plot take group 7,8,9,10
rec_Schwabe_full<- reconstruct(s,groups = list(full=c(7:10)))
sapply(rec_Schwabe_full, get_peak_distance)# 15.2
#==========
# Helper function to get period in years directly from an eigenvector index
get_vec_period <- function(ssa_obj, idx) {
  # Extract eigenvector column
  v <- ssa_obj$U[, idx]
  # Compute periodogram
  s <- stats::spectrum(v, plot = FALSE)
  top_freq <- s$freq[which.max(s$spec)] # cycles per month
  period_months <- 1 / top_freq
  return(period_months / 12)            # convert to years
}

# Scan individual eigenvectors directly
sapply(1:15, function(i) get_vec_period(s, i))
# 22.5 22.5 18.0 18.0 30.0 30.0 45.0 45.0 15.0 15.0  7.5
#  7.5  7.5 90.0 90.0
# Scan eigenvectors 16 through 35 directly
schwabe_scan <- sapply(16:35, function(i) get_vec_period(s, i))
names(schwabe_scan) <- paste0("F", 16:35)
print(schwabe_scan)
# searching the Schwabe in the residuals
# 1. Compute residuals of the solar-locked trend
res_solar <- Ocean_solar_anomaly$Anomaly - my_dat$trend.solar

# 2. Fit SSA directly on the residual series
# (Using L = 264 ~ 22 years is ideal for residuals to focus on Schwabe/Hale scales)
s_res <- ssa(res_solar, L = 264)

# 3. Direct spectral scan of the top residual eigenvectors
res_periods <- sapply(1:10, function(i) get_vec_period(s_res, i))
names(res_periods) <- paste0("F", 1:10)
print(res_periods)
