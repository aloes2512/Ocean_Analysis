Ocean_solar_anom<-readRDS("data/Ocean_solar_anom")
Ocean_solar_anom$Solar_Retained_Residuals <- Ocean_solar_anom$Anomaly - Ocean_solar_anom$Clean_Baseline_Trend
library(tidyverse)
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Solar_Retained_Residuals),col="grey")+
  geom_line(aes(y=Clean_Baseline_Trend),col=2)+
  labs(title = "Decomposition Ocean Anomalies",
       subtitle = "phaselocked to solarpower\n trend fitted to solar power")

# total harmonics
Ocean_decomp_anomalies<-Ocean_solar_anom%>%
  mutate(Anomaly.total=Solar_Retained_Residuals+har.trd)%>%
  dplyr::select(dt.mnth,Anomaly.total,trd3)
Ocean_decomp_anomalies%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Anomaly.total))+
  geom_smooth(aes(y=Anomaly.total))+
  labs(x="",title="Solar Locked Anomalies")
#-----
anoma.total<-Ocean_decomp_anomalies$Anomaly.total
dt.mnth=Ocean_decomp_anomalies$dt.mnth
last(dt.mnth) # 2026.25 Apr 2026
extend_signal <- function(x, n_back = 3000) {
  N <- length(x)
  N_ext=N+n_back
  t_ext= seq(from=-(n_back-1),N)

  idx_ext=((t_ext - 1) %% N) + 1
  x_ext <- x[idx_ext]
  return(x_ext)

}
sign_ext=extend_signal(anoma.total)
N=length(sign_ext) # 5116
library(itsmr)
Anoma.ext=tibble(
  dt.ext=seq(1600,last(dt.mnth),by=1/12),
  sign_ext=extend_signal(anoma.total)) %>%
  mutate(sign.long.hr=hr(sign_ext,N/1:20))
Anoma.ext%>%
  ggplot(aes(x=dt.ext))+
  geom_line(aes(y=sign_ext),col="grey")+
  geom_line(aes(y=sign.long.hr),col=2)
# residuals only
res.anomalies=Ocean_solar_anom$Solar_Retained_Residuals
Res.anoma.ext=tibble(
  dt.ext=seq(1600,last(dt.mnth),by=1/12),
  ext.anoma=extend_signal(res.anomalies))%>%
    mutate(long.ext.hr=hr(ext.anoma,N/1:20))
Res.anoma.ext%>%ggplot(aes(x=dt.ext))+
  #geom_line(aes(y=ext.anoma),col="grey")+
  geom_line(aes(y=long.ext.hr),col=2)
library(gsignal)
ext.pks=findpeaks(Res.anoma.ext$long.ext.hr,DoubleSided = T)
ext.pks$loc%>%diff()%>% mean() # 11.18 years
sel=seq(10,12, by=0.2)*12

tibble(
  dt.ext=seq(1600,last(dt.mnth),by=1/12),
  ext.anoma=extend_signal(res.anomalies))%>%
  mutate(long.ext.hr=hr(ext.anoma,sel))%>%
  ggplot(aes(x=dt.ext))+geom_line(aes(y=long.ext.hr
                                      ))
par(mar=c(2,2,1,1))
periodogram(sunspots,c(1,1,1,1))
length(sunspots) #2820 month
library(itsmr)

library(itsmr)

# 1. Compute periodogram vector
spec <- periodogram(sunspots, opt = 0)

# 2. Reconstruct the angular frequency grid (omega)
n <- length(sunspots)
k <- 1:length(spec)
freq <- 2 * pi * k / n   # Angular frequency in rad/time-step

# 3. Calculate 95% White Noise Significance Threshold
sigma2 <- var(sunspots)
alpha <- 0.05
M <- length(spec)
p_crit <- -sigma2 * log(1 - (1 - alpha)^(1/M))

# 4. Filter for significant frequencies/powers
sig_idx <- which(spec > p_crit)

freq_sig <- freq[sig_idx]
spec_sig <- spec[sig_idx]

# Convert significant frequencies to periods (e.g., in years/months)
periods_sig <- (2 * pi) / freq_sig

# 5. Plot the spectrum with threshold overlay
plot(freq, spec, type = "l", col = "gray",
     xlab = "Angular Frequency (rad/unit time)",
     ylab = "Periodogram Power",
     main = "Periodogram with 95% Threshold")

# Add threshold line and highlight significant peaks
abline(h = p_crit, col = "red", lty = 2, lwd = 2)
points(freq_sig, spec_sig, col = "blue", pch = 19, cex = 0.8)
1/freq_sig
