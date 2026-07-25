df=readRDS("data/NOAA.ocean.anomalies.list.rds") #2116 obs
library(tidyverse)
# median not used here
Ocean_anoma<-df$data%>%dplyr::select(dt.mnth,"Anomaly"=anoma.mean)
N=NROW(Ocean_anoma) #2116
# eliminate annual and semiannual (at equator) and linear trend
library(itsmr)
M=c("season",12,"season",6)
Ocean_anomaly=Ocean_anoma%>%  mutate(anom=Resid(Anomaly,M),
                                     trd3=trend(anom,3),
                                     res3=Resid(anom,3))
library(lubridate)

Ocean_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anom),col="grey")+
  geom_line(aes(y=trd3),col=2,linewidth=1.3)


# combine solar_monthly with ocean anomaly
solar_mnthly=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))

# combine solar and ocean data
Ocean_solar_anom=Ocean_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,Anomaly,"Solar_Variation"=SI,trd3)
# anchor timeseries solar power
library(splines)
# find zero crossings of SI
find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}
nodes <- find_nodes(solar_mnthly$SI)
# 1. Identify the time coordinates of your zero-crossings
zero_crossing_times <- solar_mnthly$dt.mnth[find_nodes(solar_mnthly$SI)]

# 2. Fit a continuous spline locked at these exact nodes
# degree = 1 creates  sawtooth a continuous line that bends at the zero-crossings
# spline degree = 3 creates a smooth, continuous curve through the zero-crossings
trend_fit <- lm(Anomaly ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                data = Ocean_solar_anom)

# 3. Extract the stable Baseline Trend and visualise
Ocean_solar_anom$Clean_Baseline_Trend <- predict(trend_fit)
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Clean_Baseline_Trend ),col=2)+
  labs(title="Clean_Baseline_Trend",
       subtitle = "extracted with spline fit solar power degree= 3")
#===========
Ocean_solar_anom=Ocean_solar_anom%>%mutate(har.trd=Clean_Baseline_Trend-trd3)
saveRDS(Ocean_solar_anom,"data/Ocean_solar_anom")
FFT.har.trd=tibble(idx=1:N-1,
                   spc=fft(Ocean_solar_anom$har.trd),
                   amp=Mod(spc))
#max periods N/2:9
library(itsmr)
Ocean_solar_anom=Ocean_solar_anom%>%
  mutate(max.hr.trd=hr(har.trd,N/2:9)) # periods manually from FFT.har.trd
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=har.trd,col="har.trd"))+
  geom_line(aes(y=Clean_Baseline_Trend,col="Baseline"))+
  geom_line(aes(y=trd3,col="trd3"))
#========
# Compare spline fitted trend with baseline trend

# eliminate season & polynomial trend
#?M3=c("season",12,"season",6,"trend",3)
#?Ocean_solar_anom=Ocean_solar_anom%>%mutate(res3=Resid(anoma.mean,M3))
# 4. Calculate residuals of solar phase locket ocean anomalies
# residual will now fully retain the solar peaks, troughs, and historical minima!
Ocean_solar_anom$Solar_Retained_Residuals <- Ocean_solar_anom$Anomaly - Ocean_solar_anom$Clean_Baseline_Trend
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Solar_Retained_Residuals),col="grey")+
  geom_line(aes(y=har.trd,col="har.trd"))+
  geom_line(aes(y=Clean_Baseline_Trend,col="Baseline"))+
  geom_line(aes(y=trd3,col="trd3"))+
  labs(x="",title = "Decomposed Solar_locked Anomalies",
       subtitle = "Trends: Baseline, polytrd3\n harmonic part of Baseline")
 # total harmonics
Ocean_decomp_anomalies<-Ocean_solar_anom%>%
  mutate(Anomaly.total=Solar_Retained_Residuals+har.trd)%>%
  dplyr::select(dt.mnth,Anomaly.total,trd3)
Ocean_decomp_anomalies%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Anomaly.total))+
  geom_smooth(aes(y=Anomaly.total))+
  labs(x="",title="Solar Locked Anomalies")
#=========================
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Solar_Retained_Residuals),col="grey")+
  geom_line(aes(y=Clean_Baseline_Trend),col=2)+
  labs(title = "Decomposition Ocean Anomalies",
       subtitle = "phaselocked to solarpower\n trend fitted to solar power")
saveRDS(Ocean_solar_anom,"data/Ocean_solar_anom.rds")
# Long_periods==Ocean_solar_anom+smth.res (Solar_Retained_Residuals,f=0.02)
#====================
Ocean_decomp_anomalies<-Ocean_decomp_anomalies%>%
  mutate(hr15=hr(Anomaly.total,N/1:15), # 11.75 years
         hr50=hr(Anomaly.total,N/1:50))

Ocean_decomp_anomalies%>%  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Anomaly.total),col="grey")+
  geom_line(aes(y=hr15),col=2)+
  geom_line(aes(y=hr50),col=4)+
  labs(x="",title = "Resids of sol.phaselocked-trend",
       subtitle="low-pass filtered hr min pperiods > 11.75 yr;> 3.5 yr")

#==================
# dominant long periods from phase locked retained resids and phase locked trend
# Analyze trends: Clean_Baseline_Trend)-poly trend 3°

# ext backwards
library(gsignal)
extend_backwards <- function(x, n_back = 3000) {
  N <- length(x)
  N_ext=N+n_back
  t_ext= seq(from=-(n_back-1),N)
  fft_orig <- fft(x) /N
  idx_ext=((t_ext - 1) %% N) + 1
  fft_ext <- fft_orig[idx_ext]
  Re(ifft(fft_ext) * N_ext)
}
# har.trd contains harmonic periods


Ocean_solar_anom=Ocean_solar_anom%>%
  mutate(har.trdres=trd3-Clean_Baseline_Trend,
         sol.res=Anomaly-Clean_Baseline_Trend,
         anomalies=har.trdres+sol.res)
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anomalies),col=2)+
  geom_line(aes(y=sol.res),col="grey")

N =length(Ocean_solar_anom$dt.mnth) # 2116
n_back=3000


Solar_total_residuals=Ocean_solar_anom%>%
  dplyr::select(dt.mnth,anomalies,sol.res)
Solar_total_residuals$dt.mnth[1] # 1850
# extend backwards
my_dominant=Solar_total_residuals$anomalies
Ext.dominant= tibble(dates_ext= seq(1600,by=1/12,length.out=N+n_back),
                      anoma.ext=hr(extend_backwards(my_dominant,n_back = 3000),N/1:20)
)

#============
Ext.dominant %>%
  ggplot(aes(x = dates_ext, y = anoma.ext)) +
  geom_line() +

  # 1. Background Rectangles (Dates must match column type)
  annotate("rect",
           xmin = 1645, xmax = 1715+11/12,
           ymin = -Inf, ymax = Inf, fill = "blue", alpha = 0.15) +
  annotate("rect",
           xmin = 1730, xmax = 1750+11/12,
           ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.15) +
  # NEW: Added Dalton Minimum Shading
  annotate("rect",
           xmin = 1790, xmax = 1830+11/12,
           ymin = -Inf, ymax = Inf, fill = "darkgreen", alpha = 0.12) +

  # 2. Reference Line for the end of LIA
  geom_vline(xintercept = 1850,
             linetype = "dashed", color = "darkred") +

  # 3. Text Labels
  annotate("text",
           x = 1680,
           y = 0,                       # Set to the middle of your y-axis scale
           label = "Maunder Minimum",
           angle = 90,
           vjust = 0.5,                 # Centers the text on the 'y' coordinate
           size = 3.5,
           color = "blue4",
           fontface = "bold") +

  annotate("text",
           x = 1740,
           y = 0,
           label = "18th C Warmth",
           angle = 90,
           vjust = 0.5,
           size = 3.5,
           color = "orange4",
           fontface = "bold") +

  # NEW: Added Dalton Minimum Text
  annotate("text",
           x = 1810,    # Centered between 1790 and 1830
           y = 0,
           label = "Dalton Minimum",
           angle = 90,
           vjust = 0.5,
           size = 3.5,
           color = "darkgreen",
           fontface = "bold") +

  theme_minimal() +
  labs(title = "Ext.backwards of harmonic Residuals",
       subtitle = "harmonic part of phaselocked trend\n smoothed retained resids",
       x = "Year", y = "Dominant Harmonic")
#=======================

