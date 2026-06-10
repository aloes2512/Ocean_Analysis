df=readRDS("data/NOAA.ocean.anomalies.rds") #2068 obs
library(tidyverse)
# median not used here
Ocean_anoma<-df$data%>%dplyr::select(-median.anom)
# eliminate annual and semianual (at equator)
library(itsmr)
M=c("season",12,"season",6)

Ocean_anomaly=Ocean_anoma%>%  mutate(anom=Resid(mean.anom,M),
                                     trd3=trend(anom,3),
                                     res3=Resid(anom,3))
# format date as numeric (year decimals)
library(lubridate)
Ocean_anomaly=Ocean_anomaly%>%
  mutate(Year=year(dt.mnth),
         mnth= (month(dt.mnth)-1)/12,
         dt.mnth=Year+mnth)
Ocean_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anom),col="grey")+
  geom_line(aes(y=trd3),col=2,linewidth=1.3)
# combine solar_monthly with ocean anomaly
solar_mnthly=readRDS("data/solar_power.variation.rds")%>%
  dplyr::select(Time,SI)

solar_mnthly=solar_mnthly%>%rename("dt.mnth"=Time)
# combine solar and ocean data
Ocean_solar_anom=Ocean_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,mean.anom,"Solar_Variation"=SI,trd3)
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
trend_fit <- lm(mean.anom ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                data = Ocean_solar_anom)

# 3. Extract the stable Baseline Trend and visualise
Ocean_solar_anom$Clean_Baseline_Trend <- predict(trend_fit)
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Clean_Baseline_Trend ),col=2)+
  labs(title="Clean_Baseline_Trend",
       subtitle = "extracted with spline fit solar power")
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth,y=Clean_Baseline_Trend))+
  geom_line(col=2)+
  geom_line(aes(y=trd3),col=4)+
  geom_line(aes(y=Clean_Baseline_Trend-trd3),col=9)+
  labs(title="Compare Ocean Anomaly Trends",
       subtitle="3rd° poly(blue),Baseline Trend(red),\ndifference (black) ")
# eliminate season & polynomial trend
M3=c("season",12,"season",6,"trend",3)
Ocean_solar_anom=Ocean_solar_anom%>%mutate(res3=Resid(mean.anom,M3))
# 4. Calculate residuals of solar phase locket ocean anomalies
# residual will now fully retain the solar peaks, troughs, and historical minima!
Ocean_solar_anom$Solar_Retained_Residuals <- Ocean_solar_anom$mean.anom - Ocean_solar_anom$Clean_Baseline_Trend
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Solar_Retained_Residuals),col="grey")+
  geom_line(aes(y=Clean_Baseline_Trend),col=2)+
  labs(title = "Decomposition Ocean Anomalies",
       subtitle = "phaselocked to solarpower\n trend fitted to solar power")
saveRDS(Ocean_solar_anom,"data/Ocean_solar_anom.rds")
#====================
Ocean_solar_anom%>%
  mutate(smth.res=smooth.fft(Solar_Retained_Residuals,f=0.02))%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Solar_Retained_Residuals),col="grey")+
  geom_line(aes(y=smth.res),col=2)+
  labs(x="",title = "Resids of sol.phaselocked-trend",
       subtitle="low-pass filtered (red); cutoff 3.4 years")
#==================
# dominant long periods from phase locked retained resids and phase locked trend
#   diff trd3-Cleaned_Baseline_trend / smth.res = smooth.fft(Solar_Retained...)
Long_periods=Ocean_solar_anom%>%
  mutate(smth.res=smooth.fft(Solar_Retained_Residuals,f=0.02))
colnames(Long_periods)
# Analyze trends: Clean_Baseline_Trend)-poly trend 3°

Long_periods%>% mutate(trd.period=trd3-Clean_Baseline_Trend,
                       long.sum=trd.period+smth.res)%>%
  ggplot(aes(x=dt.mnth))+geom_line(aes(y=long.sum),col=2)+
  geom_line(aes(y=smth.res),col=3,linetype = 2)+
  geom_line(aes(y=trd.period),col=4,linetype = 2)+
  labs(x="",title = "Sum of dominant harmonics",
       subtitle = "harmonic part of phaselocked trend \n smoothed sol. retained res")
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
#apply smth.res, trd.period, long.sum
Long_periods=Long_periods%>% mutate(trd.period=trd3-Clean_Baseline_Trend,
                       long.sum=trd.period+smth.res)
colnames(Long_periods)
N =length(Long_periods$dt.mnth) # 2068
n_back=3000
#
Long_dominant=Long_periods%>%
  dplyr::select(dt.mnth,smth.res,trd.period,long.sum)
Long_periods$dt.mnth[1] # 1854
my_dominant=Long_dominant$long.sum
Ext.dominant= tibble(dates_ext= seq(1604,by=1/12,length.out=N+n_back),
                      anoma.ext=smooth.fft(extend_backwards(my_dominant,n_back = 3000),f=0.01))
colnames(Ext.dominant) # dates_ext; anoma_ext
Ext.dominant%>%
  ggplot(aes(x = dates_ext, y = anoma.ext)) +
  geom_line()+
  # 1. Background Rectangles (Dates must match column type)
  annotate("rect",
           xmin = 1645, xmax = 1715+11/12,
           ymin = -Inf, ymax = Inf, fill = "blue", alpha = 0.15) +
  annotate("rect",
           xmin = 1730, xmax = 1750+11/12,
           ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.15) +

  # 2. Reference Line for the end of LIA
  geom_vline(xintercept = 1850,
             linetype = "dashed", color = "darkred") +
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
  theme_minimal() +
  labs(title = "Ext. sum of Dominant Harmonics",
       subtitle = "harmonic part of phaselocked trend\n smoothed retained resids",
       x = "Year", y = "Dominant Harmonic")
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
           x = 1810,                    # Centered between 1790 and 1830
           y = 0,
           label = "Dalton Minimum",
           angle = 90,
           vjust = 0.5,
           size = 3.5,
           color = "darkgreen",
           fontface = "bold") +

  theme_minimal() +
  labs(title = "Ext. sum of Dominant Harmonics",
       subtitle = "harmonic part of phaselocked trend\n smoothed retained resids",
       x = "Year", y = "Dominant Harmonic")
