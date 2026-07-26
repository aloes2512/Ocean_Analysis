df=readRDS("data/Ocean_solar_anom.rds") #2068 obs
colnames(df)
library(tidyverse)
library(itsmr)
Ocean_anoma<-Ocean_anoma<-df%>%
  dplyr::select(dt.mnth,"Anomaly"=anom,Clean_Baseline_Trend,Solar_Retained_Residuals)
N=NROW(Ocean_anoma) #2116
# eliminate annual and semianual (at equator)

M=c("season",12,"season",6)
Ocean_anomaly=Ocean_anoma%>%  mutate(anom=Resid(Anomaly,M),
                                     trd5=trend(anom,5),
                                     res5=Resid(anom,5))

Ocean_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anom),col="grey")+
  geom_line(aes(y=trd5),col=2,linewidth=1.3)

# combine solar_monthly with ocean anomaly
solar_mnthly=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))

# combine solar and ocean data
Ocean_solar_anom=Ocean_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,anom,"Solar_Variation"=SI,trd5)

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
trend_fit <- lm(anom ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                data = Ocean_solar_anom)


# 3. Extract the stable Baseline Trend and visualise
Ocean_solar_anom$Clean_Baseline_Trend <- predict(trend_fit)
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Clean_Baseline_Trend ),col=2)+
  labs(title="Clean_Baseline_Trend",
       subtitle = "extracted with spline fit solar power degree= 3")
#===========
Ocean_solar_anom=Ocean_solar_anom%>%mutate(har.trd=Clean_Baseline_Trend-trd5)

FFT.har.trd=tibble(idx=1:N-1,
                   spc=fft(Ocean_solar_anom$har.trd),
                   amp=Mod(spc))
FFT.ordrd<-FFT.har.trd%>%arrange(desc(amp))%>%subset(idx<N/2)
prds.ordrd<-FFT.ordrd%>%mutate(prds=N/(idx))%>%pull(prds) # mx.prds yrs: 58.8;22;44.1;25.2;35.3;88.2
Ocean_solar_anom=Ocean_solar_anom%>%
  mutate(max.hr.trd=hr(har.trd,prds.ordrd[1:5])) # periods manually from FFT.har.trd
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=max.hr.trd),col=2,linetype = 2)+
  geom_line(aes(y=Clean_Baseline_Trend))+
  geom_line(aes(y=Clean_Baseline_Trend-max.hr.trd))+
  labs(x="",title = "5 max. harm in sol.lckd trd")
#========

# Compare spline fitted trend with baseline trend

Ocean_solar_anom%>%ggplot(aes(x=dt.mnth,y=Clean_Baseline_Trend))+
  geom_line(col=2)+
  geom_line(aes(y=trd5),col=4)+
  geom_line(aes(y=Clean_Baseline_Trend-trd5),col=9)+
  labs(title="Compare Ocean Anomaly Trends",
       subtitle="5th° poly(blue),Baseline Trend(red),\ndifference (black) ")
# residual will now fully retain the solar peaks, troughs, and historical minima!
Ocean_solar_anom$Solar_Retained_Residuals <- Ocean_solar_anom$anom - Ocean_solar_anom$Clean_Baseline_Trend
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth)) +
  geom_line(aes(y=Clean_Baseline_Trend),col=2)+
  labs(title = "Decomposition Ocean Anomalies",
       subtitle = " trend fitted to solar power")
saveRDS(Ocean_solar_anom,"data/Ocean_solar_anom.rds")
# Long_periods==Ocean_solar_anom+smth.res (Solar_Retained_Residuals,f=0.02)
#====================
Ocean_solar_anom%>%
  mutate(smth.res.11=smooth.fft(anom,f=0.065),
         smth.res.3.4=smooth.fft(anom,f=0.02))%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anom),col="grey")+
  #geom_line(aes(y=smth.res.11),col=2)+
  geom_line(aes(y=smth.res.3.4),col=4)+
  labs(x="",title = "Resids of sol.phaselocked-trend",
       subtitle="low-pass:cutoff 11.2 years (red) 3.4 years (blue)")

#==================
# dominant long periods from phase locked retained resids and phase locked trend
#   diff trd5-Cleaned_Baseline_trend / smth.res = smooth.fft(Solar_Retained...)
Long_periods=Ocean_solar_anom%>%
  mutate(smth.res=smooth.fft(Solar_Retained_Residuals,f=0.02))
# Analyze trends: Clean_Baseline_Trend)-poly trend 3°

Long_periods%>% mutate(trd.period=trd5-Clean_Baseline_Trend,
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
# add trd.harm to anomaly to get long.sum
#apply smth.res, trd.period, long.sum
Long_periods=Long_periods%>% mutate(trd.period=trd5-Clean_Baseline_Trend,
                                    long.sum=trd.period+smth.res)
N =length(Long_periods$dt.mnth) # 2116
n_back=3000

Long_dominant=Long_periods%>%
  dplyr::select(dt.mnth,smth.res,trd.period,long.sum)
Long_periods$dt.mnth[1] # 1850
# extend backwards
my_dominant=Long_dominant$long.sum
Ext.dominant= tibble(dates_ext= seq(1604,by=1/12,length.out=N+n_back),
                     anoma.ext=smooth.fft(extend_backwards(my_dominant,n_back = 3000),f=0.01))

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

