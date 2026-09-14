df=readRDS("data/NOAA.ocean.anomalies.list.rds") #2116 obs
library(tidyverse)
# median not used here
names(df)#
Ocean_anoma<-df$data
N=NROW(Ocean_anoma) #2116
# eliminate annual and semianual (at equator)
library(itsmr)
M=c("season",12,"season",6)
Ocean_anomaly=Ocean_anoma%>%  mutate(anom=Resid(anoma.mean,M),
                                     trd7=trend(anom,7),
                                     res7=Resid(anom,7))
plt.sol.trd7=Ocean_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_point(aes(y=res7),col="grey",size=0.3)+
  geom_line(aes(y=trd7,col="trd7"))+
  labs(x="",y="anomaly",title = "Global Ocean Anomaly",
       subtitle="trd7 =polynomial degree 7 fitted ")
# combine solar_monthly with ocean anomaly
solar_mnthly=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))

# combine solar and ocean data
Ocean_solar_anom=Ocean_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,"Anomaly"=anoma.mean,"Solar_Variation"=SI,trd7)
# anchor  at time-series solar power
library(splines)
# find zero crossings of Solar Irradiation SI
find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}
# nodes <- find_nodes(solar_mnthly$SI)
# time coordinates of zero-crossings
zero_crossing_times <- solar_mnthly$dt.mnth[find_nodes(solar_mnthly$SI)]

# 2. Fit a continuous spline locked @ exact nodes
# degree = 1 creates  sawtooth a continuous line that bends at the zero-crossings
# spline: degree = 3 creates a smooth, continuous curve through the zero-crossings
trd.solar <- lm(Anomaly ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                data = Ocean_solar_anom)
Ocean_solar_anom$trend.solar<-predict(trd.solar)
Ocean_solar_anom=Ocean_solar_anom%>%
  mutate(trd7=trend(Anomaly,7),difference=trend.solar-trd7)
plt_trd7.trend.solar<-Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.solar,col="trend.solar"))+
  geom_line(aes(y=trd7,col="poly7.trnd"))+
  geom_line(aes(y=difference,col="difference"))+
  labs(x="",title="Compare Solar-locked-trend\nPolynomial Trend 7°")

#ggsave("figs/sol.poly7.trds.tiff")
saveRDS(Ocean_solar_anom,"data/Ocean_solar_anom.rds")
#--------
# 3. Extract the stable Baseline Trend and visualise
plt.decomp_solar.lckd<-Ocean_solar_anom%>%  mutate(solar.res=Anomaly-trend.solar)%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=solar.res),col="grey")+
  geom_line(aes(y=trend.solar,col="trend.solar"))+
  labs(title="Solar-locked Trend",
       subtitle = "extracted with spline fit to solar power")
#===========
# periods of trend.solar estimate
Ocean_solar_anom=Ocean_solar_anom%>%
 mutate(har.res=trend.solar-trd7)
FFT.har.trd=tibble(idx=1:N-1,
                  spc=fft(Ocean_solar_anom$har.res),
                   amp=Mod(spc))
idx.mx=FFT.har.trd%>%arrange(desc(amp))%>%subset(idx< 100)%>%head(20)%>%pull(idx)
FFT.har.trd%>%subset(idx<30)%>%
  ggplot(aes(x=idx))+
  geom_point(aes(y=amp))+
  labs(x="rel.frequency",
       title="Fourier Coef. Harmonic Trend")
id=idx.mx[1:9]+1
which.max(FFT.har.trd$amp) #4  44yrs
library(itsmr)
# 9 max harmonics
Ocean_solar_anom<-Ocean_solar_anom%>%
  mutate(max.hr.res=hr(har.res,N/idx.mx[id]))
Ocean_solar_anom<-Ocean_solar_anom%>%  mutate(diff=trend.solar-max.hr.res)
# diff=trd.solar-max.har.res)
Ocean_solar_anom%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.solar,col="trd.solar"))+
  geom_line(aes(y=diff,col="differnc"))+
  geom_line(aes(y=har.res,col="har.res"))+
  labs(title="Solar-locked-trend",
       subtitle = "low-pass-filtered extracted periods")
#========
# compare solar locked trend with low-pass filtered
Long_periods=Ocean_solar_anom%>%
  mutate(smth.anomaly=hr(Anomaly,N/1:20))

Long_periods%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=har.res,col="har.res"),linetype = 2)+
  geom_line(aes(y=smth.anomaly,col="smth.anomaly"),linetype = 2)+
  geom_line(aes(y=smth.anomaly+har.res,col="sum"))
# Compare spline fitted trend with baseline trend

Long_periods%>%ggplot(aes(x=dt.mnth,y=trend.solar))+
  geom_line(col=2)+
  geom_line(aes(y=trd7),col=4)+
  geom_line(aes(y=trend.solar-trd7),col=9)+
  labs(title="Compare Ocean Anomaly Trends",
       subtitle="5th° poly(blue),Baseline Trend(red),\ndifference (black) ")
# eliminate season & polynomial trend
#?M3=c("season",12,"season",6,"trend",3)
#?Ocean_solar_anom=Ocean_solar_anom%>%mutate(res3=Resid(anoma.mean,M3))
# 4. Calculate residuals of solar phase locket ocean anomalies
# residual will now fully retain the solar peaks, troughs, and historical minima!
Ocean_solar_anom$Solar_Retained_Residuals <- Ocean_solar_anom$Anomaly - Ocean_solar_anom$trend.solar
Long_periods=Ocean_solar_anom%>%
  mutate(smth.res=smooth.fft(Anomaly,f=0.02))
Long_periods%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Solar_Retained_Residuals),col="grey")+
  geom_line(aes(y=trend.solar),col=2)+
  labs(title = "Decomposition Ocean Anomalies",
       subtitle = "phaselocked to solarpower\n trend fitted to solar power")
saveRDS(Long_periods,"data/Ocean_solar_anom.rds")
# Long_periods==Ocean_solar_anom+smth.res (Solar_Retained_Residuals,f=0.02)
#====================
Long_periods%>%
  mutate(smth.res.11=smooth.fft(Solar_Retained_Residuals,f=0.065),
         smth.res.3.4=smooth.fft(Solar_Retained_Residuals,f=0.02))%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Solar_Retained_Residuals),col="grey")+
  geom_line(aes(y=smth.res.11),col=2)+
  geom_line(aes(y=smth.res.3.4),col=4)+
  labs(x="",title = "Resids of sol.phaselocked-trend",
       subtitle="low-pass filtered  \ncutoff 11.2 years (red) 3.4 years (blue)")

#==================
# dominant long periods from phase locked retained resids and phase locked trend
#   diff trd7-Cleaned_Baseline_trend / smth.res = smooth.fft(Solar_Retained...)
Long_periods=Ocean_solar_anom%>%
  mutate(smth.res=smooth.fft(Solar_Retained_Residuals,f=0.02))
# Analyze trends: Clean_Baseline_Trend/trend.solar)-poly trend 3°

Long_periods%>% mutate(trd.period=trd7-trend.solar,
                       long.sum=trd.period+smth.res)%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=long.sum),col=2)+
  geom_line(aes(y=smth.res),col=3,linetype = 2)+
  geom_line(aes(y=trd.period),col=4,linetype = 2)+
  labs(x="",title = "Sum of dominant harmonics",
       subtitle = "harmonic part of phaselocked trend \n smoothed sol. retained res")
#==================================
