Ocean_anoma=readRDS("data/NOAA.OCEAN.ANOMALIES.rds")$data #2120 obs
library(tidyverse)
N=NROW(Ocean_anoma) #2120
colnames(Ocean_anoma)
# eliminate annual and semianual (at equator)
library(itsmr)
M=c("season",12,"season",6)
Ocean_anomaly=Ocean_anoma%>%
  mutate(anomaly=Resid(anoma.mean,M))
fit.loess<-loess(anomaly ~ dt.mnth, span=0.25,data = Ocean_anomaly)
Ocean_anomaly<-Ocean_anomaly%>%
  mutate(secular=predict(fit.loess))
# combine solar_monthly with ocean anomaly
solar_mnthly=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))
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
trd.solar <- lm(anomaly ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                data = Ocean_anomaly)
Ocean_anomaly<-Ocean_anomaly%>%mutate(trend.solar=predict(trd.solar))
Ocean_solar_anomaly=Ocean_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,anomaly,"Solar_Variation"=SI,secular,trend.solar)%>%
  mutate(trd_hr.res=trend.solar-secular)

plt.sol.trd=Ocean_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_point(aes(y=anomaly),col="grey",size=0.3)+
  geom_line(aes(y=secular,col="trend"))+
  labs(x="",y="anomaly",title = "Global Ocean Anomaly",
       subtitle="secular = loess.smth span=0.25 ")
# combine solar_monthly with ocean anomaly
solar_mnthly=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))

# combine solar and ocean data
Ocean_solar_anomaly=Ocean_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,anomaly,"Solar_Variation"=SI,secular)
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
trd.solar <- lm(anomaly ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                data = Ocean_solar_anomaly)
Ocean_solar_anomaly$trend.solar<-predict(trd.solar)
Ocean_solar_anomaly=Ocean_solar_anomaly%>%
  mutate(trd_hr.res=trend.solar-secular)
plt_trend.solar<-Ocean_solar_anomaly%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.solar,col="trd.sol."))+
  geom_line(aes(y=secular,col="secular"))+
  geom_line(aes(y=trd_hr.res,col="har.part"))+
  labs(x="",title="Compare Solar-locked-trend\nSecular Trend")

#ggsave("figs/sol.trds.png")
#saveRDS(Ocean_solar_anomaly,"data/Ocean_solar_anomaly.rds")
#--------
# 3. Extract the stable Baseline Trend and visualise
plt.decomp_solar.lckd<-Ocean_solar_anomaly%>%
  mutate(solar.res=anomaly-trend.solar)%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=solar.res),col="grey")+
  geom_line(aes(y=trend.solar,col="trend.solar"))+
  labs(title="Solar-locked Trend",
       subtitle = "extracted with spline\n fit to solar power")
#===========
# periods of trend.solar estimate
trd_hr.res <-Ocean_solar_anomaly$trd_hr.res
FFT.har.trd=tibble(idx=1:N-1,
                   spc=fft(trd_hr.res),
                   amp=Mod(spc))
FFT.har.trd%>%subset(idx<20)%>%
  ggplot(aes(x=idx))+
  geom_point(aes(y=amp))+
  labs(x="rel.frequency",
       title="Fourier Coef. Harmonic Trend")
FFT.max<-tibble(idx=1:N-1,spc=rep(0+0i,N),amp=Mod(spc))
which.max(FFT.har.trd$amp) #8 22yrs 530 mnth
idx.mx=FFT.har.trd%>%subset(amp>5.1)%>% pull(idx)
FFT.max[idx.mx,]<-FFT.har.trd[idx.mx,]
library(gsignal)
plt.har.comp=tibble(dt.mnth=Ocean_solar_anomaly$dt.mnth,
       y.mx=Re(ifft(FFT.max$spc)))%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=y.mx))+
  labs(x=NULL,title="Harmonic Component of\n Solar Locked Trend")
library(itsmr)
# 6 max harmonics
Ocean_solar_anomaly%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.solar,col="trd.solar"))+
  geom_line(aes(y=Re(ifft(FFT.max$spc)),col="har.part"))+
  geom_line(aes(y=trd_hr.res,col="trd_hr.res"))+
  labs(title="Solar-locked-trend",
       subtitle = "low-pass-filtered extracted periods")
#========
