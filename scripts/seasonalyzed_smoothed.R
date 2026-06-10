Ocean_anomaly.phaselocked=readRDS("data/Ocean_anomaly.phaselocked.rds")
library(tidyverse)
library(itsmr)
df=Ocean_anomaly.phaselocked$Ocean_anomaly_detrend
M <- c("season", 12, "season", 6)
df=df%>%mutate(Detrnd.seasonalyzed_anomaly=Resid(Detrended_Anomaly,M))
N=NROW(df)

Df.seasonalyzed=df%>%dplyr::select(Time,Detrnd.seasonalyzed_anomaly)

my_signl=Df.seasonalyzed$Detrnd.seasonalyzed_anomaly
max.idx=function(signl,n){
  N= length(signl)
  FFT=tibble(idx=1:N-1,
             spc=fft(signl),
             amp=Mod(spc))
  idx.mx=FFT%>%arrange(desc(amp))%>%
    subset(idx<N/2)%>% pull(idx)%>% {.[[n]]}
  mx.period=N/idx.mx
  return(mx.period)
}
max.idx(my_signl,1)# 44.95
N/max.idx(my_signl,1) # 46 month  3.83333 years
Df.seasonalyzed=Df.seasonalyzed%>%
  mutate(mx.har=hr(Detrnd.seasonalyzed_anomaly,46))
Df.seasonalyzed%>%ggplot(aes(x=Time,y=mx.har))+geom_line()
Df.seasonalyzed=Df.seasonalyzed%>%
  mutate(har.2nd=hr(Detrnd.seasonalyzed_anomaly,N/max.idx(my_signl,2)))
Df.seasonalyzed%>%ggplot(aes(x=Time))+
  geom_line(aes(y=har.2nd))
Df.seasonalyzed%>%
  mutate(smth.sgnl.01=smooth.fft(Detrnd.seasonalyzed_anomaly,f= 0.01),
         smth.sgnl.02=smooth.fft(Detrnd.seasonalyzed_anomaly,f= 0.02),
         smth.sgnl.015=smooth.fft(Detrnd.seasonalyzed_anomaly,f= 0.015),
         smth.sgnl.03=smooth.fft(Detrnd.seasonalyzed_anomaly,f= 0.03),

  )%>%
  ggplot(aes(x=Time))+
    geom_line(aes(y=smth.sgnl.015))+
    geom_line(aes(y=smth.sgnl.02),col=2)+
   geom_line(aes(y=smth.sgnl.03),col=3)+
  labs(x="",y=" K",title = "Ocean Anomalies Harmonic Periods",
       subtitle=" smoothed.fft cutoff frq 0.015,0.02,0.03 " )
range(Df.seasonalyzed$Detrnd.seasonalyzed_anomaly)

