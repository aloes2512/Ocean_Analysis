# Strucchange modified:
# 1 GAM filter with 6 Basis ffunction
Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
library(tidyverse)
library(itsmr)
library(strucchange)
# Eliminate seasons
M.ssn=c("season",12,"season",6)
Ocean_solar_anomaly<-Ocean_solar_anomaly%>%mutate(anomaly=Resid(Anomaly,M.ssn))
anomaly_df<-Ocean_solar_anomaly%>%dplyr::select(dt.mnth,anomaly)
#Set k approx 22:
gam22.mdl<- gam(data=anomaly_df,anomaly ~ s(dt.mnth, k = 22,  bs = "cc"))
anomaly_df<-anomaly_df%>%mutate(anoma.filt=predict(gam22.mdl))
anomaly_df%>%ggplot(aes(x=dt.mnth))+geom_line(aes(y=anoma.filt))
#-------efp-----
mosum_anomaly.filt <- efp(anoma.filt ~ 1, data = anomaly_df, type = "OLS-MOSUM", h = 0.15)
str(mosum_anomaly.filt)
tms<-seq(0.0747,0.925,length.out=1800)
tibble(times=tms*2116/12+1850,
       values=mosum_anomaly.filt$process)%>%
  ggplot(aes(x=times,y=values))+geom_line()
bp_anomaly <- breakpoints(anoma.filt~ 1, data = anomaly_df, breaks = 1)
tms[1566]
time_seq <- seq(1850, by = 1/12, length.out = 2116)
time_seq[1566] # 1980
plt.sruch_15=anomaly_df%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anoma.filt))+
  geom_vline(aes(xintercept =time_seq[bp_anomaly$breakpoints] ),col=2) #1566

bp_anomaly2 <- breakpoints(anoma.filt~ 1, data = anomaly_df, h=0.2)
plt.sruch_15=anomaly_df%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anoma.filt))+
  geom_vline(aes(xintercept =time_seq[bp_anomaly$breakpoints] ),col=2)+ #1566
  geom_vline(aes(xintercept=time_seq[1693] ),col=3)

