Ocean_NOAA.data=readRDS("data/NOAA.ocean.anomalies.rds")$data
M=c("season",6,"season",12)
Ocean_NOAA.data=Ocean_NOAA.data%>%mutate(anoma.mean=Resid(anoma.mean,M),
                         trd3=trend(anoma.mean,3))
Ocean_anoma.mean=Ocean_NOAA.data%>%dplyr::select(dt.mnth,anoma.mean,trd3)
Areas.ts=readRDS("data/Areas.ts.rds")
AMO=Areas.ts$amo
colnames(AMO)
N=NROW(AMO) #2116
library(itsmr)
M=c("season",12,"season",6)
AMO=AMO%>%  mutate(amo.anoma=amo_ts-mean(amo_ts),
                             anomaly=Resid(amo.anoma,M),
                                     trd3=trend(anomaly,3),
                                     res3=Resid(anomaly,3),
                           Year=year(t),mnth=month(t),
                           dt.mnth=Year+(mnth-1)/12)
AMO_anomaly=AMO%>%dplyr::select(dt.mnth,anomaly,trd3)
AMO.plt=AMO_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anomaly),col="green")+
  geom_line(aes(y=trd3),col=2,linewidth=1.3)
AMO.plt+geom_line(data = Ocean_NOAA.data,
                  aes(x=dt.mnth,y=trd3),col=4)+
  geom_line(data = Ocean_NOAA.data,
             aes(x=dt.mnth,y=anoma.mean),size=0.2)+
  labs(x="",title = "Ocean-, AMO-anomalies ",
       subtitle="Global (black),AMO (green) and poly trend 3")
ggsave("figs/Global_AMO.anomalies.png")
#=========================
# lock to solar power
solar_mnthly=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))%>%dplyr::select(-TSI)
# combine solar and amo data
AMO_solar_anom=AMO_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,anomaly,"Solar_Variation"=SI,trd3)
Global_solar_anom=Ocean_anoma.mean%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  rename("anomaly"=anoma.mean)

# anchor timeseries solar power
library(splines)
# find zero crossings of SI
find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}
nodes <- find_nodes(solar_mnthly$SI)
# 1. Identify the time coordinates of your zero-crossings
zero_crossing_times <- solar_mnthly$dt.mnth[find_nodes(solar_mnthly$SI)]

# 2.a Fit a continuous spline locked at these exact nodes
# spline degree = 3 creates a smooth, continuous curve through the zero-crossings
trend_fit <- lm(anomaly ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                data = AMO_solar_anom)
# 2.b Fit splines locked at solar nodes to global
trend_fit.global <- lm(anomaly ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                data = Global_solar_anom)
Global_solar_anom$Baseline_Trend=predict(trend_fit.global)
# 3. Extract the stable Baseline Trend and visualise
AMO_solar_anom$Baseline_Trend <- predict(trend_fit)
AMO_Baseline_Trend.plt=AMO_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Baseline_Trend ),col=2)
AMO_Baseline_Trend.plt+ geom_line(data=Global_solar_anom,
            aes(x=dt.mnth,y=Baseline_Trend ),col=4)+
  labs(x="",title = "Solar Locked Baseline Trends",
       subtitle = "Global(blue)/AMO(red) ")
ggsave("figs/Global.AMO_lcked.trends.png")
