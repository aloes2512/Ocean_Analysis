library(tidyverse)
library(splines)
library(mFilter)  # Or use standard loess() below
library(ggplot2)
Ocean_solar_anom<-readRDS("data/Ocean_solar_anom.rds")
N_global<-NROW(Ocean_solar_anom)
# solar_monthly power variation
solar_mnthly=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))
# find zero crossings of Solar Irradiation SI
find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}
nodes <- find_nodes(solar_mnthly$SI)
# time coordinates of zero-crossings
zero_crossing_times <- solar_mnthly$dt.mnth[find_nodes(solar_mnthly$SI)]

# 1. Set up time bounds and interior zero-crossing knots
start_time <- min(Ocean_solar_anom$dt.mnth, na.rm = TRUE)
end_time   <- max(Ocean_solar_anom$dt.mnth, na.rm = TRUE)

# Keep interior knots strictly inside the data set range
interior_knots <- zero_crossing_times[zero_crossing_times > start_time &
                                        zero_crossing_times < end_time]
# 2. Fit the knot-anchored B-spline linear model
fit_lm <- lm(
  Anomaly ~ bs(dt.mnth,
               knots = interior_knots,
               Boundary.knots = c(start_time, end_time),
               degree = 3),
  data = Ocean_solar_anom
)
# 3. Predict the composite trend
Ocean_solar_anom$trd.solar <- predict(fit_lm)
Ocean_solar_anom<-Ocean_solar_anom%>%mutate(resids_solar=Anomaly-trd.solar)

plt.sol.lckd_trend=Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=resids_solar),col="grey")+
  geom_line(aes(y=trd.solar),col=2)+
  labs(x="",y="anomaly",
       title = "Solar-locked Decomposion",
       subtitle= "trend: phase locked splines")
# 4. Extract the Periodic Component
  ## Using LOESS to strip the smooth secular rise
Ocean_solar_anom$secular_baseline <- predict(
  loess(trd.solar ~ dt.mnth, data = Ocean_solar_anom, span = 0.5)
)
#Ocean_solar_anom$periodic.trnd <- Ocean_solar_anom$trd.solar - Ocean_solar_anom$secular_baseline
Ocean_solar_anom<-Ocean_solar_anom%>%mutate(periodic.trnd=trd.solar-secular_baseline)
# plotting

plt.sol.trend=Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=secular_baseline,col="secular"))+
  geom_line(aes(y=trd.solar,col="trd.sol.lckd"))+
  geom_line(aes(y=secular_baseline-trd.solar,col="periods.sol.lckd"))+
  labs(x="",y="anomaly trend",title = "Solar Locked Trend",
       subtitle = "secular trend separated periods")
#---------
plt.sol.periods=Ocean_solar_anom%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=periodic.trnd))+
  labs(x="",y="period.trend",
       title="Solar Locked Trend \nPeriodic Component")
#-----------
plt.decommp=Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=resids_solar,col= "resids"))+
  geom_line(aes(y=trd.solar,col="trd.solar"),linewidth=1.3)+
  geom_line(aes(y=periodic.trnd,col="periodic trend"),linewidth =1.3)+
  labs(x="",y= "anomalies",
       title = " Solar Locked Anomalies",
       subtitle = "decomposed trend and residuals")
# secular locked resiuduals
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.solar,col="trend.solar"))+
  geom_line(aes(y=secular_baseline,col="sec.base"))+
  geom_line(aes(y=periodic.trnd,col="periodic.trnd"))
Ocean_selcd.sol<-Ocean_solar_anom%>%dplyr::select(dt.mnth,Anomaly,secular_baseline)%>%
  mutate(total.resids=Anomaly-secular_baseline)
Ocean_selcd.sol%>%ggplot(aes(x=dt.mnth))+geom_line(aes(y=total.resids))
library(mgcv)
library(itsmr)
# GAM extend preWWII 1938
full.resds=Ocean_selcd.sol%>%dplyr::select(dt.mnth,total.resids,Anomaly)
train.resds=Ocean_selcd.sol[1:floor(N_global/2),]
# 2. GAM fit  cyclic  Splines (bs = "cc")
gam_raw=gam(total.resids ~ s(dt.mnth,k=180,bs="cc"),data = full.resds)

# k = 90 fits the  1905er trough perfect : cyclical smooth (bs = "cc"),
# fitting with Penalized Maximum Likelihood Estimation 90 parameters say beta 1 t0 beta 90
gam_train <- gam(total.resids ~ s(dt.mnth, k = 90, bs = "cc"), data = train.resds)
gam_predct <- predict(gam_train, newdata = full.resds)
gam_anoma <- gam(Anomaly ~ s(dt.mnth,k=80,bs= "cc"),data = full.resds)
full.resds$gam_Anomaly<-predict(gam_anoma)
# 3. PREDICT = extension of pre-industrial
full.resds$gam_predct <- predict(gam_train, newdata = full.resds)
colnames(full.resds) # "dt.mnth" "total.resids" "gam_predct"  "Anomaly"
full.resds%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Anomaly,col="Anomaly"))+
  geom_line(aes(y=total.resids,col="total"))+
  geom_line(aes(y=gam_predct,col="predct"))
colnames(full.resds)
full.resds%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=gam_Anomaly,col="gm.anoma"))+
  geom_line(aes(y=gam_predct,col="gm.prdct"))
