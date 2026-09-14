Ocean_solar_anom<-readRDS("data/Ocean_solar_anomaly.rds")%>%
  rename("Clean_Baseline_Trend"=trend.solar)%>%
  dplyr::select(dt.mnth,Clean_Baseline_Trend)
# 1. Non-parametric extraction of the purely secular baseline (span 0.9 - 1.0)
library(tidyverse)
library(minpack.lm)
library(broom)
loess_secular <- loess(Clean_Baseline_Trend ~ dt.mnth,
                       data = Ocean_solar_anom,
                       span = 0.9,
                       degree = 2)

# Save the smooth secular curve and multidecadal residuals
Ocean_solar_anom<-Ocean_solar_anom%>%mutate(Secular_LOESS=predict(loess_secular),
                          Multidecadal_Residuals=Clean_Baseline_Trend - Secular_LOESS)
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Clean_Baseline_Trend,col="basel.trd"))+
  geom_line(aes(y=Secular_LOESS,col="secu.loess"))+
  geom_line(aes(y=Multidecadal_Residuals,col="decad.res"))
# 2. Fit Sigmoid/Richards/Tanh model directly to the clean LOESS secular curve
# Multidecadal noise is stripped=>nlsLM converges instantly
# without parameter drift
fit_secular_tanh <- nlsLM(
  Secular_LOESS ~ A0 + (dT / 2) * (1 + tanh((dt.mnth - M) / tau)),
  data = Ocean_solar_anom,
  start = list(A0 = -0.2, dT = 1.0, M = 1985, tau = 40)
)
summary(fit_secular_tanh)
AIC(fit_secular_tanh)
tanh.fit=augment(fit_secular_tanh)
tanh.fit%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=.fitted))+
  labs(x="",title = "Tanh-Clean_Baseline_Trend")
#--------------
Ocean_solar_anom<-Ocean_solar_anom%>%mutate(tanh.trd=tanh.fit$.fitted)
# 1. Extract coefficients from your successful tanh fit
coef_tanh <- coef(fit_secular_tanh) # expects names: A0, dT, M, tau
tanh.data=augment(fit_secular_tanh)
# 2. Derive exact initial parameters for the Richards curve
tanh.reslts <- list(
  A0 = unname(coef_tanh["A0"]), # -0.1257
  K  = unname(coef_tanh["dT"]), # 0.7164489
  M  = unname(coef_tanh["M"]),  # M  = 1998.002
  B  = 2 / unname(coef_tanh["tau"]) # 0.05857891
)
start_richards<-tanh.reslts
start_richards$nu<- 1
# 3. Fit the 5-parameter Richards curve directly to the clean .fitted column

fit_richards_final <- nlsLM(
  .fitted ~ A0 + K * (1 + exp(-B * (dt.mnth - M)))^(-1 / nu),
  data = tanh.data, # or the data frame containing .fitted
  start = start_richards,
  lower = c(A0 = -1.0, K = 0.2, M = 1995, B = 0.001, nu = 0.05),
  upper = c(A0 =  0.5, K = 5.0, M = 2050, B = 0.500, nu = 20.0),
  control = nls.lm.control(maxiter = 200)
)
AIC(fit_richards_final)
coefinal<-coef(fit_richards_final)
# 4. View extracted physical parameters
summary(fit_richards_final)
augment(fit_richards_final)%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=.fitted,col= "Richards"))+
  geom_vline(xintercept = coefinal["M"],linetype = 2)+
  labs(x="",y="lck.-trend",title = "Solar-locked-trend",
       subtitle="limited as sigmoid")
