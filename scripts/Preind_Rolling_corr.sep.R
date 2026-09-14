library(tidyverse)
Ocean_solar_anom<-readRDS("data/Ocean_solar_anom.rds")
N_global<-NROW(Ocean_solar_anom)
Ocean_solar_anom=Ocean_solar_anom%>%rename("Clean_Baseline_Trend"=trend.solar)
modern.trend=subset(Ocean_solar_anom,dt.mnth>1975)%>%
  pull(Clean_Baseline_Trend)
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Clean_Baseline_Trend))+labs(title = "Solar Locked Basline Trend")
# Define your control parameters
nls_control <- nls.control(maxiter = 1000) # Increase max iterations to 1000
# estimate M
M=1970#  1950 1970 1980 1990
fit_logistic <- nls(
  Clean_Baseline_Trend  ~ A0+K*(1+exp(-0.12*(dt.mnth-M)))^(-1/nu),
  data = Ocean_solar_anom,
  start=list(A0=0.000,K=1.2,M= 1970,nu=1.6),
  control = nls_control
)
coefs=coef(fit_logistic)

#second iteration
A0=coefs["A0"] # - 0.1297
K=coefs["K"]   # 0.804
M=coefs["M"]   # 2023.47
nu=coefs["nu"] # 3.312
B= 0.12
# optim M of richards
Rich.mdl=tibble(dt.mnth=Ocean_solar_anom$dt.mnth,
                  Rch_t= A0+K*(1+exp(-B*(dt.mnth-1995)))^(-1/nu))
Ocean_solar_anom$Rch_t<-Rich.mdl$Rch_t
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Rch_t))+
  geom_line(aes(y=Clean_Baseline_Trend))
fit_M.logistic <- lm(
  Clean_Baseline_Trend   ~ Rch_t,
  data = Ocean_solar_anom)
library(broom)
fit_M.logistic%>%augment()%>%ggplot(aes(x=1:2116))+
  geom_line(aes(y=.fitted))
tidy(fit_logistic)
summary(fit_M.logistic)

library(slider)
richards_res <- Ocean_solar_anom$Anomaly - predict(fit_M.logistic)
Ocean_solar_anom<-Ocean_solar_anom%>%
  mutate(richards_res=Anomaly-predict(fit_M.logistic),
         vector_res=hr(richards_res,N_global/1:20))
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=vector_res))+
  labs(x="",title = "Low-Pass filtered Difference",
       subtitle = " Difference=Anomaly-Logistic_Trend")
vector_res <- hr(richards_res, N_global/1:20)
mean(vector_res) # e-05
vector_obs<-hr(Ocean_solar_anom$Anomaly,N_global/1:20)

vector_obs
vector_res[1:100]
rolling_corr2 <- slide_index_dbl(
  .x = 1:length(vector_obs),
  .i = 1:length(vector_obs),
  .f = function(idx) {
    # If the window doesn't have a full 11 years of data, return NA
    if(length(idx) < 132) return(NA)

    # Calculate simple vector correlation on the current slice of indices
    cor(vector_obs[idx], vector_res[idx])
  },
  .before = 239,   # Looks 131 steps backward + current step = 132 months
  .complete = TRUE # Automatically returns NA for incomplete windows at the start
)
RollingCor=tibble(Time=Ocean_solar_anom$dt.mnth,rolling_corr2=as.numeric(rolling_corr2))#
preind.brk=Ocean_solar_anom$dt.mnth[which(RollingCor$rolling_corr2>0.4& rolling_corr2<=0.5)
]
RollingCor%>%ggplot(aes(x=Time))+
  geom_line(aes(y=rolling_corr2))+
  geom_vline(xintercept = preind.brk,col=2)
# plot optim Richards
Rich.mdl=tibble(dt.mnth=Ocean_solar_anom$dt.mnth,
Rch_t= A0+K*(1+exp(-B*(dt.mnth-1995)))^(-1/nu))
Rich.mdl%>%
  ggplot(aes(x=dt.mnth,y=Rch_t))+geom_line()

