library(tidyverse)
Ocean_solar_anom<-readRDS("data/Ocean_solar_anom.rds")
N_global<-NROW(Ocean_solar_anom)
Ocean_solar_anom=Ocean_solar_anom%>%rename("Clean_Baseline_Trend"=trend.solar)

plt.sol.lckd=Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Clean_Baseline_Trend))+
  labs(title = "Solar Locked Basline Trend")
# Define your control parameters
nls_control <- nls.control(maxiter = 1000) # Increase max iterations to 1000
# estimated M
M=1997#  1950 1970 1980 1997
# estimated B
#B=0.12
fit_logistic <- nls(
  Clean_Baseline_Trend  ~ A0+K*(1+exp(-B*(dt.mnth-M)))^(-1/nu),
  data = Ocean_solar_anom,
  start=list(A0=0.000,K=1.2,B= 0.12,nu=1.6),
  control = nls_control
)
# parameters fitted with M= 1970
coefs=coef(fit_logistic)
AIC(fit_logistic) # -5539.545

A0=coefs["A0"] # - 0.125
K=coefs["K"]   # 1.88
nu=coefs["nu"] # 0.3865237
B= coefs["B"]  # 0.02501659
# A optim M B of richards
Rich.mdl=tibble(dt.mnth=Ocean_solar_anom$dt.mnth,
                Rch_t= A0+K*(1+exp(-B*(dt.mnth-1970)))^(-1/nu))
Rich.mdl%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Rch_t))
Ocean_solar_anom$Rch_t<-Rich.mdl$Rch_t
plt.1st.iteration=Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Rch_t))+
  geom_line(aes(y=Clean_Baseline_Trend))+
  labs(x="",title = "Logistic Model & Solar Baseline Trend",
       subtitle= "1st. step: M set 1970  B calculated 0.02,")

fit_M.logistic <- lm(
  Clean_Baseline_Trend   ~ Rch_t,
  data = Ocean_solar_anom)
AIC(fit_M.logistic) # -5152.097
library(broom)
Adapted.logistic<-fit_M.logistic%>%augment()
Adapted.logistic=Adapted.logistic%>%mutate(dt.mnth=Ocean_solar_anom$dt.mnth)
plt.logstc.fit=Adapted.logistic%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=.fitted,col="fittd.logist"))+
  geom_line(aes(y=Clean_Baseline_Trend))+
  labs(x="",y= "K",title = "Solar Locked Trend Fitted",
       subtitle = "two steps fitted by nls")
#B filter loess Clean Baseline trend
mdl.loess.base.trd=loess(Clean_Baseline_Trend ~ dt.mnth,data =Ocean_solar_anom,span=0.9 )
Ocean_solar_anom$loes.baseline<-predict(mdl.loess.base.trd)
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Clean_Baseline_Trend,col="lckd.base.trd"))+
  geom_line(aes(y=loes.baseline,col="secul.lckd"))+
  geom_line(data=Adapted.logistic,aes(x=dt.mnth,y=.fitted,col="logistic"))+
  labs(x="",y="K",title = "Decomposed Locked Basel.Trd",
       subtitle = "secular trend / logistic trd")

Ocean_solar_anom<-Ocean_solar_anom%>%mutate(sol.lckd.res=Anomaly-loes.baseline)
richards_res <- Ocean_solar_anom$Anomaly - predict(fit_M.logistic)
Ocean_solar_anom<-Ocean_solar_anom%>%
  mutate(richards_res=Anomaly-predict(fit_M.logistic),
         vector_res=hr(richards_res,N_global/1:20))
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=vector_res))+
  labs(x="",title = "Low-Pass filtered Difference",
       subtitle = " Difference=Anomaly-Logistic_Trend")
#-----------
library(slider)
vector_res1 <- hr(richards_res, N_global/1:20)
mean(vector_res1) # e-17
vector_obs<-hr(Ocean_solar_anom$Anomaly,N_global/1:20)
rolling_corr2 <- slide_index_dbl(
  .x = 1:length(vector_obs),
  .i = 1:length(vector_obs),
  .f = function(idx) {
    # If the window doesn't have a full 11 years of data, return NA
    if(length(idx) < 132) return(NA)

    # Calculate simple vector correlation on the current slice of indices
    cor(vector_obs[idx], vector_res1[idx])
  },
  .before = 239,   # Looks 131 steps backward + current step = 132 months
  .complete = TRUE # Automatically returns NA for incomplete windows at the start
)
RollingCor=tibble(Time=Ocean_solar_anom$dt.mnth,rolling_corr2=as.numeric(rolling_corr2))#
preind.brk=Ocean_solar_anom$dt.mnth[which(RollingCor$rolling_corr2>0.45& rolling_corr2<=0.5)
]
RollingCor%>%ggplot(aes(x=Time))+
  geom_line(aes(y=rolling_corr2))+
  geom_vline(xintercept = preind.brk,col=2)+
  labs(x="",title="Rolling Correlation of",
       subtitle = "Long Periods SST Anomalies and\nResiduals of adapted Richards")
library(zoo)
brks=as.yearmon(range(preind.brk))
messg=paste0("cor breaks from ",brks[1]," to ",brks[2])
print(messg)

