Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
N=NROW(Ocean_solar_anomaly) #2116
library(itsmr)
M.ssn=c("season",12,"season",6)
df_full<-Ocean_solar_anomaly%>%
  mutate(anomaly=Resid(Anomaly,M.ssn))%>%dplyr::select(dt.mnth,anomaly)
df_train<-df_full[1:floor(N/2),]
gam_train <- gam(anomaly ~ s(dt.mnth, k = 90, bs = "cc"), data = df_train)
gam_filt <- gam(anomaly~s(dt.mnth,k= 180,bs="cs"),data= df_full)
df_full[floor(N/2),] # 1938. , -0.123

extrapolation_pre1938<-predict(gam_train, newdata = df_full)
df_results <- df_full%>%
  mutate(
    # Fit with  ALL Data (adapted k to length)
    signl = predict(gam_filt),
    poly7 = trend(signl,7),
    #
    extrapolation_pre1938 = predict(gam_train, newdata = df_full),

    difference = signl - extrapolation_pre1938
  )
df_results%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=extrapolation_pre1938,col="extrapol"))+
  geom_line(aes(y=difference,col="differnc"))+
  labs(x="",title = "Anomaly Pre 1938 @ Difference",
       subtitle="method  extrapolated GAM ")
df_results%>%ggplot(aes(x=dt.mnth,y=signl))+
  geom_line(col="grey")+
  geom_line(aes(y=poly7))+
  geom_vline(xintercept = 1997)
# fit logistic
# 1. Fit standard 4-parameter logistic model to the 1-2 SSA trend
fit_log <- nlsLM(
  difference ~ L / (1 + exp(-k * (dt.mnth - t0))) + b,
  data = df_results,
  start = list(L = 1.0, k = 0.03, t0 = 1990, b = -0.2)
)
round(coef(fit_log),3) #  L       k       t0        b
                      # 0.629    0.107 1988.588   -0.003
df_results$logist<-predict(fit_log)
df_results%>%ggplot((aes(x=dt.mnth)))+
  geom_line(aes(y=logist))+
  geom_line(aes(y=difference))+
  geom_line(data=df_full,aes(x=dt.mnth,y=anomaly),col="grey")
# fit richards
## step 1 estim t0 to from plot, and optim  with t0 fixed
t0_fixed= 1988.588 # estimated from logistic
fit_anomaly <- nls(
  sgnl ~ L + Delta / (1 + exp(-B * (dt.mnth - t0_fixed)))^(1 / nu),
  data =df_full,
  start = list(L = 0, Delta = 0.9, B = 0.1, nu = 1)
)
coef1=coef(fit_anomaly)
## step2 t0 variable
L=coef1["L"] # -0.1251381
Delta=coef1["Delta"] # 1.8828
B=coef1["B"] # 0.0250
nu=coef1["nu"] # 0.3868
fit_anomaly_2<- nls(
  sgnl ~ L + Delta / (1 + exp(-B * (dt.mnth - t0)))^(1 / nu),
  data = df_full,
  start = list(t0=t0_fixed)
)
coef(fit_anomaly_2) # 1996.973
Rich_final=tibble(dt.mnth=df_full$dt.mnth,
  richrds=predict(fit_anomaly_2),
  sgnl=predict(gam_full))
ceilng=L+(Delta^-nu)# 0.6031968
Rich_final%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=sgnl),col="grey")+
  geom_line(aes(y=richrds,col="Richard"))+
  geom_vline(xintercept = coef(fit_anomaly_2),
             linetype = 2)+
  geom_hline(yintercept = ceilng,linetype= 3)+
  labs(x="",y="trends",title="Anomalies fitted to Richards")
ceilng=L+(Delta^-nu)# 0.6031968
#------------
Rich_final%>%mutate(Resds=sgnl-richrds)%>%
  ggplot(aes(x=dt.mnth))+geom_line(aes(y=Resds))
