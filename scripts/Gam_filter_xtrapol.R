library(mgcv)
library(tidyverse)
library(itsmr)
# daten
df=readRDS("data/NOAA.ocean.anomalies.rds")
summary(df) # 1850 : 2026 ; Anomaly, Median
global.ts=df$data%>%dplyr::select(dt.mnth,"Anomaly"=anoma.mean)



N=NROW(global.ts) # 2116
global.ts$dt.mnth[N/2] # 1938.083
# 1. Daten aufteilen in Kalibrierung (vor 1940) und Zukunft
df_train <- global.ts %>% subset(dt.mnth <= global.ts$dt.mnth[N/2])
df_full <- global.ts
gam_raw=gam(Anomaly ~ s(dt.mnth,k=180,bs="cc"),data = global.ts)
# 2. GAM fit  cyclic  Splines (bs = "cc")

# k = 90 fits the  1905er trough perfect : cyclical smooth (bs = "cc"),
# fitting with Penalized Maximum Likelihood Estimation 90 parameters say beta 1 t0 beta 90
gam_train <- gam(Anomaly ~ s(dt.mnth, k = 90, bs = "cc"), data = df_train)
gam_full <- predict(gam_train, newdata = df_full)
# 3. predict = extension of pre-industrial
df_full$gam_full <- predict(gam_train, newdata = df_full)
df_full%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Anomaly,col="Anomaly"))+
  geom_line(aes(y=gam_full,col="gam_full"))+
  geom_line(aes(y=Anomaly-gam_full,col="difference"))
df_full<-df_full%>%mutate(lowp.anoma=hr(Anomaly,N/1:20),
                 lowp.diffc=hr((Anomaly-gam_full),N/1:20))
df_full%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=lowp.anoma,col="Anomaly"))+
  geom_line(aes(y=gam_full,col="gam_full"))+
  geom_line(aes(y=lowp.anoma-gam_full,col="difference"))+
  labs(x="",title="Anomaly Decomposed",
       subtitle="low pass filtered,pre 1940 gam extended ")

# Compare raw gam fitted with train gam fitted
# 2. Basis-k train span (90)
k_train <- 90
# 3. length in month training and full ts
n_months_train <- nrow(df_train) # 1058
n_months_full  <- nrow(df_full)  # 2116
# 4. set k_full to be equal to training
k_full <- round(k_train * (n_months_full / n_months_train)) # 180
# 5. both GAMs fit with adapted curvature
#gam_train <- gam(Anomaly ~ s(dt.mnth, k = k_train, bs = "cc"), data = df_train)
library(broom)
# 6. extract forecast
gam_raw%>%augment()%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=.resid,col="residuals"))+
  geom_line(aes(y=.fitted,col="fitted"))+
  labs(x="",title = "Gobal Anomalies \nGAM k=180 fitted")
df_results <- df_full %>%
  mutate(
    # Fit with  ALL Data (adapted k to length)
    fit_total = predict(gam_raw),
    #
    extrapolation_pre1938 = predict(gam_train, newdata = df_full),

    difference = fit_total - extrapolation_pre1938
  )
df_results%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=extrapolation_pre1938,color="GAM xtrpol "), linewidth = 0.4) +
  geom_line(aes(y = fit_total, color = "fit_total ")) +
  geom_vline(xintercept = 1938, linetype = "dashed", color = 1) +
  geom_line(aes(y=difference,colour = "Difference"),linewidth =0.5)+
  theme_minimal() +
  labs(title = "GAM-filterd & extrapolated", y = "Anomalie K", color = "Typ")
#===============
library(itsmr)
# aprox difference with polynomial 5°
df_results=df_results%>%mutate(trd.dif=trend(difference,5))
# turning point = max slope
diff(df_results$trd.dif)%>%which.max(.) # 1750
df_results$dt.mnth[1750] # 1995.75
library(zoo)
t0_fixed=df_results$dt.mnth[1750] # Oct 1997 == 1995.75

df_results%>%  ggplot(aes(x=dt.mnth))+
                geom_line(aes(y=trd.dif))+
                geom_vline(xintercept=1995.75,linetype="dashed")+
  labs(x="",y="anomaly K",title = "Polynomial 5°-Fit to Difference",
       subtitle="difference is total GAM fit - GAM extrapolation TP ~ 1998 ")
# fitting with Richards
# Hardcoding t0 based on your observed anomaly minimum == TP difference

trd.data_GAM=df_results%>%dplyr::select(dt.mnth,difference)
# generalized logistic = Richards curve (Wikipedia)
fit_richards <- nls(
  difference ~ L + Delta / (1 + exp(-B * (dt.mnth - t0_fixed)))^(1 / nu),
  data = trd.data_GAM,
  start = list(L = 0, Delta = 0.9, B = 0.1, nu = 1)
)
coef(fit_richards) #      L          Delta            B           nu
                   #-0.003889104   0.635245515  0.120588599  1.592539744
TP=global.ts$dt.mnth[which.max(diff(predict(fit_richards)))] # 1991.833
yearmon(TP)# Nov 1991 < 1997.75 (estimated t0_fixed)
trd.data_GAM$richard<-predict(fit_richards)
trd.data_GAM%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=difference,col= "difference"))+
  geom_line(aes(y=richard,col="Richards"))+
  geom_vline(xintercept = TP, linetype = "dashed", color = 1) +
  labs(x="",title="GAM-Diff( Observed Anomaly, \nPreindustrial Extended)",
       subtitle = "Trend-Curve (red) fitted with 'Richards'")
saveRDS(trd.data_GAM,"data/trd.data_GAM.rds")
ggsave("figs/Trend_GAM_difference.png")
#--------
# Anomalies Deviation from Richards
global.ts<-global.ts%>%
  mutate(rich.resds=Anomaly-predict(fit_richards),
         rich.resds=rich.resds-mean(rich.resds))
rich.resds<-global.ts$rich.resds

N=length(rich.resds)
global.ts%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=rich.resds))+
  labs(x="",title = "Difference Anomalies Richards Trend")
FFT.rich.rs= tibble(idx=1:N-1,
                    spc=fft(rich.resds),
                    amp=Mod(spc))
which.max(FFT.rich.rs$amp) # idx=2  => 88.16 yrs
FFT.rich.rs%>%subset(idx<140)%>%
  ggplot(aes(x=idx))+geom_point(aes(y=amp))
which(FFT.rich.rs$amp>25&FFT.rich.rs$idx< N/2)
mx.har=FFT.rich.rs%>%subset(idx<140)%>%arrange(desc(amp))%>%pull(idx)
N/(12*mx.har[1:6])

#rm(list=setdiff(ls(),"trd.data_GAM"))
global.ts=global.ts%>%
  mutate(hr1_10=hr(rich.resds,N/1:10),
         hr6_10=hr(rich.resds,N/6:10))
plt.rich.res=global.ts%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=rich.resds),col="grey")+
  geom_line(aes(y=hr1_10,col="perds > 18yr"))+
  geom_line(aes(y=hr6_10,col="18<perds<30yr"))+
  geom_line(aes(y=hr1_10-hr6_10,col="perds < 45yr"))+
  labs(x="",title = "Difference: Anomaly/Richards-Trend",
       subtitle="Long Periods > 18yrs;18/30 yrs from Difference")
print(plt.rich.res)
