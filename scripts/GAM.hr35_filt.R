library(mgcv)
library(tidyverse)
# daten
global.ts=readRDS("data/Ocean_solar_anom.rds")
N=NROW(global.ts)
global.hr35<-global.ts|>mutate(hr35=hr(Anomaly,N/1:35))|>
  dplyr::select(dt.mnth,hr35)

plt.hr35<-global.hr35|>ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=hr35))+
  labs(title = "Long Periods of Anommaly",
       subtitle= "periods > 5 years")

# 1. Daten aufteilen in Kalibrierung (vor 1940) und Zukunft
df_train <- global.hr35 %>% subset(dt.mnth <= global.hr35$dt.mnth[N/2])
global.hr35 <- global.hr35
plt.hr35+geom_line(data=df_train,aes(x=dt.mnth,y=hr35,col="train"))+
  geom_line(data=subset(global.hr35,dt.mnth>1938),
            aes(x=dt.mnth,y=hr35,col="full"))+
  labs(x="",title = "Observed Anomaly Filtered & Split")
#--------
# 2. Basis-k train span (90)
k_train <- 90
# 3. length in month training and full ts
n_months_train <- nrow(df_train) # 1058
n_months_full  <- nrow(global.hr35)  # 2116
# 4. set k_full to be equal to training
k_full <- round(k_train * (n_months_full / n_months_train)) # 180

# 2. GAM fit  cyclic  Splines (bs = "cc")
# k = 90 fits the  1905er trough perfect
gam_train <- gam(hr35 ~ s(dt.mnth, k = 90, bs = "cc"), data = df_train)
# 3. predict = extension of pre-industrial
global.hr35$gam_fit <- predict(gam_train, newdata = global.hr35)
global.hr35%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=gam_fit,col="GAM-fit"))+
  geom_line(aes(y=hr35,col="long_periods"),linetype = 2)+
  labs(x="",title="GAM fit Extrapolate vs long Periods\nobserved Anomalies")

# Compare raw gam fitted with train gam fitted
# 2. Basis-k train span (90)
k_train <- 90
# 3. length in month training and full ts
n_months_train <- nrow(df_train) # 1058
n_months_full  <- nrow(global.hr35)  # 2116
# 4. set k_full to be equal to training
k_full <- round(k_train * (n_months_full / n_months_train)) # 180
# 5. both GAMs fit with adapted curvature
gam_model <- gam(hr35 ~ s(dt.mnth, k = k_train, bs = "cc"), data = df_train)

gam_train <- gam(hr35 ~ s(dt.mnth, k = k_train, bs = "cc"), data = df_train)
gam_full  <- gam(hr35 ~ s(dt.mnth, k = k_full,  bs = "cc"), data = global.hr35)

# 6. extract forecast
df_results <- global.hr35 %>%
  mutate(
    # Fit with  ALL Data (adapted k to length)
    fit_total = predict(gam_full, newdata = global.hr35),
    #
    extrapolation_pre1938 = predict(gam_train, newdata = global.hr35),

    difference = fit_total - extrapolation_pre1938
  )
df_results%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=extrapolation_pre1938,color="GAM xtrpol "), linewidth = 0.4) +
  geom_line(aes(y = fit_total, color = "fit_total ")) +
  geom_vline(xintercept = 1938, linetype = "dashed", color = 1) +
  geom_line(aes(y=difference,colour = "Difference"),linewidth =0.7)+
  theme_minimal() +
  labs(title = "Difference: GAM-filtered & extrapolated", y = "Anomalie K", color = "Typ")
#===============
library(itsmr)
# approx difference with polynomial 5°
df_results=df_results%>%mutate(trd.dif=trend(difference,5))
# turning point = max slope
diff(df_results$trd.dif)%>%which.max(.) # 1740
df_results$dt.mnth[1740] # 1994.917
library(zoo)
t0_fixed=df_results$dt.mnth[1740] # was with shorter ts: Oct 1997

df_results%>%  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trd.dif))+
  geom_vline(xintercept=t0_fixed,linetype="dashed")+
  labs(x="",y="hr35 K",title = "Polynomial 5°-Fit to Difference",
       subtitle="difference is total GAM fit - GAM extrapolation TP ~ 1998 ")
# fitting with sigmoid
# Hardcoding t0 based on your observed hr35 minimum == TP difference

trd.data_GAM.hr35=df_results%>%dplyr::select(dt.mnth,difference)
fit_richards <- nls(
  difference ~ L + Delta / (1 + exp(-B * (dt.mnth - t0_fixed)))^(1 / nu),
  data = trd.data_GAM.hr35,
  start = list(L = 0, Delta = 0.9, B = 0.1, nu = 1)
)
coef(fit_richards) #      L          Delta            B           nu
                #-0.004816888  0.622218082  0.127004737   1.615469513
TP=global.hr35$dt.mnth[which.max(diff(predict(fit_richards)))] # 1991.083
yearmon(TP)# Feb 1991 < 1994.917 (estimated t0_fixed)
trd.data_GAM.hr35$pred=predict(fit_richards)
trd.data_GAM.hr35%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=difference),col= "blue")+
  geom_line(aes(y=pred),col="red")+
  geom_vline(xintercept = TP, linetype = "dashed", color = 1) +
  labs(x="",title="Trend: Diff( Observed hr35, Harmonic Forecast)",
       subtitle = "Trend-Curve (red) fitted with 'Richards'")
saveRDS(trd.data_GAM.hr35,"data/trd.hr_GAM.rds")
ggsave("figs/Trend_GAM.hr35_difference.png")
rm(list=setdiff(ls(),"trd.data_GAM.hr35"))
trd.data_GAM.hr35%>%head(2)
N<-NROW(trd.data_GAM.hr35)
FFT.gam.hr=tibble(idx=1:N-1,
                 spc=fft(trd.data_GAM.hr35$difference),
                 amp=Mod(spc),
                 perid=c(0,N/idx[-1])/12)
FFT.gam.hr%>%subset(idx<200)%>%arrange(desc(amp))%>%head(20)
# 9 mx harm idx= 1:9
FFT.maxhr=FFT.gam.hr%>% mutate(spc= ifelse (idx %in% 1:20,spc,0+0i))
library(dplyr)
library(gsignal)
Y.gam.hr=tibble(dt.mnth=trd.data_GAM.hr35$dt.mnth,
                y.max9=Re(ifft(FFT.maxhr$spc)))
Y.gam.hr%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=y.max9))

