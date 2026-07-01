library(mgcv)
library(tidyverse)
# daten
global.ts=readRDS("data/NOAA.ocean.anomalies.rds")
summary(global.ts) # 1850 : 2026 ; Anomaly, Median


N=NROW(global.ts) # 2116
global.ts$dt.mnth[N/2]
# 1. Daten aufteilen in Kalibrierung (vor 1940) und Zukunft
df_train <- global.ts %>% subset(dt.mnth <= global.ts$dt.mnth[N/2])
df_full <- global.ts

# 2. GAM fit  cyclic  Splines (bs = "cc")
# k = 90 fits the  1905er trough perfect
gam_model <- gam(Anomaly ~ s(dt.mnth, k = 90, bs = "cc"), data = df_train)

# 3. predict = extension of pre-industrial
df_full$gam_fit <- predict(gam_model, newdata = df_full)

# Compare raw gam fitted with train gam fitted
# 2. Basis-k train span (90)
k_train <- 90
# 3. length in month training and full ts
n_months_train <- nrow(df_train) # 1058
n_months_full  <- nrow(df_full)  # 2116
# 4. set k_full to be equal to training
k_full <- round(k_train * (n_months_full / n_months_train)) # 180
# 5. both GAMs fit with adapted curvature
gam_train <- gam(Anomaly ~ s(dt.mnth, k = k_train, bs = "cc"), data = df_train)
gam_full  <- gam(Anomaly ~ s(dt.mnth, k = k_full,  bs = "cc"), data = df_full)

# 6. extract forecast
df_results <- df_full %>%
  mutate(
    # Fit with  ALL Data (adapted k to length)
    fit_total = predict(gam_full, newdata = df_full),
    #
    extrapolation_pre1940 = predict(gam_train, newdata = df_full),

    difference = fit_total - extrapolation_pre1940
  )
df_results%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=extrapolation_pre1940,color="GAM xtrpol "), linewidth = 0.4) +
  geom_line(aes(y = fit_total, color = "fit_total ")) +
  geom_vline(xintercept = 1938, linetype = "dashed", color = 1) +
  geom_line(aes(y=difference,colour = "Difference"),linewidth =1)+
  theme_minimal() +
  labs(title = "GAM-filterd & extrapolated", y = "Anomalie K", color = "Typ")
#===============
library(itsmr)
# aprox difference with polynomial 5°
df_results=df_results%>%mutate(trd.dif=trend(difference,5))
# turning point = max slope
diff(df_results$trd.dif)%>%which.max(.) # 1774
df_results$dt.mnth[1774] # 1997.75
library(zoo)
t0_fixed=df_results$dt.mnth[1774] # Oct 1997

df_results%>%  ggplot(aes(x=dt.mnth))+
                geom_line(aes(y=trd.dif))+
                geom_vline(xintercept=1997.75,linetype="dashed")+
  labs(x="",y="anomaly K",title = "Polynomial 5°-Fit to Difference",
       subtitle="difference is total GAM fit - GAM extrapolation TP ~ 1998 ")
# fitting with sigmoid
# Hardcoding t0 based on your observed anomaly minimum == TP difference

trd.data_GAM=df_results%>%dplyr::select(dt.mnth,difference)
fit_richards <- nls(
  difference ~ L + Delta / (1 + exp(-B * (dt.mnth - t0_fixed)))^(1 / nu),
  data = trd.data_GAM,
  start = list(L = 0, Delta = 0.9, B = 0.1, nu = 1)
)
coef(fit_richards) #      L          Delta            B           nu
                   #-0.006840323  0.871218806  0.122671543  1.810572179
TP=global.ts$dt.mnth[which.max(diff(predict(fit_richards)))] # 1992.833
yearmon(TP)# Nov 1992 < 1997.75 (estimated t0_fixed)
trd.data_GAM$pred=predict(fit_richards)
trd.data_GAM%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=difference),col= "blue")+
  geom_line(aes(y=pred),col="red")+
  geom_vline(xintercept = TP, linetype = "dashed", color = 1) +
  labs(x="",title="Trend: Diff( Observed Anomaly, Harmonic Forecast)",
       subtitle = "Trend-Curve (red) fitted with 'Richards'")
saveRDS(trd.data_GAM,"data/trd.data_GAM.rds")
ggsave("figs/Trend_GAM_difference.png")
rm(list=setdiff(ls(),"trd.data_GAM"))
