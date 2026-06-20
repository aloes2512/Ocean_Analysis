library(mgcv)
library(tidyverse)
# daten
ts.list=readRDS("data/Area10+.ts.rds")
global.ts=ts.list$global.ts
N=NROW(global.ts) # 2116
global.ts$dt.mnth[N/2]
# 1. Daten aufteilen in Kalibrierung (vor 1940) und Zukunft
df_train <- global.ts %>% subset(dt.mnth <= 1938.083+1/12)
df_full <- global.ts

# 2. GAM fit mit zyklischen Splines (bs = "cc")
# Ein hohes k fängt das 1905er Loch lokal perfekt ab
gam_model <- gam(anomaly ~ s(dt.mnth, k = 90, bs = "cc"), data = df_train)

# 3. Vorhersage für den gesamten Zeitraum
df_full$gam_fit <- predict(gam_model, newdata = df_full)

# Compare raw gam fitted with train gam fitted




# 2. Basis-k für den Trainingszeitraum festlegen (z.B. 90)
k_train <- 90

# 3. Berechne die Anzahl der Monate in beiden Datensätzen
n_months_train <- nrow(df_train) # 1058
n_months_full  <- nrow(df_full)  # 2116

# 4. Skaliere k proportional für die Gesamtdaten

k_full <- round(k_train * (n_months_full / n_months_train)) # 180


# 5. GAMs fit with adapted curvature
gam_train <- gam(anomaly ~ s(dt.mnth, k = k_train, bs = "cc"), data = df_train)
gam_full  <- gam(anomaly ~ s(dt.mnth, k = k_full,  bs = "cc"), data = df_full)

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
  labs(title = "GAM-filtrd & extrpltd", y = "Anomalie K", color = "Typ")
#===============
library(itsmr)
# aprox difference with polynomial 5°
df_results=df_results%>%mutate(trd.dif=trend(difference,5))
# turning point = max slope
diff(df_results$trd.dif)%>%which.max(.) # 1774
df_results$dt.mnth[1774] # 1997.75
library(zoo)
yearmon(1997.75) # Oct 1997
df_results%>%  ggplot(aes(x=dt.mnth))+
                geom_line(aes(y=trd.dif))+
                geom_vline(xintercept=1997.75,linetype="dashed")+
  labs(x="",y="anomaly K",title = "Polynomial Fit to Difference",
       subtitle="difference is total GAM fit - GAM extrapolation TP ~ 1998 ")
# fitting with sigmoid
# Hardcoding t0 based on your observed anomaly minimum (e.g., late 1998)
t0_fixed <- 1997.75
trd.data=df_results%>%dplyr::select(dt.mnth,difference)
fit_richards <- nls(
  difference ~ L + Delta / (1 + exp(-B * (dt.mnth - t0_fixed)))^(1 / nu),
  data = trd.data,
  start = list(L = 0, Delta = 0.9, B = 0.1, nu = 1)
)
coef(fit_richards) #      L          Delta            B           nu
                   #-0.006840323  0.871218806  0.122671543  1.810572179
TP=global.ts$dt.mnth[which.max(diff(predict(fit_richards)))] # 1992.833
yearmon(TP)# Nov 1992
trd.data$pred=predict(fit_richards)
trd.data%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=difference),col= "blue")+
  geom_line(aes(y=pred),col="red")+
  geom_vline(xintercept = TP, linetype = "dashed", color = 1) +
  labs(x="",title="Observed Anomaly - Harmonic Forecast",
       subtitle = "fitted with Richards' Curve (red)")
