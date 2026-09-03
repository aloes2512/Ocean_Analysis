sol.power<- Ocean_solar_anomaly$Solar_Variation
log.resids<-Ocean_solar_anomaly$Anomaly-Ocean_solar_anomaly$logistic.trd
log.resids<-log.resids-mean(log.resids)
pow.fit<-lm(log.resids~sol.power,data=Ocean_solar_anomaly,method = "qr")
pow.loess<-loess(log.resids~sol.power,data=Ocean_solar_anomaly,method = "loess")
Ocean_solar_anomaly$pow.fit=predict(pow.fit)
Ocean_solar_anomaly$pow.loess=predict(pow.loess)
Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=pow.loess,col="loess"))+
  geom_line(aes(y=pow.fit,col="lm"))+
  labs(x="",y="ocean anomaly",
       title = "Ocean Anomaly Response to\n Solar Forcing (1850–2020)")
ggsave("figs/Anomaly.solarforced.png")
# comparing AIC
# Extract residual sum of squares and effective parameters
n <- pow.loess$n #2112
rss <- sum(residuals(pow.loess)^2) # 23.81074
enp <- pow.loess$enp # 4.747 effective degrees of freedom

# Gaussian non-parametric AIC
aic_loess <- n * log(rss / n) + 2 * (enp + 1)
aic_loess # -9461.361

# AIC for linear model
rss_lm <- sum(residuals(pow.fit)^2)
k_lm <- length(coef(pow.fit)) # k = 2 (intercept + slope)
aic_lm <- n * log(rss_lm / n) + 2 * (k_lm + 1)

# Compare
c(LM = aic_lm, LOESS = aic_loess)
AIC(pow.fit)

library(mgcv)

# Fit GAM with smooth term on solar power
pow.gam <- gam(log.resids ~ s(sol.power, bs = "cr"),
               data = Ocean_solar_anomaly,
               method = "REML")
Ocean_solar_anomaly$pow.gam=predict(pow.gam)
# 2. Extract components for GAM
n_obs    <- nrow(Ocean_solar_anomaly)# 2112
rss_gam  <- sum(residuals(pow.gam, type = "response")^2) #23.67268
edf_gam  <- sum(pow.gam$edf) # Effective degrees of freedom for GAM 7.385471

# 3. Compute unified non-parametric AIC for GAM
aic_gam  <- n_obs * log(rss_gam / n_obs) + 2 * (edf_gam + 1)

# 4. Compare all three on the exact same scale
aic_results <- c(
  LM    = aic_lm,     # from earlier
  LOESS = aic_loess,  # from earlier
  GAM   = aic_gam
)
aic_results
# compare filtered variables
Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=pow.gam,col="gam"))+
  geom_line(aes(y=pow.loess,col="loess"))+
  geom_line(aes(y=pow.fit,col="lm"))+
  labs(x="",y="ocean anomaly",
       title = "Ocean Anomaly Response to\n Solar Forcing (1850–2020)")
ggsave("figs/Anomaly.solarforced.png")
Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=pow.gam,col="gam"))+
  geom_line(aes(y=pow.loess,col="loess"),linetype = 4)+
  labs(x="",y="ocean anomaly",
       title = "Ocean Anomaly Response to\n Solar Forcing (1850–2020)",
       subtitle= " comparing: gam and loess fit")
ggsave("figs/Anomaly.solar.loess.gam.png")
