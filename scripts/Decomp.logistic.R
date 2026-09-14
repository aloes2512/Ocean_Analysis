library(tidyverse)
library(itsmr)
library(strucchange)

#data
Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
# Eliminate seasons
M.ssn=c("season",12,"season",6)
Ocean_solar_anomaly<-Ocean_solar_anomaly%>%mutate(anomaly=Resid(Anomaly,M.ssn))
anomaly_df<-Ocean_solar_anomaly%>%dplyr::select(dt.mnth,anomaly)
#Set k approx 22:
gam22.mdl<- gam(data=anomaly_df,anomaly ~ s(dt.mnth, k = 22,  bs = "cc"))
anomaly_df<-anomaly_df%>%mutate(anoma.filt=predict(gam22.mdl))
Ocean_solar_anomaly<-Ocean_solar_anomaly%>%
  mutate(anoma.filt=predict(gam22.mdl))
Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anoma.filt))
anomaly_df<-anomaly_df%>%mutate(t_vec=dt.mnth-1850)
anomaly_df<-anomaly_df%>%rename("anoma_trend"=anoma.filt,"observed_anomaly"=anomaly)
# 1. Scalar Objective Function (Residual Sum of Squares)
logistic_rss <- function(parms, t_vec, observed, harmonic) {
  T0 <- parms[1]
  r  <- parms[2]
  K  <- parms[3]

  # Safeguard for bounds inside objective
  if (T0 <= 0 || r <= 0 || K <= 0) return(1e10)

  exp_rt <- exp(r * t_vec)
  denom  <- K + T0 * (exp_rt - 1)

  # Safeguard against zero division
  if (any(denom <= 0) || any(!is.finite(denom))) return(1e10)

  pred <- (K * T0 * exp_rt) / denom

  # Calculate RSS
  rss <- sum((observed - pred)^2, na.rm = TRUE)

  if (!is.finite(rss)) return(1e10)
  return(rss)
}

# 2. Starting Parameters and Physical Bounds
start_parms <- c(T0 = 0.005, r = 0.025, K = 1.5)
lower_bounds <- c(T0 = 1e-4,  r = 1e-3,  K = 0.3)
upper_bounds <- c(T0 = 0.10,   r = 0.10,   K = 5.0)

# 3. Fit via L-BFGS-B
fit_optim <- optim(
  par = start_parms,
  fn = logistic_rss,
  t_vec = anomaly_df$t_vec,
  observed = anomaly_df$observed_anomaly,
  harmonic = anomaly_df$anoma_trend,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(maxit = 1000, factr = 1e7)
)

# Extract optimized parameters (order as with par=start_prms)
opt_parms <- fit_optim$par
names(opt_parms) <- c("T0", "r", "K")
print(opt_parms)
#------

model.f=function(t_vec,T0,r,K){
  exp_rt <- exp(r * t_vec)
  denom  <- K + T0 * (exp_rt - 1)
   f=(K * T0 * exp_rt) / denom
   return(f)
}
t_vec<-anomaly_df$t_vec
T0=opt_parms["T0"] # 1e-04
r=opt_parms["r"] # 0.03461522
K=opt_parms["K"]
my.model=tibble(dt.mnth=anomaly_df$dt.mnth,
                logistic=model.f(t_vec,T0,r,K))
my.model%>%ggplot(aes(x=dt.mnth))+geom_line(aes(y=logistic))
Ocean_solar_anomaly<-Ocean_solar_anomaly%>%full_join(my.model)
anomaly_DF=anomaly_df%>%full_join(my.model)
Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anoma.filt,col="anom.trd"))+
  geom_line(aes(y=logistic,col="logistic"))+
  geom_line(aes(y=anoma.filt-logistic,col="diffrnc"))+
  labs(x="",y="anoma",title ="Decomposed Trend" )
# Residuals logistic
Ocean_solar_anomaly%>%mutate(logis.resd=anoma.filt-logistic)%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anomaly-logistic),col="grey")+
  geom_line(aes(y=logis.resd))+
  geom_line(aes(y=trend.solar-logistic),col=2)
