# "Draft_Report_files/WiP260811.odt"
library(itsmr)
library(tidyverse)
Ocean_solar_anom=readRDS("data/Ocean_solar_anom.rds")
N <- nrow(Ocean_solar_anom) #2116
periods <- N / (1:20)
Mssn=c("season",12,"season",6)
Ocean_solar_anom<-Ocean_solar_anom%>% mutate(anomaly=Resid(Anomaly,Mssn))
# 1. Total Low-Pass Ocean Anomaly Baseline
Ocean_solar_anom<-Ocean_solar_anom%>%mutate(fit_anom_dec=hr(anomaly,periods))
Ocean_solar_anom <- Ocean_solar_anom %>%
  mutate(
    # 1. Total low-pass ocean anomaly baseline
    anom_dec = as.numeric(hr(anomaly, periods)),

    # 2. Low-pass solar trend spectrum (TSI)
    solar_dec = as.numeric(hr(trend.solar, periods)),

    # 3. Low-pass residuals (internal ocean dynamics)
    resids_dec = as.numeric(hr(anomaly - trend.solar, periods)),

    # 4. Detrended solar cycles (TSI minus Richards curve)
    solar_cycles_dec = as.numeric(hr(trend.solar - richards, periods))
  )
# 2. Path A: Low-Pass Solar Trend Spectrum
Ocean_solar_anom<-Ocean_solar_anom%>%
  mutate(solar_dec=
           hr(Ocean_solar_anom$trend.solar, periods))

# 3. Path B: Low-Pass Residuals (Internal Ocean Dynamics)
Ocean_solar_anom<-Ocean_solar_anom%>%
  mutate(resids_dec=hr(anomaly - trend.solar, periods))

# 4. Path C: Detrended Solar Cycles (TSI Trend minus Richards Curve)
Ocean_solar_anom<- Ocean_solar_anom%>%
  mutate(solar_dec=hr(trend.solar - richards, periods))
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=solar_dec,col="solar.har"))+
  geom_line(aes(y=trend.solar,col="trend.sol"))+
  geom_line(aes(y=richards,col="richards"))+
  labs(x="",title = "Decomposition: Solar locked Trend")
plt.residuals=Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=resids_dec,col="resids_dec"))+
  geom_line(aes(y=solar_dec,col="solar_dec"))+
  geom_line(aes(y=solar_dec+resids_dec,col="sum"))+
  labs(x="",title="Sum of Residuals (period > 8.8 y)",
       subtitle = "Res. of sol.lckd + Har. of sol.lckd")
print(plt.residuals)
Ocean_solar_anom<-Ocean_solar_anom%>%mutate(resds.sum=solar_dec+resids_dec)
# --- Variance Partitioning ---
v_total  <- var(Ocean_solar_anom$anom_dec, na.rm = TRUE)
v_solar  <- var(Ocean_solar_anom$solar_dec, na.rm = TRUE)
v_resids <- var(Ocean_solar_anom$resds.sum, na.rm = TRUE)
cov_term <- 2 * cov(Ocean_solar_anom$solar_dec, Ocean_solar_anom$resids_dec, use = "complete.obs")

# Print Summary Table
data.frame(
  Component = c("Solar Trend (TSI)", "Internal Residuals", "Phase Coupling (2*Cov)"),
  Variance = c(v_solar, v_resids, cov_term),
  Pct_Decadal_Var = c(v_solar/v_total, v_resids/v_total, cov_term/v_total) * 100
)
Ocean_solar_anom<-Ocean_solar_anom%>%
  mutate(long.prd.model=richards+resds.sum)
plt.model=Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=long.prd.model))+
  geom_line(aes(y=richards,col="richards"))+
  geom_line(aes(y=resds.sum,col="resds.sum"))+
  labs(x="",y="K",title = "Decomposed Solar Locked Trend ",
       subtitle= "logistic trend (Richards-Curve)+\nlow-passed Residuals")
print(plt.model)
