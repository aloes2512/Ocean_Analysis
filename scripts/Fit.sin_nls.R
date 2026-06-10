#{r fit secular Har + linear}
#data used for extension
list.files(here("data"), recursive = TRUE)
Ocean_anoma=readRDS("data/ocean.averaged_monthly.temperatures.rds")
Ocean_anoma%>% subset(dates>"1981-12-01")%>%
  dplyr::select(ocean_anom)%>%range() # -0.0008980202  0.0044235446
readRDS("data/NOAA.ocean.anomalies.rds")%>%
  dplyr::select(mean.anom)%>%range # -0.48  0.56
cor.fact= -0.48/ -0.0008980202
# scale with previous timeseries values from NOAA.oceananomalies:-0.48 /-0.0008980202
Ocean_anoma=Ocean_anoma%>%mutate(ocean_anom.cor=ocean_anom*(-0.48 /-0.0008980202))

# data from 1981 to 2023

colnames(Ocean_anoma) # dates  ocean_anom "ocean_anom.cor"
M.ssn=c("season",12)
summary(Ocean_anoma)
library(itsmr)
Ocean_anoma=Ocean_anoma%>%mutate(yr.anom=Resid(ocean_anom.cor,M.ssn),
                                 yr.smth=smooth.fft(yr.anom,0.01),
                                 res3=Resid(yr.smth,3),
                                 trd3=trend(yr.smth,3))
Ocean_anoma%>%ggplot(aes(x=dates,y=yr.smth))+geom_line(col=2)
#--------
colnames #"dates"          "ocean_anom"     "ocean_anom.cor" "yr.anom"
          #[5] "yr.smth"        "res3"           "trd3"
Ocean_trends=Ocean_anoma%>%
  mutate(yr.anom=Resid(ocean_anom.cor,M.ssn),
         yr.smth=smooth.fft(yr.anom,0.1),
         trd1=trend(yr.anom,1),
         trd2=trend(yr.anom,2),
         trd3=trend(yr.anom,3),
         trd4=trend(yr.anom,4))

Ocean_trends%>%ggplot(aes(x=dates,y=yr.smth))+geom_line()
#_______
library(tidyverse)
N= NROW(Ocean_trends)
trd.slope=diff(range(Ocean_trends$trd1))/N
Amp=sd(Ocean_trends$yr.anom)*sqrt(2) # 0.0014
T.trd=Amp*2*pi/trd.slope
t.mnth0=N/2 # -1803
dates=Ocean_trends$dates
dates[N/2] # 1953-01-15
phi1=0

fit_secular <- nls(
  sgnl ~ Amp* sin((2*pi/T) * t + phi) + C+ d*t,
  data = list(sgnl = Ocean_anoma$yr.smth, t = (1:N) - 1),
  start = list( T=T.trd, phi = phi1, C = 0,d=3e-7),
  control = nls.control(maxiter = 1000)
)
AIC(fit_secular) # --523.2194
T1=coef(fit_secular)["T"] # 3769
phi1=coef(fit_secular)["phi"] # 2.615
coefs=coefficients(fit_secular)
coefs
coefs[[1]]# 3769 month-> 314 years

fit_amp <- nls(
  sgnl ~ -0.47* sin((2*pi/T) * t + phi1) + C+ d*t,
  data = list(sgnl = Ocean_trends$yr.smth, t = (1:N) - 1),
  start = list(T=T1,  C = coef(fit_secular)["C"],d=coef(fit_secular)["d"]),
  control = nls.control(maxiter = 1000)
  )
coefs.A=coefficients(fit_amp)
coefs.A
coefs.A[[1]]# - 0.473 old was
AIC(fit_amp) # -63.6


