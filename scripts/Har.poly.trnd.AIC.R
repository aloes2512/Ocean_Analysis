# Har opti from:
path="~/R-Studio Directory/Ocean_Periods_median.qmd"
Har.opti%>%ggplot(aes(x=dates))+geom_line(aes(y=y.opti),col=2,linewidth=1.5)+
  geom_line(data=Ocean_anoma,aes(y=yr.anom),col=3)+
  labs(title = "Long Period Harmonic Fit",
       subtitle="long period 3637.885 month")

#{r nls}
# BEST:model with linear component # -22459.87
# has smalest AIC # -22459.87
fit_secular <- nls(
  sgnl ~ Amp* sin((2*pi/T) * t + phi) + C+ d*t,
  data = list(sgnl = Ocean_anoma$yr.anom, t = (1:N) - 1),
  start = list( T=T.trd, phi = phi1, C = 0,d=3e-7),
  control = nls.control(maxiter = 1000)
)
AIC(fit_secular) # -22459.87
coefs=coefficients(fit_secular)
coefs
coefs[[1]]#  2924.556 244 years
#====
# 2nd best:no linear term but + constant value C
#model with C ## -21721.12
fit_secular2 <- nls(
  sgnl ~ Amp* sin((2*pi/T) * t + phi)+C ,
  data = list(sgnl = Ocean_anoma$yr.anom, t = (1:N) - 1),
  start = list( T=T.trd, phi = phi1,C=0),
  control = nls.control(maxiter = 1000)
)
AIC(fit_secular2)# -21721.12
coef2=coefficients(fit_secular2)
coef2[[1]]/12 # 303.157 years
#================
#====

# 3rd model (only harmonic)
# w/o constant # AIC -21186.3
fit_secular3 <- nls(
  sgnl ~ Amp* sin((2*pi/T) * t + phi) ,
  data = list(sgnl = Ocean_anoma$yr.anom, t = (1:N) - 1),
  start = list( T=T.trd, phi = phi1),
  control = nls.control(maxiter = 1000)
)
coefs3=coefficients(fit_secular3)
coefs3
coefs3[[1]]/12#  4881.688 = 406 years
AIC(fit_secular3) # -21186.3
#====
# iteration with optim values from 2nd as constant
# and optim amplitude only
fit_secular2 <- nls(
  sgnl ~ Amp* sin((2*pi/T) * t + phi)+C ,
  data = list(sgnl = Ocean_anoma$yr.anom, t = (1:N) - 1),
  start = list( T=T.trd, phi = phi1,C=0),
  control = nls.control(maxiter = 1000)
)
AIC(fit_secular2)# -21721.12
coef2=coefficients(fit_secular2)
coef2[[1]]/12 # 303.157 years
#---------
Amp=sd(LinRes$yr.anom)*sqrt(2) # 0.0014
fit_secular.amp <- nls(
  sgnl ~ Amp* sin((2*pi/coef2[[1]]) * t + coef2[[2]])+coef2[[3]] ,
  data = list(sgnl = Ocean_anoma$yr.anom, t = (1:N) - 1),
  start = list( Amp=0.0014),
  control = nls.control(maxiter = 1000)
)
AIC(fit_secular.amp) # -21751.7
coefficients(fit_secular.amp)# 0.001551118
#============
#as before but having C also as variable to optim
fit_secular.amp.C <- nls(
  sgnl ~ Amp* sin((2*pi/coef2[[1]]) * t + coef2[[2]])+C ,
  data = list(sgnl = Ocean_anoma$yr.anom, t = (1:N) - 1),
  start = list( Amp=0.0014,C=0),
  control = nls.control(maxiter = 1000)
)
AIC(fit_secular.amp.C )# -21763.49
# adjusting phi also
fit_secular.amp.phi.C <- nls(
  sgnl ~ Amp* sin((2*pi/coef2[[1]]) * t + phi)+C ,
  data = list(sgnl = Ocean_anoma$yr.anom, t = (1:N) - 1),
  start = list( Amp=0.0014,phi=coef2[[2]],C=0),
  control = nls.control(maxiter = 1000)
)
AIC(fit_secular.amp.phi.C) # -21799.74
coef.2nd=coefficients(fit_secular.amp.phi.C)
# try additiona iteration to immprove T singular gradient matrixat initial parameters
fit_secular.3rd <- nls(
  sgnl ~ coef.2nd[[1]]* sin((2*pi/T) * t + phi)+coef.2nd[[3]],
  data = list(sgnl = Ocean_anoma$yr.anom, t = (1:N) - 1),
  start = list( T=3000,Amp=0.0014,phi=0),
  control = nls.control(maxiter = 1000)
)
N=length(Ocean_anoma$yr.anom) # 1754
FFT=tibble(idx=1:N,spc=Ocean_anoma$yr.anom%>%fft(),amp=Mod(spc))
idx.ord=FFT%>%arrange(desc(amp))%>%pull(idx)
idx.mx= idx.ord[1:12]
library(gsignal)

FFT_filtered <- FFT %>%
  mutate(spc = if_else(idx %in% idx.mx, spc, 0))
tibble(t=1:N,y.mx=Re(ifft(FFT_filtered$spc)))%>%
  ggplot(aes(x=t,y=y.mx))+geom_line()
tibble(t=1:2*N,y.mx=Re(ifft(FFT_filtered$spc)))%>%
  ggplot(aes(x=t,y=y.mx))+geom_line()
#extending the function 1.5 times in the past
N_ext <- round(1.5 * N)  # Total length: original N + 50% past
t_past <- seq( -0.5 * N + 1,0, by = 1)  # Past: 0, -1, ..., -(0.5N-1)
t_ext <- c(t_past, 1:N)  # Combined: past + original
library(dplyr); library(ggplot2)

y_ext <- Re(ifft(FFT_filtered$spc)[(t_ext - 1) %% N + 1])
tibble(t = t_ext, y.mx_ext = y_ext) %>%
  ggplot(aes(x = t, y = y.mx_ext)) +
  geom_line()
