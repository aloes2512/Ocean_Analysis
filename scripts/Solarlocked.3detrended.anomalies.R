# fit richards to solar phase locked trend:
# data from line 49:
"~/projects/Ocean_Analysis/scripts/SolarPower.phaselocked_Anomalies.R"
library(tidyverse)
Ocean_solar_anom<-readRDS("data/Ocean_solar_anom.rds")
N_global<-NROW(Ocean_solar_anom)
# check mean
Ocean_solar_anom$Clean_Baseline_Trend%>%mean() # e-16
modern.trend=subset(Ocean_solar_anom,dt.mnth>1975)%>%pull(Clean_Baseline_Trend)
estim.TP=which.min(diff(modern.trend))%>%as.numeric()
estim.TP=1975+estim.TP/12 # 1994.083
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Clean_Baseline_Trend))+
  geom_vline(xintercept=estim.TP,linetype=2,col=2)
# fit logistic to the Clean_Baseline_Trend
baseline.trd=Ocean_solar_anom$Clean_Baseline_Trend
t.ctr=Ocean_solar_anom$dt.mnth-1850
M.ctr=estim.TP-1850
df.ctr=tibble(t.ctr=t.ctr,
              baseline.trd=baseline.trd)
# Define your control parameters
nls_control <- nls.control(maxiter = 1000) # Increase max iterations to 1000
fit_logistic <- nls(
  baseline.trd  ~ A0+K*(1+exp(-B*(t.ctr-M.ctr)))^(-1/nu),
  data = df.ctr,
  start=list(A0=0.000,K=1.2,B=0.12,nu=1.6),
  control = nls_control
)
coefs=coef(fit_logistic)
#second iteration
A0=coefs["A0"]
K=coefs["K"]
B=coefs["B"]
nu=coefs["nu"]
# optim M of richards
fit_M.logistic <- nls(
  baseline.trd  ~ A0+K*(1+exp(-B*(t.ctr-M.ctr)))^(-1/nu),
  data = df.ctr,
  start=list(M.ctr=144),
  control = nls_control
)
AIC(fit_M.logistic) # -5548.318
coef(fit_M.logistic) # 144.0585
M.opt=as.numeric(coef(fit_M.logistic))
Richards=tibble(dt.mnth=Ocean_solar_anom$dt.mnth,
                t.ctr=t.ctr,
  richards =  A0+K*(1+exp(-B*(t.ctr-M.opt)))^(-1/nu) )
plt.rich=Richards%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=richards,col="logistic"))+
  geom_vline(xintercept=1850+M.opt,linetype=2,col=2)
plt.rich+geom_line(data=Ocean_solar_anom,
                   aes(y=Clean_Baseline_Trend),col=4)+
  labs(x="",title = "Solar Phase Locked Trend",
       subtitle="approx with generalized logistic")

richards_res <- Ocean_solar_anom$Anomaly - predict(fit_M.logistic)
vector_res <- hr(richards_res, N_global/1:20)
mean(vector_res) # e-05
vector_obs<-hr(Ocean_solar_anom$Anomaly,N_global/1:20)
tibble(dt.mnth=Ocean_solar_anom$dt.mnth,
       richards_res=richards_res,
       vector_res=vector_res,
       vector_obs=vector_obs)%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=richards_res,col="rich_res"))+
  geom_line(aes(y=vector_res,col="vect_res"))+
  geom_line(aes(y=vector_obs,col="vector_obs"))+
  labs(x="",title = "Global Long-Periods",
       subtitle="richards res & long-periods of richards res")
# rolling cor
# Calculate the rolling correlation over a 132-month (11-year) window
# 'slide_index' passes a slice of the indices to your function
library(slider)
rolling_corr <- slide_index_dbl(
  .x = 1:length(vector_obs),
  .i = 1:length(vector_obs),
  .f = function(idx) {
    # If the window doesn't have a full 11 years of data, return NA
    if(length(idx) < 132) return(NA)

    # Calculate simple vector correlation on the current slice of indices
    cor(vector_obs[idx], vector_res[idx])
  },
  .before = 131,   # Looks 131 steps backward + current step = 132 months
  .complete = TRUE # Automatically returns NA for incomplete windows at the start
)
# plot result
RollingCor=tibble(Time=Ocean_solar_anom$dt.mnth,rolling_corr=as.numeric(rolling_corr))
which(RollingCor$rolling_corr< 0.51 &RollingCor$rolling_corr> 0.49)
Ocean_solar_anom$dt.mnth[1661] # 1988.333
RollingCor%>%subset(Time>1950&Time<1975)%>%pull(rolling_corr)%>%which.min()
1950+90/12
RollingCor%>%subset(Time>1975 & Time<2000)%>%pull(rolling_corr)%>%which.min()
1975+261/12 # 1996.75
RollingCor%>%ggplot(aes(x=Time,y=rolling_corr))+
  geom_line(col="blue")+
  geom_vline(xintercept = 1988.8333)+
  geom_vline(xintercept = 1950+90/12)+
  labs(y="11-Year Rolling Correlation",
       title="Correlation of harmonic fitted observations\n  and logistic detrended residuals")
#======== extrapolation backwards
# alternative breaks
brk1=1950+90/12 #1957
brk2=1988+10/12 #1988
brk3=1975+261/12 # 1996.75
# Preind 1 & 2
Preind1=Ocean_solar_anom%>%subset(dt.mnth<brk1)
N1=NROW(Preind1) # 1290
Preind2=Ocean_solar_anom%>%subset(dt.mnth<brk2)
N2=NROW(Preind2) #1666
Preind3=Ocean_solar_anom%>%subset(dt.mnth<brk3)
N3=NROW(Preind3) #1761
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Anomaly))+
  geom_vline(xintercept = brk1)+
  geom_vline(xintercept = brk2)+
  geom_vline(xintercept = brk3)+
  labs(title="Anomaly Preindustrial Limits",
       subtitle=paste("limits:",round(brk1),",",round(brk2),",",round(brk3)))
# extend back:
# ext backwards
library(gsignal)
extendx_backwards <- function(x, n_back = 3000) {
  N <- length(x)
  N_ext=N+n_back
  t_ext= seq(from=-(n_back-1),N)

  idx_ext=((t_ext - 1) %% N) + 1

  x.ext <- x[idx_ext]

}
signal1=Preind1$Anomaly
N1=length(signal1) # 1290
signal2=Preind2$Anomaly
N2=length(signal2) #1666
signal3=Preind3$Anomaly
sgnx1=extendx_backwards(x=signal1,n_back=3000)
summary(sgnx1)
dates1.ext=seq(1600,by=1/12,length.out=N1+3000)
sgnx2=extendx_backwards(x=signal2,n_back=3000)
dates2.ext=seq(1600,by=1/12,length.out=N2+3000)
sign3.ext=extendx_backwards(x=signal3,n_back=3000)
dates3.ext=seq(1600,by=1/12,length.out=N3+3000)
Y1.ext=tibble(dates=seq(1600,by=1/12,length.out=N1+3000),
       sgn1=extendx_backwards(x=signal1,n_back=3000))%>%
        mutate(sgn1=sgn1-mean(sgn1),
               hrsg1=hr(sgn1,N/1:5))
Y1.ext%>%ggplot(aes(x=dates))+
  geom_line(aes(y=sgn1),col="grey")+
  geom_line(aes(y=hrsg1),col=2)
#======
sgn2=extendx_backwards(x=signal2,n_back=3000)
N2=length(sgn2)
Y2.ext=tibble(dates=seq(1600,by=1/12,length.out=N2+3000),
              sgn2=extendx_backwards(x=signal2,n_back=3000))%>%
  mutate(sgn2=sgn2-mean(sgn2),
         hrsg2=hr(sgn2,N2/1:5))
Y3.ext=tibble(dates=seq(1600,by=1/12,length.out=N3+3000),
              sgn3=extendx_backwards(x=signal3,n_back=3000))%>%
  mutate(sgn3=sgn3-mean(sgn3),
         hrsg3=hr(sgn3,N3/1:5))


Y1.ext%>%ggplot(aes(x=dates))+
  geom_line(aes(y=sgn1),col="grey")+
  geom_line(aes(y=hrsg1),col=2)
Y2.ext%>%ggplot(aes(x=dates))+
  geom_line(aes(y=sgn2),col="grey")+
  geom_line(aes(y=hrsg2),col=2)+
  labs(x="",title = "1988.333 Backcast Anomaly",
       subtitle="periods: 28 to 139 yrs")
Y3.ext%>%ggplot(aes(x=dates))+
  geom_line(aes(y=sgn3),col="grey")+
  geom_line(aes(y=hrsg3),col=2)
