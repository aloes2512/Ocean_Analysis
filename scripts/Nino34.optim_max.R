require(tidyverse)
require(itsmr)
require(zoo)
NINO.34=readRDS("data/NINO3.4.rds")
# moving average as low-pass filter
NINO.34=NINO.34%>%mutate(T.smth=smooth.ma(tmp_C,12))
NINO.34%>%ggplot(aes(x=yr.mnth))+
  geom_line(aes(y=T.smth))
# itsmr model detrend and remove seasonal
M= c("trend",4,"season",12)
NINO.res=NINO.34%>%
  mutate(lin.res=Resid(T.smth,M))
NINO.res%>%ggplot(aes(x=yr.mnth))+
  geom_line(aes(y=lin.res))+
  labs(x="",y= "residuals [K]",
       title = "NINO 3.4 Temper.- Anomaly",
       subtitle = "12 mnth. moving average, polyn. trend 4,
  seasonal trend removed")
summary(NINO.res)
N=NROW(NINO.res)
f0=1/N
colnames(NINO.res)
#===============
#Fourier transform
FFT=tibble(idx=1:N,spc=fft(NINO.res$lin.res),amp=Mod(spc))
idx.max.4=FFT%>%arrange(desc(amp))%>%head(8)%>%pull(idx) # 13 28 18
FFT[idx.max.4,]
FFT.filt=tibble(idx=1:N,spc= rep(0i,N),amp=Mod(spc))
FFT.filt[idx.max.4,]=FFT[idx.max.4,]
tibble(t=NINO.34$yr.mnth,y.4max=ifft(FFT.filt$spc))%>%
  ggplot(aes(x=t,y=y.4max))+geom_line()+
  labs(x="",y= "sum 4 max [ΔΚ]",
       title = "NINO 3.4 Periodic Trend",
       subtitle = "deviation polynom 4. degree & seasonal")
# select the 2 longest periods
NINO.detrend=NINO.res%>% mutate(hr1=hr(lin.res,N),
                                  hr2=hr(lin.res,N/2)
                                  )
NINO.detrend%>%ggplot(aes(x=yr.mnth))+
  geom_line(aes(y=hr1+hr2),col=4)+
  geom_line(aes(y=hr1),col=2)+
  geom_line(aes(y=hr2),col=3)+
  labs(title ="NINO 3.4 2 longest periods + sum  ")
#smooth by annual average
NINO.yr=NINO.res%>%mutate(yr=floor(as.numeric(yr.mnth)))%>%
group_by(yr)%>%
mutate(Tmean.yr=mean(lin.res))%>%ungroup()
NINO.yr%>%ggplot(aes(x=yr.mnth))+
geom_line(data = NINO.yr,aes(y=Tmean.yr),col="red")+
labs(x="mnth",y = " T [°C]",
     title = "NINO 3.4 Mean Annual Temperatures")
NINO.trnd.yr=NINO.yr%>%mutate(yr.trnd=trend(Tmean.yr,4))
NINO.trnd.yr%>%
  ggplot(aes(x=yr.mnth,y=yr.trnd))+geom_line()
# smooth by moving average on noisy signal
noisy.signl=NINO.res$lin.res
NINO.res$signl.smth=NINO.res$lin.res%>%smooth.ma(12)
chk.lin=NINO.res%>%mutate(linear=trend(signl.smth,1))
chk.lin$linear%>%range() # -0.005428084  0.002435980
chk.lin$linear%>%mean() # -0.001496052
#  dominant frequencies by fourier
## from smoothed signal
FFT=fft(NINO.res$signl.smth)
amp.mx=FFT[1:(N%/%2)]%>%Mod()%>%sort(decreasing = T)%>%head(4) #[1] 216.0108 203.4824
idx.mx=which(Mod(FFT[1:(N%/%2)])%in% amp.mx) # 9 13 18 28
frq.mx=idx.mx-1 #  8 12 17 27

# periods of dominant frequencies
N/frq.mx #  231.62500  154.41667 109.00000  68.62963 month
         #  19.302083  12.868056   9.083333  5.719136  years
N/frq.mx[1] # 231.625  month
N/frq.mx[2] # 154.4167 month
#============
## filter dominant periods from signl.smoothed
NINO.smth.mx=NINO.res%>% mutate(hr8=hr(signl.smth,N/8),
                                hr12=hr(signl.smth,N/12),
                                hr17=hr(signl.smth,N/17),
                                hr27=hr(signl.smth,N/27))

NINO.smth.mx%>% ggplot(aes(x=yr.mnth))+
                  geom_line(aes(y=hr27),col=4)+
                  geom_line(aes(y=hr17),col=3)+
                  geom_line(aes(y=hr12),col=2)+
                  geom_line(aes(y=hr8+hr12+hr17+hr27))
#=========
# sum of two longest and two strongest harmonics
# periods N, N/2. N/27, N/41
NINO.res=NINO.res%>% mutate(hr1=hr(signl.smth,N),
                   hr2= hr(signl.smth,N/2),
                   hr8=hr(signl.smth,N/8),
                   hr12=hr(signl.smth,N/12),
                   hr17=hr(signl.smth,N/17),
                   hr27=hr(signl.smth,N/27),
                   sum_hr=hr1+hr2+hr8+hr12+hr17+hr27
)
NINO.res%>%ggplot(aes(x=yr.mnth))+
            geom_line(aes(y=sum_hr))+
            labs(x="",y= "ΔΚ",title="NINO 3.4 Sum of Harmonics",
            subtitle="2 longest periods + 4 strongest periods")
#===============================
# calculate max ampl periods by maximum likelihood
#========= model function

mdl.fun= function(frequency,phase,time,amplitude){amplitude * cos(2 * pi * frequency * time + phase)}
# ============ loglike function
log.lik= function (signl,f,phi){
  signl=(signl-mean(signl)/sd(signl)) # normalize to μ=0 and sd=1
  time= as.numeric(NINO.34$yr.mnth)
  RR=signl- mdl.fun(frequency=f,phase=phi,signl,amplitude = sd(signl))      #resids == diff of norm.signl and model function
  ll=-sum(dnorm(RR,0,1,log = T))
  return(ll)
}

log.lik(signl = NINO.res$lin.res,f=frq.mx[1]/N,phi = 0)
# range for f and phi optimization
frq.mx # 8 12 17 27
frq.mx=frq.mx/N # 27/N  estimated max harmonic
f.seq=c(seq(0.8*frq.mx[1],frq.mx[1],length=100)[-100],
        seq(frq.mx[1],frq.mx[1]*1.2,length=100)) # seq of f symmetric to frq.mx[1]
lke.grid=expand.grid(f=f.seq,
                     ph=seq(-pi-.5,pi,length=100))
range(lke.grid$f) #  0.003453859 0.005180788
lke.grid$f[100]== frq.mx[1]# TRUE;  0.01457097
frq.mx[1] # start value 0.004317323
signl=NINO.res$signl.smth
#=========
Likes=lke.grid%>%mutate(like=map2_dbl(f,ph,log.lik,signl=signl))
colnames(Likes) # "f"    "ph"   "like"
#============
range(Likes$like) # 1889.245 2088.361
which.max(Likes$like) # 19702
Likes[which.max(Likes$like),]
max.prms=Likes[which.max(Likes$like),]%>%as.numeric()
'          f             ph      like
 19702 0.003453859 3.141593 2088.361'

# re-construct dominant frequency
#mdl.fun(frequency,phase,time,amplitude)
Y.recst=tibble(t=as.numeric(NINO.34$yr.mnth),
             y=mdl.fun(frequency = max.prms[1],
                       phase=max.prms[2],
                       time=t,amplitude = 1))
Y.recst%>%ggplot(aes(x=t,y=y))+geom_line()+
  labs(x="",y = "ΔΚ", title = "NINO 3.4 Dominant Period",
       subtitle = paste( round(1/max.prms[1]), "month"))
#======================
# View Likes in grid
Likes%>%colnames() #[1] "f"    "ph"   "like"
range(Likes$like) #
Likes[which.max(Likes$like),] # 19702.......
range(lke.grid$f) #
f1<-Likes[which.max(Likes$like),1] #0.003453859

L2=subset(Likes,f!=f1)
L2[which.max(L2$like),] # 19603
f2=L2[which.max(L2$like),1] # 0.00346258
Max2=bind_rows(Likes[19702,],L2[19603,])
Likes%>%
  ggplot(aes(x=f,y=ph,col=like)) +geom_point(shape=1,size=0.2)+
  geom_point(data = Max2,aes(x=f,y=ph),size=2,col="red")+
  geom_vline(xintercept = frq.mx[1],col=2)+
  labs(x= "frequency [1/month]", y= "phase [radians]",
       title = "Most Likeli Harmonic Parameters",
       subtitle = "Model Function: cos(2pi*f*t+ph)")
#==============
require(bbmle)
# Example negative log-likelihood function
harmonic_nll <- function(frequency, phase, amplitude = 1) {
  predicted <- amplitude * sin(2 * pi * frequency * time + phase)
  residuals <- observed - predicted
  -sum(dnorm(residuals, mean = 0, sd = sd(sin(time)), log = TRUE))  # Gaussian likelihood
}
# start values
library(bbmle)
#  2 dominant frequencies by fourier
FFT=fft(NINO.res$lin.res)
amp.mx=FFT[1:(N%/%2)]%>%Mod()%>%sort(decreasing = T)%>%head(2) #[1] 216.0108 203.4824
idx.mx=which(Mod(FFT[1:(N%/%2)])%in% amp.mx) # 28 42
frq.mx=idx.mx-1 # 27  41



# Define starting values and bounds
start_vals <- list(frequency = frq.mx[1]/N, phase = -1)
lower_bounds <- c(frequency = frq.mx[1]*0.9/N, phase = -2.2)  # Lower bounds
upper_bounds <- c(frequency = frq.mx[1]*1.1/N, phase = 1)   # Upper bounds

# Fit the model with bounded optimization
fit <- mle2(
  minuslogl = harmonic_nll,
  start = start_vals,
  method = "L-BFGS-B",  # Optimizer that supports bounds
  lower = lower_bounds,
  upper = upper_bounds,
  data = list(time = time, observed = observed, sigma = sd(sin(time)))
)
Max2=bind_rows(Likes[13440,],L2[13638,])
opt=coef(fit)
opt.prms=c(opt[1],opt[2],log.lik(signl,f=opt[1],phi=opt[2]))
names(opt.prms)<-c("f","ph","like")
Mx.prms=bind_rows(Max2,opt.prms)
log.lik(signl = NINO.res$lin.res,f=frq.mx[1]/N,phi = 0)
opt=coef(fit)
Likes%>%
  ggplot(aes(x=f,y=ph,col=like)) +geom_point(shape=1,size=0.2)+
  geom_point(data = Max2,aes(x=f,y=ph),size=2,col="red")+
  geom_vline(xintercept = frq.mx[1],col=2)+
  labs(x= "frequency [1/month]", y= "phase [radians]",
       title = "Most Likeli Harmonic Parameters",
       subtitle = "Model Function: sin(2pi*f*t+ph)")
