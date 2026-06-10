require(lubridate)
t=seq.Date(from = ymd("1870-01-01"), to= ymd("2024-12-31"),by="1 days")
tnum= (1:length(t))-1
n=length(t)
mdl.fun= function(f,ph,n){
  t= 0:(n-1)
  f0=1/n
  fun= sin(2*pi*f0*f*t+ph)
  fun=fun/sd(fun)
  return(fun)
}
Y=tibble(tnum=tnum,
         y=mdl.fun(f=1,ph=pi,n=n))
summary(Y)
Y%>%ggplot(aes(x=tnum,y=y))+geom_line()
Y.sub=subset(Y,tnum<n/2)
NROW(Y.sub)# 28307
Y.sub%>%ggplot(aes(x=tnum,y=y))+geom_line()
fft(Y.sub$y)%>%NROW() # 28307
which.max(Mod(fft(Y.sub$y))) # 1
require(gsignal)
y.rslt=ifft(fft(Y.sub$y))
length(y.rslt)#
tibble(t=(1:length(y.rslt))-1,y.rslt)%>%
  ggplot(aes(x=t,y=y.rslt))+geom_line()       
# check ampl of fft
FFT.sub=tibble(idx=1:NROW(Y.sub),spc=fft(Y.sub$y),amp=Mod(spc))
FFT.sub%>%subset(idx<10)%>%
  ggplot(aes(x=idx,y=amp))+geom_point()
FFT.sub%>%subset(amp> as.numeric(FFT.sub[6,3]))
new.FFT=tibble(idx=1:NROW(FFT.sub),spc=0i,amp=Mod(spc))
spc.max=FFT.sub%>%subset(amp> as.numeric(FFT.sub[6,3]))%>%pull(idx)
new.FFT[spc.max,]<-FFT.sub[spc.max,]
y.recstr=ifft(new.FFT$spc)

tibble(t=(1:length(y.recstr))-1,y.recstr)%>%
  ggplot(aes(x=t,y=y.recstr))+geom_line()       
#=================
#original data
Y=tibble(tnum=tnum,
         y=mdl.fun(f=1,ph=pi,n=n))
tibble(
y.long=c(y.recstr,-y.recstr),
t=(1:length(y.long))-1)%>%ggplot(aes(x=t,y=y.long))+
  geom_line(col="red")+
  geom_line(data=Y,aes(x=tnum,y))
# with added noise
