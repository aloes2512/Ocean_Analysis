my.source <- readRDS("~/projects/NOAA TS Analysis/NOAA.psl.rds")
my.source$source
URL="https://psl.noaa.gov/data/timeseries/month/"
browseURL(URL)
# 2 timeseries available starting 1870 or 1950
nino.long=readRDS("~/projects/Global_Temp/data/NINO3.4.rds")
my.source$data%>%summary()
library(tidyverse)
my.source$data%>%map_dbl(~ NROW(.x)) # range 564 to 2064
nino.short=readRDS("~/projects/Global_Temp/data/NOAA.psl.monthly.rds")%>%
  {.$NINO34}
library(itsmr)
library(zoo)
M= c("trend",1,"season",12)
nino.short=nino.short%>%
  mutate(res= Resid(value,M),
         res.s=smooth.ma(res,12),
         yr.mnth=yearmon(year+(month-1)/12))
dat.shrt=nino.short%>%dplyr::select(yr.mnth,res.s)
dat.shrt%>%ggplot(aes(x=yr.mnth,y=res.s))+geom_line()
dat.lng=nino.long%>%mutate(res=Resid(tmp_C,M),res.l=smooth.ma(res,12))
dat.lng=dat.lng%>%dplyr::select(yr.mnth,res.l)
dat.lng%>%ggplot(aes(x=yr.mnth,y=res.l)) +geom_line()
dat.lng%>%left_join(dat.shrt)%>%
  ggplot(aes(x=yr.mnth))+
  geom_line(aes(y=res.l),col=2,linewidth = 2)+geom_line(aes(y=res.s),col=3)
FFT.shrt=tibble(idx=seq_along(dat.shrt$yr.mnth),
                          spc.s=fft(dat.shrt$res.s),
                          amp=Mod(spc.s))
FFT.lng=tibble(idx=seq_along(dat.lng$yr.mnth),
                spc.l=fft(dat.lng$res.l),
                amp=Mod(spc.l))
FFT.shrt%>%arrange(desc(amp))%>%subset(idx< 100)
FFT.lng%>%arrange(desc(amp))%>%subset(idx< 100)
# max harmonics look very different
#try length of longer integer multible of shorter
del=NROW(dat.lng)-2*NROW(dat.shrt) #29
dat.lng.new=dat.lng[-c(1:del),]
FFT.lng_new=tibble(idx=seq_along(dat.lng.new$yr.mnth),
               spc.l=fft(dat.lng.new$res.l),
               amp=Mod(spc.l))
FFT.lng_new%>%arrange(desc(amp))
# select hr Nl/1.20
Nl=NROW(dat.lng.new) # 1824
Ns=NROW(dat.shrt) # 912
Ns/(1:10)==
Nl/(2*(1:10))
lng.har=dat.lng.new%>%mutate(lng.har=hr(res.l,Nl/(2*(1:10))))
colnames(lng.har)
lngplt=lng.har%>%ggplot(aes(x=yr.mnth,y=lng.har))+
  geom_line()+labs(x="",y="sum amplitude",title = "Sum 10 Harm longest periods")
shrt.har=dat.shrt%>%mutate(shrt.har=hr(res.s,Nl/(2*(1:10))))
shrt.har%>%ggplot(aes(x=yr.mnth,y=shrt.har))+geom_line()
lngplt+geom_line(data=shrt.har,aes(x=yr.mnth,y=shrt.har),col=2)
Nl/(2*(1:10))
Ns/1:10
#  restrict longer series to 2 *shorter
dat.lng.new=dat.lng[-c(1:del),]
FFT.new= tibble(idx=1:NROW(dat.lng.new),
                spc.l= fft(dat.lng.new$res.l),
                amp=Mod(spc.l))
FFT.new%>%arrange(desc(amp))

