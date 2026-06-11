# sunspots periodicity
require(tidyverse)
SP_path<-"http://www.sidc.be/silso"
browseURL(SP_path)
# download from silso
#`=====`
#Contents:
#Column 1: Gregorian Year
#Column 2: Gregorian Month
#Column 3: Gregorian Day
#Column 4: Decimal date
#Column 5: Estimated Sunspot Number
#Column 6: Estimated Standard Deviation
#Column 7: Number of Stations calculated
#Column 8: Number of Stations available
#v.names<-c("Jahr","Monat","Tag","dat.dec","SP","SD","N.calc","N.station")
Source="WDC-SILSO, Royal Observatory of Belgium, Brussels"
browseURL("https://www.sidc.be/SILSO/datafiles")
# Install if needed: install.packages("tidyverse")
sn_data <- read_csv2("https://www.sidc.be/SILSO/DATA/SN_d_tot_V2.0.csv",
                    col_names = c("year", "month", "day", "dec_date", "SN", "SN_err", "N_obs"),
                    na = "-1")  # safe better than converting dec_date

sn_data=sn_data%>% mutate(month=as.integer(month),
                  day=as.integer(day),
                  SN_err=as.double(SN_err))%>%
            mutate(dt.mnth=year+(month-1)/12,dec.yr=dt.mnth+day/364.25)%>%
            dplyr::select(dt.mnth,dec.yr,month,SN)
# proportion 0 counts
(which(sn_data$SN<1)>0) %>%sum()/NROW(sn_data) # 15%
is.na(sn_data$SN) %>%sum() #3247
(is.na(sn_data$SN) %>%sum())/NROW(sn_data) # 5%
range(sn_data$dt.mnth) # 1818.000 2026.333"
summary(sn_data) # NA 3247
NROW(sn_data) # 76091
#SP_daily_raw<-read.csv2(sn.file.path,col.names = v.names,na.strings = c("-1",NA))
#SP_daily_raw=SP_daily_raw%>%mutate(date = glue::glue("{Jahr}-{Monat}-{Tag}"),
#                                   date = as.Date(date, format = "%Y-%m-%d"))%>%
#  dplyr::select(date,SP)
#SP_dailyobs=SP_daily_raw%>% mutate(SP = if_else(SP == -1, NA_real_, SP))

#summary(SP_dailyobs) # 3246 NA
#NROW(SP_dailyobs) # 75605
#'==='
#'select only observations with more than zero
#SP_daily_1min=SP_dailyobs%>% subset(SP>0 & !is.na(SP))
#SP_daily_min= sn_data%>%dplyr::select(date,SN)%>% subset(SN>0 & !is.na(SN))
#======================
SP_daily=sn_data%>%subset(SN>0 & !is.na(SN)) # first 7 SN in 1801-01 are NA
SP_daily%>%NROW() # 61477 168.234 years 2018.815 month
colnames(SP_daily) #dt.mnth dec.yr month SN
# set NA to zero
SP_no.obs=SP_daily%>% mutate(SP=ifelse(is.na(SN), 0,SN))
(which(sn_data$SN<1)>0) %>%sum() # 11398
NROW(SP_no.obs) # 61477 168.3195 years 2019.834 month
#=========

#==========
FFT=tibble(idx=1:NROW(SP_no.obs)-1,spc=fft(SP_no.obs$SN),amp=Mod(spc))
FFT%>%subset(idx<NROW(SP_no.obs)/2)%>%arrange(desc(amp))
# interpolate by splines to get regular time series
#------------
SP_intpol.spl=spline(x=SP_daily$dec.yr,y=SP_daily$SN,
                     xout = seq(from = first(SP_daily$dec.yr),
                                to= last(SP_daily$dec.yr),by= 1/365.24))
length(SP_intpol.spl$x) # !!spline cuts off the leading NA!! 76115
length(SP_intpol.spl$y) # 76115
SP_intpol= tibble(dec.yr=SP_intpol.spl$x,
                  SP.intpl=SP_intpol.spl$y)

#-------------

#spline interpolation produced " SP  negative counts" which needs to be eliminated
SP_intpol=SP_intpol%>%mutate(SP.intpl= if_else(SP.intpl<0, 0, SP.intpl))
colnames(SP_intpol) #dec.yr  SP.intpl
SP_intpol%>%ggplot(aes(x=dec.yr,y=SP.intpl))+
  geom_point(size=0.2,col="grey")+
  geom_smooth()

require(TSA)
#___
FFT.SP=tibble(freq=1:length(SP_intpol$SP.intpl)-1,
              spc= fft(SP_intpol$SP.intpl),
              amp=Mod(spc))
n=length(SP_intpol$SP.intpl) # 76115
FFT.SP%>%arrange(desc(amp))%>%subset(freq<n/2)%>%
  mutate(Period.d= n/freq,Period.yr=round(Period.d/365.24,2))
range(SP_daily$dec.yr)#  "1818.022" "2026.418"


#================

#  fit with neg.binomial
library(MASS)
negative.binomial(theta = 2,link="log")
library(mgcv)
SP_mdl.nbin <- mgcv::gam(data=SP_intpol,formula = SP.intpl~s(dec.yr,k=60),family=negbin(1),scale=1)
summary(SP_mdl.nbin) # intercept 4.176396

AIC(SP_mdl.nbin)# [1] 789681.4
require(broom)
SP_NB_fit<-SP_mdl.nbin%>%augment()%>%
  dplyr::select(dec.yr,SP.intpl,.fitted)%>%
  mutate(SP_fit=exp(.fitted))
head(SP_NB_fit)


summary(SP_NB_fit)
saveRDS(SP_NB_fit,"data/SP_NB_fit.rds")
rm(SP_NB_fit)
SP_NB_fit=readRDS("data/SP_NB_fit.rds")
head(SP_NB_fit)
SP_NB_plt<-SP_NB_fit%>%
  ggplot(aes(x=dec.yr,y= SP_fit))+
  geom_line()+
  ggtitle("Sunspots daily count",
          subtitle = "fitted(red) with NegBinom")+
  labs(x="",y= "SP counts",caption = "data: www.sidc.be/silso")
SP_NB_plt+geom_smooth(aes(y= SP_fit),col="red")

# ===================================
FFT.SP=tibble(idx=1:NROW(SP_NB_fit),
              spc=SP_NB_fit$SP_fit%>%fft(),
              amp=Mod(spc))
FFT.SP%>%
  subset(freq< NROW(SP_NB_fit)/2)%>%
  arrange(desc(amp))
N=NROW(SP_NB_fit) # 76115 days
(N/19)/365.24 # 10.96 years
round((N/c(19,20,21,18,26,17))/365.24,2) # 10.97 10.42  9.92 11.58 8.02 12.26
library(dplyr)
library(lubridate)
library(MASS)

# 1. Add numeric dates
SP_NB_fit  <- SP_NB_fit %>%
  mutate(date = date_decimal(SP_NB_fit$dec.yr))  # Remove negatives first
colnames(SP_NB_fit ) #"date" "SN_interpolated" "nmdt"
#============================
require(xts)
SP_NB.fitted.xts<-xts(SP_NB_fit$SP_fit,order.by = SP_NB_fit$date)
SP_period.tb<-TSA::periodogram(y=coredata(SP_NB.fitted.xts))
SP_NB_periods<-tibble(Freq=SP_period.tb$freq,
                        Spec=SP_period.tb$spec)
dim(SP_NB_fit) #76115   5
SP_NB_periods.<-SP_NB_periods%>%subset(Spec>11000)
head(SP_NB_periods.,15)
max.Freq.Spec<-SP_NB_periods.[which.max(SP_NB_periods.$Spec),]#19
(1/max.Freq.Spec$Freq)/365.24# 11.06698 yrs
SP_periods<-SP_NB_periods.%>% mutate(Period.yr=(1/Freq)/364.5)%>%arrange(desc(Spec))
SP_NB_periods.%>%ggplot(aes(x=Freq,y=Spec))+
  geom_point()

SP_NB_fit
# relation of sunspots to Total Solar Irradiance (TSI) Ambelu et al.
browseURL("http://www.lajpe.org/dec11/LAJPE_576_Ambelu_Tebabal_prreprint_corr.pdf")
a= 1365.66 # W/m^2
b= 0.0136  # W/(SP*m^2)
c=-0.0000076 # W/(SP*m)^2
#TSI= a+ b*SP+c*SP^2
S_power=SP_NB_fit%>% mutate(TSI=a+b*SP_fit+c*SP_fit^2,TSI=1000*TSI/a)%>%dplyr::select(date,dec.yr,TSI)
saveRDS(S_power,"data/S_power.rds")
S_power%>%ggplot(aes(x=date,y=TSI))+geom_line()+
  labs(y="Solar Peak Power [W/m^2]", title = "Solar Power Variation",
       subtitle = "@ sea level,estimate from sunspot smoothed count",
       caption = "http://www.lajpe.org/dec11/LAJPE_576_Ambelu_Tebabal..")
 ggsave("figs/Solar_Irradiation_equator.png")


