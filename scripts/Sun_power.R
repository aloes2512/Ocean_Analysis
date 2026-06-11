# sunspots periodicity
require(tidyverse)
SP_path<-"http://www.sidc.be/silso"
#browseURL(SP_path)
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
#=========
sn_data <- read_csv2("https://www.sidc.be/SILSO/DATA/SN_d_tot_V2.0.csv",
                     col_names = c("year", "month", "day", "dec_date", "SN", "SN_err", "N_obs"),
                     na = "-1")  # safe better than converting dec_date

sn_data=sn_data%>% mutate(month=as.integer(month),
                          day=as.integer(day),
                          SN_err=as.double(SN_err))%>%
  mutate(dt.mnth=year+(month-1)/12,dec.yr=dt.mnth+day/364.25)%>%
  dplyr::select(dt.mnth,dec.yr,month,SN)
#======================
SP_daily=sn_data%>%subset(SN>0 & !is.na(SN)) # first 7 SN in 1801-01 are NA
SP_intpol.spl=spline(x=SP_daily$dec.yr,y=SP_daily$SN,
                     xout = seq(from = first(SP_daily$dec.yr),
                                to= last(SP_daily$dec.yr),by= 1/365.24))
length(SP_intpol.spl$x) # !!spline cuts off the leading NA!! 76115
length(SP_intpol.spl$y) # 76115
SP_intpol= tibble(dec.yr=SP_intpol.spl$x,
                  SP.intpl=SP_intpol.spl$y)

SP_intpol=SP_intpol%>%mutate(SP.intpl= if_else(SP.intpl<0, 0, SP.intpl))
library(MASS)
negative.binomial(theta = 2,link="log")
library(mgcv)
SP_mdl.nbin <- mgcv::gam(data=SP_intpol,formula = SP.intpl~s(dec.yr,k=60),family=negbin(1),scale=1)
require(broom)
SP_NB_fit<-SP_mdl.nbin%>%augment()%>%
  dplyr::select(dec.yr,SP.intpl,.fitted)%>%
  mutate(SP_fit=exp(.fitted))



S_power.new=SP_NB_fit
colnames(S_power.new)
NROW(S_power.new) # 76115 daily values
num.yrmnth=function(dec.yr){
  yr=floor(dec.yr)
  mnth=floor((dec.yr-yr)*12)
  dt.mnth=yr+mnth/12
  return(dt.mnth)
}
S_power.new=S_power.new%>%mutate(dt.mnth=num.yrmnth(dec.yr))
SP_mnthly=S_power.new%>%group_by(dt.mnth)%>%summarise(SP_mnth=mean(SP_fit))
# relation of sunspots to Total Solar Irradiance (TSI) Ambelu et al.
#browseURL("http://www.lajpe.org/dec11/LAJPE_576_Ambelu_Tebabal_prreprint_corr.pdf")
a= 1365.66 # W/m^2
b= 0.0136  # W/(SP*m^2)
c=-0.0000076 # W/(SP*m)^2
#TSI= a+ b*SP+c*SP^2
S_power.mnthly=SP_mnthly%>%mutate(TSI=a+b*SP_mnth+c*SP_mnth^2,TSI=1000*TSI/a)%>%dplyr::select(dt.mnth,TSI)
S_power.mnthly%>%ggplot(aes(x=dt.mnth,y=TSI))+geom_line()+
  labs(y="Solar Peak Power [W/m^2]", title = "Solar Power Variation",
       subtitle = "@ sea level,estimate from sunspot smoothed count",
       caption = "http://www.lajpe.org/dec11/LAJPE_576_Ambelu_Tebabal..")
ggsave("figs/Solar_Irradiation_equator.png")
saveRDS(S_power.mnthly,"data/S_power.rds")
