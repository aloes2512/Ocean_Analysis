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
                    na = "-1") |>
  mutate(date = ymd(paste(year, month, day))) # safe better than converting dec_date
# proportion 0 counts
(which(sn_data$SN<1)>0) %>%sum()/length(sn_data$SN) # 15%
is.na(sn_data$SN) %>%sum() #3247
(is.na(sn_data$SN) %>%sum())/length(sn_data$SN) # 5%
range(sn_data$date) # "1818-01-01" "2026-04-30"
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
SP_daily%>%NROW() # 61446 168.234 years 2018.815 month
colnames(SP_daily)
# set NA to zero
SP_no.obs=SP_daily%>% mutate(SP=ifelse(is.na(SN), 0,SN))
(which(sn_data$SN<1)>0) %>%sum() # 11398
NROW(SP_no.obs) # 61446 168 years 2018.815 month
FFT=tibble(idx=1:NROW(SP_no.obs),spc=fft(SP_no.obs$SN),amp=Mod(spc))
FFT%>%subset(idx<NROW(SP_no.obs)/2)%>%arrange(desc(amp))
# interpolate by splines to get regular time series
SP_intpol.spl=spline(x=SP_daily$date,y=SP_daily$SN,
                     xout = seq.Date(from = first(SP_daily$date),
                                      to= last(SP_daily$date),by= "1 days"))
length(SP_intpol.spl$x) # !!spline cuts off the leading NA!! 76084
length(SP_intpol.spl$y) # 76084
# format interpolated data into  tibble
SP_intpol= tibble(date=SP_intpol.spl$x,
                  SP.intpl=SP_intpol.spl$y)
#spline interpolation produced " SP  negative counts" which needs to be eliminated
SP_intpol=SP_intpol%>%mutate(SP.intpl= if_else(SP.intpl<0, 0, SP.intpl))
colnames(SP_intpol) #date SP.intpl
SP_intpol%>%ggplot(aes(x=date,y=SP.intpl))+
  geom_point(size=0.2,col="grey")+
  geom_smooth()
require(TSA)
SP_prd=tibble(freq=1:length(SP_intpol$SP.intpl)-1,
              spc= fft(SP_intpol$SP.intpl),
              amp=Mod(spc))
range(SP_daily$date)#  "1818-01-08" "2026-04-30"
n=NROW(SP_prd) # 76084 days
dt_SP_periods <- tibble(freq=SP_prd$freq,amp=SP_prd$amp)%>%
  arrange(desc(amp))%>%
  mutate(Period.d= n/freq,Period.yr=round(Period.d/365.24,2))
dt_SP_periods%>%subset(freq<n/2& freq>0)

#================
range(SP_daily$date) # "1818-01-02" "2024-12-31"
SP_intpol=SP_intpol%>%mutate(nmdt=1:NROW(SP_intpol))
#  fit with neg.binomial
library(MASS)
negative.binomial(theta = 2,link="log")
library(mgcv)
SP_mdl.nbin <- mgcv::gam(data=SP_intpol,formula = SP.intpl~s(nmdt,k=60),family=negbin(1),scale=1)
summary(SP_mdl.nbin) # intercept 4.170944

AIC(SP_mdl.nbin)# [1] 789054.4
require(broom)
SP_NB_fit<-SP_mdl.nbin%>%augment()%>%
  dplyr::select(nmdt,SP.intpl,.fitted)%>%
  mutate(date=seq.Date(as.Date("1818-01-02"),as.Date("2026-04-30"), length.out = 76084),
         SP_fit=exp(.fitted))
head(SP_NB_fit)


summary(SP_NB_fit)
saveRDS(SP_NB_fit,"data/SP_NB_fit.rds")
rm(SP_NB_fit)
SP_NB_fit=readRDS("data/SP_NB_fit.rds")
head(SP_NB_fit)
SP_NB_plt<-SP_NB_fit%>%
  ggplot(aes(x=date,y= SP_fit))+
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
  subset(idx< NROW(SP_NB_fit)/2)%>%
  arrange(desc(amp))
N=NROW(SP_NB_fit) # 76084 days
(N/19)/365.24 # 10.96 years
round((N/c(19,20,21,18,26,17))/365.24,2) # 10.96 10.42  9.92 11.57 8.01 12.25
library(dplyr)
library(lubridate)
library(MASS)

# 1. Add numeric dates
SP_intpol <- SP_intpol %>%
  mutate(date = SP_NB_fit$date)  # Remove negatives first
colnames(SP_intpol) #"date" "SN_interpolated" "nmdt"
# 2. Create complete daily response variable (0s for missing days)

full_data <- SP_intpol%>%rename("SN_count"=SP.intpl)
full_data=full_data%>%mutate(SN=round(SN_count)%>%as.integer())
library(MASS)
# 3. Fit negative binomial with smooth time trend
nb_rigid_formula <- as.formula(paste("SN ~ nmdt "))

fit_nb_rigid <- glm.nb(nb_rigid_formula, data = full_data)

augment(fit_nb_rigid )%>%pull(.fitted)%>%range()
# 4. Predict fitted values (non-negative guaranteed)
full_data$SN_fitted <- predict(fit_nb_rigid, type = "response")

# Verify
summary(fit_nb_rigid)
plot(full_data$date, full_data$SN_fitted, type = "l",
     main = "Negative Binomial Smoothed Sunspots")







#============================
require(xts)
SP_NB.fitted.xts<-xts(SP_NB_fit$SP_fit,order.by = SP_NB_fit$date)
SP_period.tb<-TSA::periodogram(y=coredata(SP_NB.fitted.xts))
SP_NB_periods<-tibble(Freq=SP_period.tb$freq,
                        Spec=SP_period.tb$spec)
dim(SP_NB_fit)
SP_NB_periods.<-SP_NB_periods%>%subset(Spec>11000)
head(SP_NB_periods.,15)
max.Freq.Spec<-SP_NB_periods.[which.max(SP_NB_periods.$Spec),]#19
(1/max.Freq.Spec$Freq)/365.24# 11.06698 yrs
SP_periods<-SP_NB_periods.%>% mutate(Period.yr=(1/Freq)/364.5)%>%arrange(desc(Spec))
SP_NB_periods.%>%ggplot(aes(x=Freq,y=Spec))+
  geom_line()

SP_NB_fit
# relation of sunspots to Total Solar Irradiance (TSI) Ambelu et al.
browseURL("http://www.lajpe.org/dec11/LAJPE_576_Ambelu_Tebabal_prreprint_corr.pdf")
a= 1365.66 # W/m^2
b= 0.0136  # W/(SP*m^2)
c=-0.0000076 # W/(SP*m)^2
#TSI= a+ b*SP+c*SP^2
S_power=SP_NB_fit%>% mutate(TSI=a+b*SP_fit+c*SP_fit^2,TSI=1000*TSI/a)%>%dplyr::select(date,nmdt,TSI)
saveRDS(S_power,"data/S_power.rds")
S_power%>%ggplot(aes(x=date,y=TSI))+geom_line()+
  labs(y="Solar Peak Power [W/m^2]", title = "Solar Power Variation",
       subtitle = "@ sea level,estimate from sunspot smoothed count")
 ggsave("figs/Solar_Irradiation_equator.png")
