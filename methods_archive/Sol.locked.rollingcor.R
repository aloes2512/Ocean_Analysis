# copied from Solar_locked_anomaly.R
library(tidyverse)
library(itsmr)
Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
my_df<-Ocean_solar_anomaly
mytrend.sol<-my_df$trend.solar
my_df<-my_df%>%mutate(sol.trd7=trend(trend.solar,7),
                      dif.sol.trd=trend.solar-sol.trd7)

dt.mnth=my_df$dt.mnth
sol.trd.secular=loess(mytrend.sol~ dt.mnth,span = 0.5) %>% predict()
# separate seasonal and semiseasonal
# trend.solar from locking Anomaly to solar power
M=c("season",12,"season",6)
my_DF=tibble(dt.mnth=dt.mnth,sec.trd=sol.trd.secular,
             mytrend.sol=my_df$trend.solar,
             anomaly=Resid(my_df$Anomaly,M))
my_DF=my_DF%>%mutate(difer=mytrend.sol-sol.trd.secular,
                     sol.residuals=anomaly-sol.trd.secular,
                     sol.res.hr= hr(sol.residuals,N/1:20))
N=NROW(my_df)
my_DF%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=sol.residuals))+
  geom_line(aes(y=sol.res.hr,col="sol.res.hr"))

vector_obs<-hr(my_DF$anomaly,N/1:20)
vector_preind<-my_DF$sol.res.hr

rolling_corr.sol <- slide_index_dbl(
  .x = 1:length(vector_obs),
  .i = 1:length(vector_obs),
  .f = function(idx) {
    # If the window doesn't have a full 11 years of data, return NA
    if(length(idx) < 132) return(NA)

    # Calculate simple vector correlation on the current slice of indices
    cor(vector_obs[idx], vector_preind[idx])
  },
  .before = 132,   # Looks 131 steps backward + current step = 132 months
  .complete = TRUE # Automatically returns NA for incomplete windows at the start
)
RollingCor=tibble(Time=Ocean_solar_anom$dt.mnth,rolling_corr.sol=as.numeric(rolling_corr.sol))#
RollingCor%>%ggplot(aes(x=Time))+
  geom_line(aes(y=rolling_corr.sol))+
  geom_vline(xintercept = 1997,linetype = 2)+
  geom_vline(xintercept = 2008,linetype = 2)+
  geom_vline(xintercept = 1977,linetype = 3)+
  geom_vline(xintercept = 1988,linetype = 3)+
  geom_vline(xintercept = 2015,linetype = 4)+
  geom_vline(xintercept = 2026,linetype = 4)
#A: 1977 to 1988 loss of correlation related to
# PDO (Pacific Decadal Oscillation)
#`1976–1977, the Pacific Ocean underwent a sudden,
# well-documented regime shift from a prolonged negative (cool)
#PDO phase to a strong positive (warm) phase.`
#B. 2nd break  1997/1998 to 2008 Shift:
#`The record-breaking 1997/1998 El Niño event
#triggered another major PDO phase reorganization,
#coinciding with accelerated global ocean warming.
#This lines up directly with the sharp 1998 drop
#(centered in your 1997–2008 rolling window)
#where phase coherence completely collapsed.`
#C. `Post-2010 Oscillations:
## The quick recovery back to $+0.95$ around 2018 (window 2007–2018)
##followed by another sharp crash toward 2026
##indicates that the observed lowpass variability is
##no longer in stable phase alignment
##with the fixed pre-1938 solar periodicities.
#It is beating against MODERN FORCED SECULAR TREND





