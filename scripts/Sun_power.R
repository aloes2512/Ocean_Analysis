# Wikipedia "Radiative Forcing"
# sunspots
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
Source<-"WDC-SILSO, Royal Observatory of Belgium, Brussels"
#=========
library(tidyverse)

sn_data.new <- read_csv2(
  "https://www.sidc.be/SILSO/DATA/SN_d_tot_V2.0.csv",
  col_names = c("year", "month", "day", "dec_date", "SN", "SN_err", "N_obs", "Nr_stations"),
  col_types = "iiidddii",
  na = c("-1", "-1.0")
) %>%
  mutate(
    dt.mnth = year + (month - 1) / 12,
    dec.yr  = dt.mnth + day / 365.25
  ) %>%
  dplyr::select(dt.mnth, dec.yr, month, SN)
#============
SN_m_tot_V2.0.txt
# Column 1-2: Gregorian calendar date
#- Year
#- Month
#Column 3: Date in fraction of year for the middle of the corresponding month
#Column 4: Monthly mean total sunspot number.
#Column 5: Monthly mean standard deviation of the input sunspot numbers from individual stations.
#Column 6: Number of observations used to compute the monthly mean total sunspot number.
#Column 7: Definitive/provisional marker. A blank indicates that the value is definitive. A '*' symbol indicates that the monthly value is still provisional and is subject to a possible revision (Usually the last 3 to 6 months)

#Line format [character position]:
 # - [1-4] Year
#- [6-7] Month
#- [9-16] Decimal date
#- [19-23] Monthly total sunspot number
#- [25-29] Standard deviation
#- [32-35] Number of observations
#- [37] Definitive/provisional indicator
sp_nms <- c("Year", "mth", "dec.dt", "SP_N", "SP_sd", "N.obs", "prov_Fin")

sp_nms <- c("Year", "mth", "dec.dt", "SP_N", "SP_sd", "N.obs", "prov_Fin")

sn_month <- read.table(
  file = file.path(Sys.getenv("HOME"), "Downloads", "SN_m_tot_V2.0.csv"),
  header = FALSE,
  sep = ";",          # The file is using semicolons as delimiters here
  col.names = sp_nms,
  comment.char = "#",
  fill = TRUE,
  stringsAsFactors = FALSE
)
library(zoo)
# Clean up any trailing empty columns or semicolons if necessary
sn_month <- sn_month %>%
            dplyr::select(1,2,4)%>%mutate(dt.mnth=Year+(mth-1)/12)
#======================
SP_daily=sn_data.new%>%subset(SN>0 & !is.na(SN)) # first 7 SN in 1801-01 are NA
SP_intpol.spl=spline(x=SP_daily$dec.yr,y=SP_daily$SN,
                     xout = seq(from = first(SP_daily$dec.yr),
                                to= last(SP_daily$dec.yr),by= 1/365.24))
SP_intpol= tibble(dec.yr=SP_intpol.spl$x,
                  SP.intpl=SP_intpol.spl$y) # NROW 76236

SP_intpol=SP_intpol%>%mutate(SP.intpl= if_else(SP.intpl<0, 0, SP.intpl))
library(MASS)
negative.binomial(theta = 2,link="log")
library(mgcv)
SP_mdl.nbin <- mgcv::gam(data=SP_intpol,formula = SP.intpl~s(dec.yr,k=60),family=negbin(1),scale=1)
require(broom)
S_power=tibble(dec.yr=SP_intpol.spl$x,
                 SP_fit=exp(predict(SP_mdl.nbin)))
dim(S_power) #76236
library(zoo)
S_power<-S_power%>% mutate(dt.mnth=as.yearmon(dec.yr),dt.mnth=as.numeric(dt.mnth))

SP_mnthly=S_power%>%group_by(dt.mnth)%>%summarise(SP_mnth=mean(SP_fit))
NROW(SP_mnthly)
# relation of sunspots to Total Solar Irradiance (TSI) Ambelu et al.
#browseURL("http://www.lajpe.org/dec11/LAJPE_576_Ambelu_Tebabal_prreprint_corr.pdf")
a= 1365.66 # W/m^2
b= 0.0136  # W/(SP*m^2)
c=-0.0000076 # W/(SP*m)^2
#TSI= a+ b*SP+c*SP^2
S_power.mnthly=SP_mnthly%>%mutate(TSI=a+b*SP_mnth+c*SP_mnth^2,TSI=1000*TSI/a)%>%dplyr::select(dt.mnth,TSI)
plt.sol.power=S_power.mnthly%>%ggplot(aes(x=dt.mnth,y=TSI))+geom_line()+
  labs(y="Solar Peak Power [W/m^2]", title = "Solar Power Variation",
       subtitle = "@ sea level,estimate from sunspot smoothed count",
       caption = "http://www.lajpe.org/dec11/LAJPE_576_Ambelu_Tebabal..")
saveRDS(S_power.mnthly,"data/S_power.rds")
#======
S_power_ext=sn_month%>%mutate(TSI=a+b*SP_N+c*SP_N^2,TSI=1000*TSI/a)%>%
  dplyr::select(dt.mnth,TSI)%>%
  mutate(S.irr=TSI-mean(TSI,na.rm=T))
S_power_variation<-S_power.mnthly%>%mutate(SI=TSI-mean(TSI,na.rm=T))
S_power_variation%>%ggplot(aes(x=dt.mnth,y=SI))+geom_line()+
  labs(y="Solar Peak Power [W/m^2]", title = "Solar Power Variation",
       subtitle = "@ sea level,estimate from sunspot smoothed count",
       caption = "http://www.lajpe.org/dec11/LAJPE_576_Ambelu_Tebabal..")
library(ggplot2)
S_power_ext%>%ggplot(aes(x=dt.mnth))+geom_line (aes(y=S.irr))
S_power_variation %>%
  ggplot(aes(x = dt.mnth, y = SI)) +
  # Originale Solar-Leistung (Schwabe- & Hale-Zyklen)
  geom_line(color = "grey30", alpha = 0.8, linewidth = 0.4) +
  # Gleissberg-Zyklus als langsamer Trend (Einhüllende mittels LOESS)
  geom_smooth(method = "loess", formula = y ~ x, span = 0.25,
              color = "firebrick", linewidth = 1.1, se = FALSE) +
  labs(
    y = "Solar Peak Power [W/m^2]",
    title = "Solar Power Variation ( Gleissberg-Trend)",
    subtitle = "@ sea level, estimate from sunspot smoothed count",
    caption = "http://www.lajpe.org/dec11/LAJPE_576_Ambelu_Tebabal.."
  ) +
  theme_minimal()
ggsave("figs/sol.power.png")
# see Wikipedia on CO2 effect (Absorption) to be logarithmic
#browseURL("https://de.wikipedia.org/wiki/Klimasensitivit%C3%A4t")
#browseURL("https://gml.noaa.gov/aggi/aggi.html")

