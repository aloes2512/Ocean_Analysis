# Schwabe 132 month
n_shrt=1058
n_shrt/132 # 8.015152
f_schwabe=132/n_shrt # 0.125; period 132 month
f.0=1/n_shrt # 0.000925
f.ny=f.0*n_shrt/2 # 0.5
f_cutoff<-f_schwabe/f.ny # 0.244
#-----
ts.list=readRDS("data/Area10+.ts.rds")
global.ts=ts.list$global.ts
global.ts$dt.mnth[1058] #1938.083
glob.shrt.ts <- global.ts %>% subset(dt.mnth <= 1938.083+1/12)
NROW(glob.shrt.ts )
# season 6,12,remove
M=c("season",12,"season",6)
glob.1938=glob.shrt.ts %>% mutate(deca.anoma=Resid(anomaly,M))
glob.1938%>%ggplot(aes(x=dt.mnth))+geom_line(aes(y=deca.anoma))
FFT.1938=tibble(idx=1:n_shrt,
                spc=fft(glob.1938$deca.anoma-mean(glob.1938$deca.anoma)),
                amp=Mod(spc))
FFT.max=FFT.1938%>%arrange(desc(amp))%>%subset(idx< 540&amp>10)

range(FFT.1938$amp) # 3.330669e-15 7.564163e+01
idx.max=FFT.max$idx # 2 89  5 16  7  8 24  6 11 90 15 29 34 12 20 27 22 14 28 31
amp.min=range(FFT.max$amp)[1] #   8.85
FFT_filt <- FFT.1938|>
  mutate(
    spc = if_else(amp <10, 0 + 0i, spc),
    amp = Mod(spc)
  )
which(FFT_filt$amp>=10)
dim(FFT_filt)
colnames(FFT_filt)
library(gsignal)
Har.anomaly = tibble(dt=global.ts$dt.mnth[1:n_shrt],
                        y.1058=Re(ifft(FFT_filt$spc)))
Har.anomaly%>%ggplot(aes(x=dt,y=y.1058))+geom_line()
#_________
FFT_filt%>% subset(idx < 1058/2)%>%arrange(desc(amp))
#---------
# ext backwards
my_signal=Har.anomaly$y.1058
library(gsignal)
extend_backwards <- function(x, n_back = 3000) {
  N <- length(x)
  N_ext=N+n_back
  t_ext= seq(from=-(n_back-1),N)
  fft_orig <- fft(x) /N
  idx_ext=((t_ext - 1) %% N) + 1
  fft_ext <- fft_orig[idx_ext]
  y_numeric=Re(ifft(fft_ext) * N_ext)
  df=tibble(t_ext=t_ext,
            y=y_numeric)
  return(df)
}
#====================
Y.ext=extend_backwards(my_signal)
Y.ext%>%ggplot(aes(x=t_ext))+geom_line(aes(y=y))
# compare with solar phase locked trend;
#from "SolarPower.phaselocked_Anomalies.R"
df=readRDS("data/NOAA.ocean.anomalies.rds") #2068 obs
global.ts<-ts.list$global.ts # NROW 2116
Ocean_anoma<-df$data
library(itsmr)
library(lubridate)
M=c("season",12,"season",6)
Ocean_anomaly=Ocean_anoma%>%  mutate(anom=Resid(anoma.mean,M),
                                     trd3=trend(anom,3),
                                     res3=Resid(anom,3),
                                     Year=year(date),
                                     mnth=month(date))
# combine solar_monthly with ocean anomaly
solar_mnthly=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))
# combine solar and ocean data
Ocean_solar_anom=Ocean_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,anoma.mean,"Solar_Variation"=SI,trd3)
# anchor timeseries solar power
library(splines)
# find zero crossings of SI
find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}
nodes <- find_nodes(solar_mnthly$SI)
# 1. Identify the time coordinates of your zero-crossings
zero_crossing_times <- solar_mnthly$dt.mnth[find_nodes(solar_mnthly$SI)]
# 2. Fit a continuous spline locked at these exact nodes
trend_fit <- lm(anoma.mean ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                data = Ocean_solar_anom)
# 3. Extract the stable Baseline Trend and visualise
Ocean_solar_anom$Clean_Baseline_Trend <- predict(trend_fit)
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Clean_Baseline_Trend ),col=2)+
  labs(title="Clean_Baseline_Trend",
       subtitle = "extracted with spline fit solar power degree= 3")
range(Ocean_solar_anom$dt.mnth)
Ocean_anom<-Ocean_solar_anom%>%left_join(Low.pass,by="dt.mnth")%>%
  dplyr::select(dt.mnth,
                "BaselineTrend"=Clean_Baseline_Trend,
                y.ext)
Ocean_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=BaselineTrend))+
  geom_line(aes(y=y.ext),col=4)+
  geom_line(aes(y=BaselineTrend-y.ext),col=2,linewidth = 1.5)+
  labs(x="",y="ocean anomalies K",title = "Difference (red) : Solar-Power locked trend (black)\n & Extrapolated Periodic Anomalies (blue)")
