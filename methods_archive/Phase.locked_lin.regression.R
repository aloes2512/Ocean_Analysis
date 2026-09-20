#Phase-Locked Regression  Solar data: 260522 days
# data from "data/S_power.rds" == processed SILSO counts
# and updated "data/NOAA.ocean.anomalies.rds" monthly
#
#1. Anchor Alignment:identifying the local extrema (minima/maxima) of your
#Solar Power timeseries calculated with SP_negbin.power.R
library(tidyverse)
library(zoo)
library(lubridate)
df=readRDS("data/S_power.rds")
df_monthly <- df %>%
  mutate(
    date = yearmon(dt.mnth),
    Time=dt.mnth,
    Year = year(date),
    Month = month(date)
  ) %>%
  group_by(Year, Month,Time) %>%
  summarize(
    TSI = mean(TSI, na.rm = TRUE),
    .groups = "drop" # Automatically ungroups for you
  )
solar_mnthly=df_monthly%>%mutate(SI=TSI-mean(TSI))
solar_mnthly<-solar_mnthly%>%mutate(dt.mnth=Year+(Month-1)/12)
#=======
Ocean=readRDS("data/NOAA.ocean.anomalies.rds")$data
Ocean_anomaly=Ocean%>%dplyr::select(-anoma.median)

Ocean_solar_anomaly=Ocean_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,"Ocean_Anomaly"=anoma.mean,"Solar_Variation"=SI)

# --- 1. Data Preparation (Mockup of your data structure) ---
# Ocean_solar_anom has columns: dt.mnth, Ocean_Anomaly, Solar_Variation
# Solar_Variation: quadratic-transformed negative binomial model

# Identify Solar Nodes (Minima) to define the 10.9-year windows
# find LOCAL TRENDS of Solar Variation series
find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}

nodes <- find_nodes(Ocean_solar_anomaly$Solar_Variation)
mean(diff(nodes),na.rm=T) # 130.73 month
# add group index ~ Schwabe cycle nr
Ocean_solar_anomaly$Solar_Phase_Group <- cut(1:nrow(Ocean_solar_anomaly), breaks = c(0, nodes, nrow(Ocean_solar_anomaly)), labels = FALSE)

# --- 2. Range-Limited Trend Calculation ---
calc_segmented_trends <- function(data) {
  fit <- lm(Ocean_Anomaly ~ dt.mnth, data = data)

  data.frame(
    dt.mnth = data$dt.mnth,
    Segment_Trend = predict(fit),
    Slope = coef(fit)[2]
  )
}

local_trends <- Ocean_solar_anomaly %>%
  group_by(Solar_Phase_Group) %>%
  do(calc_segmented_trends(.)) %>%
  ungroup()
# calculate a trend for each solar cycle 'bin'

# Apply the function across each solar cycle window
local_trends <- Ocean_solar_anomaly%>%
  group_by(Solar_Phase_Group) %>%
  do(calc_segmented_trends(.))

# --- 3. Detrending and Final Output ---
Ocean_solar_anom_final <- Ocean_solar_anomaly %>%
  left_join(local_trends, by = c("dt.mnth", "Solar_Phase_Group")) %>%
  mutate(
    Detrended_Anomaly = Ocean_Anomaly - Segment_Trend
  )

# --- 4. Visualization ---

Ocean_solar_anom_final%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Ocean_Anomaly),col="grey")+
  geom_line(aes(y=Segment_Trend,col="Segment_Trnd"))+
  labs(x="",title = "Ocean Temp Anomalies \nwith Solar-Anchored Segmented Trends")

#3.`Residual Analysis:
library(itsmr)

Ocean_solar_anom_final%>%ggplot(aes(x=dt.mnth,y=Detrended_Anomaly))+
  geom_line(col="grey")+
  labs(title = "Ocean: Detrended Anomaly",
       subtitle = "Linear Trends dependent on Schwabe Cyles")
N.mnth=NROW(Ocean_solar_anom_final) #2116
Ocean_solar_anom_final=Ocean_solar_anom_final%>%
  mutate(res.smth=smooth.fft(Detrended_Anomaly,f=0.0085),
         hr.lowfrq=hr(Detrended_Anomaly,N.mnth/1:10))

Ocean_solar_anom_final%>%  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=res.smth,col="f=0.0085"))+
  geom_line(aes(y=hr.lowfrq,col="hr.lowfrq"))+
  labs(title = "Ocean: Local Detrended Anomaly",
       subtitle = "phase_locked to solar power\nsmoothed f & hr: 10 lowest f")
FFT.detrnd.loc=tibble(idx=1:N.mnth-1,
                      spc=fft(Ocean_solar_anom_final$Detrended_Anomaly),
                      amp=Mod(spc))

FFT.4mx=FFT.detrnd.loc%>%arrange(desc(amp))%>%{.[1:8,]} # 1 year, 6 mnth, 3.8 years
idx.4mx=FFT.4mx%>%pull(idx) # 353 1763 1940  176  352 1764 1939  177
FFT.0=FFT.detrnd.loc%>% subset(!idx %in% idx.4mx)%>% mutate(spc=0i,amp=0)
FFT.4mx=bind_rows(FFT.4mx,FFT.0)%>%arrange(idx)
library(gsignal)
Y.detrnd.4=tibble(Time=Ocean_solar_anom_final$dt.mnth,
                  y.4mx=Re(ifft(FFT.4mx$spc)),
                  y.smth=smooth.fft(y.4mx,f= 0.06)) # low.pass cutoff period 10.34 years
Y.detrnd.4%>%ggplot(aes(x=Time,y=y.4mx))+
  geom_line(col="grey")+
  geom_line(aes(y=y.smth),col=2)+
  labs(y="harmonics amp",title = "Resids Phase locked regression",
       subtitle = "locked to Solar Power Schwabe cycle (red)")
N.mnth*0.06 # cutt-off 124.08 month 10.3 years
#.4 Trend analysis
Ocean_solar_anom_final%>%
  ggplot(aes(x=dt.mnth,y=Segment_Trend))+
  geom_line()+
  geom_smooth()+
  labs(x="month",
  title = "Ocean_anomaly trends locked to\n Solar Power Schwabe cycle ",
  subtitle= " gam smoothed formula: y ~ s(x, bs =`cs`)")
# Integrate the linear trends
Integr.trnds <- Ocean_solar_anom_final %>%
  mutate(
    # True integration: Cumulative sum of (Rate * TimeStep)
    # Since step is 1 month, we just sum the slopes
    Integrated_trends = cumsum(Segment_Trend)
  )
Integr.trnds%>%
  ggplot(aes(x=dt.mnth,y=Integrated_trends))+
  geom_line()+labs(x="",title="Ocean Anomaly",
                   subtitle =  "integrated segmented Slopes\n locked to Schwabe Cycle")
# needs to be scaled
# Assuming 'recon_signal' is your cumsum(gam_slopes)
# and 'Ocean_Anomaly' is your original 1881-2026 data
integrated=Integr.trnds$Integrated_trends
target=Ocean_solar_anom_final$Ocean_Anomaly
scale_integrated_signal <- function(integrated, target) {
  # 1. Zero-center the integrated signal
  integrated <- integrated - mean(integrated, na.rm = TRUE)

  # 2. Scale it so its 'spread' matches the target data
  scaled <- integrated * (sd(target, na.rm = TRUE) / sd(integrated, na.rm = TRUE))

  # 3. Shift it so its 'average' matches the target data
  final <- scaled + mean(target, na.rm = TRUE)

  return(final)
}

Ocean_solar_anom_final$Integrated_Scaled <- scale_integrated_signal(
  Integr.trnds$Integrated_trends,
  Ocean_solar_anom_final$Ocean_Anomaly
)

Ocean_anomaly_trend=Ocean_solar_anom_final%>%
  dplyr::select(-Solar_Phase_Group)
Ocean_anomaly_trend%>%
  ggplot(aes(x=dt.mnth,y=Integrated_Scaled))+
  geom_line()+
  labs(x="",title = "Ocean Anomaly Trend ",
       subtitle = "integrated from lin trends locked to \nsolar power Schwabe Cycle")
FFT.seg_trd=tibble(idx=1:N.mnth,
                   spc=fft(Ocean_anomaly_trend$Integrated_Scaled),
                   amp= Mod(spc))
FFT.seg_trd%>%subset(idx<100)%>%arrange(desc(amp))
N.mnth/2/12 # 88 years 2116 month
FFT.seg_trd=FFT.seg_trd%>%mutate(period.yr=N.mnth/(idx*12))

#============
# Filter for the post-1975 'Double Warming'?
df_modern=Ocean_anomaly_trend%>%subset(dt.mnth>1975)
df_modern$dt.mnt%>%range() # 1975 2026.25

# 1. Calculate the Solar Slope (Rate of natural recovery)
modern_slope <- lm(Integrated_Scaled ~ dt.mnth, data = df_modern)$coefficients[2]
modern_slope # 0.01155
total_trend.slope <- lm(Integrated_Scaled ~ dt.mnth, data=Ocean_anomaly_trend)$coefficients[2]
total_trend.slope # -0.00266
slope.trend_change=modern_slope-total_trend.slope #0.014
ocean_slope <- lm(Ocean_Anomaly~dt.mnth,data=Ocean_solar_anomaly)$coefficients[2]
ocean_slope # 0.0028
Ocean_anomaly.phaselocked= list(Ocean_anomaly_detrend=Ocean_anomaly_trend,
                                periods=FFT.seg_trd,
                                time.range=range(Ocean_anomaly$Time))
saveRDS(Ocean_anomaly.phaselocked,"data/Ocean_anomaly.phaselocked.rds")
