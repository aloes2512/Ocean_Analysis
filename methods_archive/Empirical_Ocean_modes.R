Ocean_solar_anom=readRDS("data/Ocean_solar_anom.rds")
Ocean_solar_anom=Ocean_solar_anom%>%
  mutate(smth.res=smooth.fft(Solar_Retained_Residuals,f=0.02))
library(itsmr)
N=NROW(Ocean_solar_anom)
M=c("season",12,"season",6)
colnames(Ocean_solar_anom)
Ocean_trends=Ocean_solar_anom%>%
  mutate(trend_difference=trd3-Clean_Baseline_Trend)%>%
  dplyr::select(dt.mnth,"y.predict"=Clean_Baseline_Trend,trd3,trend_difference)
#visualise surprising fft ordered
y.predict=Ocean_trends$y.predict
FFT.resd.locked=tibble(idx=1:N-1,
                       spc=fft(y.predict),
                       amp=Mod(spc))
FFT.resd.locked%>%subset(idx<N/2)%>%arrange(desc(amp))%>% head(20)
#-----
Ocean_trends=Ocean_trends%>%mutate(trend_diff_macro=smooth.fft(trend_difference, f = 0.008))
trend_diff_macro=Ocean_trends$trend_diff_macro
phase_locked_residuals=smooth.fft(Ocean_solar_anom$Solar_Retained_Residuals,f=0.008)
# add the two residuals to include solar locked part
total_backcast_residuals <- phase_locked_residuals + trend_diff_macro
Ocean_residuals.complete=tibble(dt.mnth=Ocean_trends$dt.mnth,
                                total_backcast_residuals <- phase_locked_residuals + trend_diff_macro
)
Ocean_residuals.complete%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=total_backcast_residuals))+
  labs(title="Complete Ocean Residuals",
       subtitle = "added: solar phase locked / trend-difference macro")
#=============
#' analyze the modes'
# analyse the modes
library(EMD)
colnames(Ocean_residuals.complete)<-c("dt.mnth","total_backcast_residuals")
# Extract the Intrinsic Mode Functions (IMFs) from your complete residuals
emd_modes <- emd(Ocean_residuals.complete$total_backcast_residuals,
                 tt = Ocean_residuals.complete$dt.mnth,
                 max.imf = 3) # Limit to the top 2-3 broad physical scales

# Look at the separated modes
tibble(tt=Ocean_residuals.complete$dt.mnth,
       y.22year=emd_modes$imf[,1])%>%
  ggplot(aes(x=tt))+geom_line(aes(y=y.22year)) #will be your faster ~22-year Hale wave

y.60year=emd_modes$imf[,2] #will be your slower ~60-year Ocean Current wave
tibble(tt=Ocean_residuals.complete$dt.mnth,
       y.60year=emd_modes$imf[,2])%>%
  ggplot(aes(x=tt,y=y.60year))+geom_line()

tibble(tt=Ocean_residuals.complete$dt.mnth,
       y.xyear=emd_modes$imf[,3])%>%
  ggplot(aes(x=tt,y=y.xyear))+geom_line()
#=====
# Create a data frame for plotting the separated gears
Oceans_modes <- tibble(
  dt.mnth = Ocean_residuals.complete$dt.mnth,
  Raw     = Ocean_residuals.complete$total_backcast_residuals,
  Mode_1  = emd_modes$imf[, 1],  # The faster wave component (~22-year Hale cycle)
  Mode_2  = emd_modes$imf[, 2],  # The slower wave component (~60-year Ocean current)
  Mode_3  = emd_modes$imf[, 3],  #  100 years wave
  Trend   = emd_modes$residue    # Any remaining non-periodic baseline
)

# Plot Mode 2 to see your pure ocean oscillation
Oceans_modes %>%
  ggplot(aes(x = dt.mnth, y = Mode_2)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_line(aes(y=Mode_1),col=2)+
  labs(title = "Isolated Deep Ocean Mode (~60-year wave)\n and ~ 22-year wave")
# cross corelation of 60 year wave and 100 year wave
Oceans_modes%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Mode_1),col=2)+
  geom_line(aes(y=Mode_2),col=3)+
  geom_line(aes(y=Mode_3),col=4)+
  geom_line(aes(y=Trend))
Oceans_modes%>%mutate(long.period=Mode_2+Mode_3)%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Mode_1),col=2)+
  geom_line(aes(y=long.period),col=4)+ # long.period == Mode_2+Mode_3
  geom_line(aes(y=Trend))
# 1. Combine the phase-split quadrature modes into a single physical engine
Ocean_reconstructed <- Oceans_modes %>%
  mutate(
    Unified_MultiDecadal = Mode_2 + Mode_3,
    Hale_Cycle           = Mode_1
  )

# 2. Plot the unified components to check the final structural layers
Ocean_reconstructed %>%
  ggplot(aes(x = dt.mnth)) +
  geom_line(aes(y = Unified_MultiDecadal), color = "blue", linewidth = 1) +
  geom_line(aes(y = Hale_Cycle), color = "red", alpha = 0.7) +
  geom_line(aes(y = Trend), color = "black", linewidth = 1.2) +
  labs(
    title = "Unified Physical Engines \nof the Background Signal",
    y = "Anomalies",
    x = "Timeline"
  )
saveRDS(Ocean_reconstructed, "data/Ocean_modal_modes.rds")
######
