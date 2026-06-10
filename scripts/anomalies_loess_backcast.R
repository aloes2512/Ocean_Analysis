Ocean_anomaly.phaselocked=readRDS("data/Ocean_anomaly.phaselocked.rds")
library(tidyverse)
library(itsmr)
df=Ocean_anomaly.phaselocked$Ocean_anomaly_detrend
M <- c("season", 12, "season", 6,"trend",1)
df=df%>%mutate(Detrnd.seasonalyzed_anomaly=Resid(Detrended_Anomaly,M))
#'long-term signal is not behaving like a neat,
# rigid sine wave, Fourier methods are the wrong tool'
df=df%>%dplyr::select(Time,"anomaly"=Detrnd.seasonalyzed_anomaly)%>%
  mutate(time_index=1:NROW(df)-1)
# 1. Broad trend (looks at 50% of the entire data span at a time)
loess_macro <- loess(anomaly ~ time_index, data = df, span = 0.50, degree = 2)

# 2. Medium trend (looks at 20% of the data span, closer to multi-decadal cycles)
loess_med <- loess(anomaly ~ time_index, data = df, span = 0.20, degree = 2)

# Calculate the smooth paths
df$macro_trend <- predict(loess_macro)
df$med_trend   <- predict(loess_med)
# 3. Plot to inspect the underlying geometry
df%>%ggplot(aes(x=Time))+
  geom_line(aes(y=macro_trend),col="darkblue", lwd = 1.5)+
  geom_line(aes(y=med_trend),col=2,lwd=1.2)
#=======================
FFT.loess= tibble(idx=1:N-1,
                  spc=fft(sign.loess),
                  amp=Mod(spc))

FFT.loess%>%arrange(desc(amp))%>%subset(idx<N/2)%>%{.[1:20,]}
# select 10 strongest harmonics using hr
df%>%dplyr::select(-macro_trend) %>%
  mutate( mx.period=hr(med_trend,N/c(1:10)),
          men.per=hr(med_trend,N/4) ) %>%
  ggplot(aes(x=Time,y=mx.period))  +geom_line()+
  geom_line(aes(y=men.per),col=2)+
  labs(title="loess smoothed trends",
       subtitle = "Long range periods;related to 43 years period")
#-----------
# medium long periods
df%>% dplyr::select(-macro_trend) %>%
  mutate( period.2nd=hr(med_trend,N/c(12:15)), # 11.5 to 14.4 years
          schwabe=hr(med_trend,N/15)) %>%
  ggplot(aes(x=Time,y=period.2nd))  +
  geom_line()+geom_line(aes(y=schwabe),col=2)+
  labs(title = "Medium Periods",subtitle = "includes Schwabe-cycle")
#------------------------
# decadal periods ~8 years
df%>% dplyr::select(-macro_trend) %>%
  mutate( period.3rd=hr(med_trend,N/c(20:24)),
          med.per=hr(med_trend,N/20)) %>%
  ggplot(aes(x=Time,y=period.3rd))  +geom_line()+
  geom_line(aes(y=med.per),col=2)+
  labs(title="Decadal Cycles",subtitle= "reference 8.6 years")

#=======================
# 1. Gather parameters from the 25 clean, verified FFT coefficients
# N is the length of your original data vector
colnames(df) # "Time"       "anomaly"    "time_index" "med_trend"
N <- length(df$med_trend) #2068
FFT.loess=tibble(idx=1:N-1,spc=fft(df$med_trend),amp=Mod(spc))
verified_modes <- tibble(idx = 1:25) %>%
  mutate(
    # Pull directly from your clean loess spectrum
    spc   = FFT.loess$spc[idx],
    amp   = Mod(spc) * (2 / N),
    phase = Arg(spc),
    # Frequency per time step (months/years depending on your df structure)
    freq  = (idx - 1) / N
  )
# 2. Build a stable timeline from the year 1600 to 2026
# (Assuming monthly steps; adjust the conversion to match your time index)
start_idx <- 1 - (250 * 12) # Go back 250 years before your 1850 start
extended_timeline <- start_idx:N # -2999 to 2068

# 3. Project the waves backward stably
pure_backcast <- numeric(length(extended_timeline))

for(i in 1:nrow(verified_modes)) {
  pure_backcast <- pure_backcast +
    verified_modes$amp[i] * cos(2 * pi * verified_modes$freq[i] * extended_timeline + verified_modes$phase[i])
}

# 4. Map back to calendar years for a clean plot
plot_years <- 1850 + (extended_timeline / 12)
plot(plot_years, pure_backcast, type="l", col="black", xlab="Year", ylab="True Harmonic Backcast")
abline(v=1850, col="red", lty=2) # Marks where your real data begins
# Historical data
# Option A: Downloading a standard historical reconstruction text file from NOAA
proxy_url <- "https://www.ncei.noaa.gov/pub/data/paleo/contributions_by_author/mann2008/recons/nh-cps-multi-proxy.txt"

# Read the historical data (skipping metadata headers depending on the file layout)
historical_proxies <- read.table(proxy_url, header = TRUE, skip = 80)

# Option B: Using the paleoclimate repository API (Lipd)
# install.packages("remotes")
# remotes::install_github("LinkedEarth/LiPD-Utilities", sub = "R")
library(lipdR)

# This allows you to pull down marine sediment or coral records
# matching your specific oceanic regions of interest
