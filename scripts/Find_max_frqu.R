NINO.34=readRDS("data/NINO_day.spline.rds")%>%
  mutate(yr.mnth=as.yearmon(date))%>%
  group_by(yr.mnth)%>%
  summarise(T.mnth=mean(tmp_C))
mdl=c("season",12,"trend",1)
NINO.34=NINO.34%>%  mutate(T=Resid(T.mnth,mdl))
n=NROW(NINO.34)
sst_deviations=NINO.34$T
#======
compute_fft_peaks <- function(sst_deviations, sampling_rate, n_peaks=5) {
  n <- length(sst_deviations)
  fft_result <- fft(sst_deviations)
  frequencies <- seq(0, sampling_rate/2, length.out = n/2)
  magnitude <- Mod(fft_result)[1:(n/2)]^2 / n
  
  # Find top peaks
  peak_indices <- order(magnitude, decreasing = TRUE)[1:n_peaks]
  data.frame(
    frequency = frequencies[peak_indices],
    amplitude = magnitude[peak_indices]
  )
}
#========
frqs=compute_fft_peaks(sst_deviations, 1, n_peaks=5)$frequency
1/n # 0.0005396654
frqs*n
FFT=tibble(idx=1:n,spc=fft(NINO.34$T),amp=Mod(spc))
require(zoo)
idx.mx=FFT%>%arrange(desc(amp))%>%head(12)%>%subset(idx<n/2)%>%pull(idx)
(idx.mx-1)/n
NINO.34%>%ggplot(aes(x=yr.mnth,y=T))+geom_line()
