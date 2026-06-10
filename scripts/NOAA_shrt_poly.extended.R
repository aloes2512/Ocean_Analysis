library(tidyverse)
NOAA.Ocean.anomalies%>%colnames()#dt.mnth mean.anom median.anom
N=NOAA.Ocean.anomalies%>%NROW() # 494
M.ssn=c("season",12)
M1= c("season",12,"trend",1)
library(itsmr)
Ocean_mean.anoma=NOAA.Ocean.anomalies%>%mutate(yr.anom=Resid(mean.anom,M.ssn),yr.anom=smooth.fft(yr.anom,0.1),
                                               trd1=trend(yr.anom,1),
                                               trd2=trend(yr.anom,2),
                                               trd3=trend(yr.anom,3),
                                             trd4=trend(yr.anom,4))

head(dt.mnth)
dt.mnth=Ocean_mean.anoma$dt.mnth
seq.Date(first(dt.mnth)-30,by ="-1 month",length=741)%>%range()
ord.seq.dts=seq.Date(as.Date("1920-03-01"),by = "1 month",length=741)
range(ord.seq.dts)

ext.dates=c(ord.seq.dts,dt.mnth)
length(ext.dates)
Ocean_mean.anoma%>%ggplot(aes(x=dt.mnth))+
                        geom_line(aes(y=trd1),col=1)+
                        geom_line (aes(y=trd2),col=3)+
                        geom_line (aes(y=trd3),col=5)+
                        geom_line (aes(y=trd4),col=2)+
                        geom_line(data=Mean.Res,aes(y=linres),col=4)+
                        geom_line(data=Mean.Res,aes(y=smth.res))
Mean.Res= Ocean_mean.anoma%>% mutate(linres=Resid(mean.anom,M1),
                              smth.res=smooth.fft(linres,f= 0.01))
Mean.Res%>%ggplot(aes(x=dt.mnth,y=yr.anom))+geom_line()+
                    geom_line(aes(y=smth.res),col=2)
FFT=tibble(idx=1:N,spc=fft(Mean.Res$smth.res),amp= Mod(spc))
idx.ord=FFT%>%arrange(desc(amp))%>%pull(idx)
idx.mx=idx.ord[1:4]
N/3 # 164 month 13.66667 years
N/2 # 247 month 20.58333 years
library(gsignal)
FFT_filtered <- FFT %>%
  mutate(spc = if_else(idx %in% idx.mx, spc, 0))
tibble(t=dt.mnth,y.mx=Re(ifft(FFT_filtered$spc)))%>%
  ggplot(aes(x=t,y=y.mx))+geom_line()
# perplexity example
N <- NROW(NOAA.Ocean.anomalies)
spc <- FFT$spc%>%as.vector()
ext_factor <- 2.5
past_len <- round(N * (ext_factor - 1)) : 741
t_ext <- seq(-past_len + 1, N)  # Backward + original
head(t.mnth)
idx_ext <- ((t_ext - 1) %% N) + 1
spc_ext=spc[idx_ext]
#using gsignal
y.ext=ifft(spc_ext)
k <- (0:(N-1)) / N  # Normalized freqs
signal_ext <- Re( sapply(idx_ext, \(i) sum(spc * exp(1i * 2 * pi * k * (i - 1)) / N ) ))
 length(signal_ext)
 tibble(dts=ext.dates,y=signal_ext)  %>% ggplot(aes(x=dts,y=y))  +geom_line()


