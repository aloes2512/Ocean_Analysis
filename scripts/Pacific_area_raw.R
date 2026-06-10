library(terra)
library(tidyverse)
ersst= terra::rast("~/projects/Ocean_analysis/data/gistemp1200_GHCNv4_ERSSTv5.nc")
dim(ersst) # 90 80 1754
dates=time(ersst) # rangge "1880-01-15" "2026-02-15"
N= length(dates)
# multible polygons 3 polygons (from matrices)
coords1 <- rbind(c(-179, 5), c(-120, 5), c(120, -5), c(-179, -5), c(-179, 5))
coords2 <- rbind(c(-150, 0), c(-100, 0), c(-100, -5), c(-150, -5), c(-150, 0))
coords3 <- rbind(c(50, 2), c(100, 2), c(100, -3), c(50, -3), c(50, 2))
# Combine into one SpatVector
polys <- vect(list(coords1, coords2, coords3),
              type = "polygons",
              crs = crs(ersst))
poly_trop_pac <- vect("MULTIPOLYGON(((110 20, 180 20, 180 -20, 110 -20, 110 20)),
                               ((-180 20, -90 20, -90 -20, -180 -20, -180 20)))",
                      crs = crs(ersst))
# crop trop_pac
ersst_trop.pac=crop(ersst, poly_trop_pac)
temp_anom.trop.pac=terra::global(ersst_trop.pac,"mean",na.rm=T)$mean
Trop.Pac_tmp.anom=tibble(dates=time(ersst),temp_anom=temp_anom.trop.pac)
library(itsmr)
M.trd= c("season",12,"trend",1)
Trop.Pac_tmp.anom=Trop.Pac_tmp.anom%>%mutate(res.pac=Resid(temp_anom,M.trd),
                                             lin.trnd=trend(temp_anom,1),
                                             anom.smth=smooth.fft(temp_anom,f= 0.02),
                                             res.smth=smooth.fft(res.pac,f= 0.02))
Trop.Pac_tmp.anom%>%ggplot(aes(x=dates,y=anom.smth))+
  geom_line(col="grey")+
  geom_line(aes(y=res.smth),col=2)+
  labs(x= "month", y="tempanomaly K",
       title =" Tropical Pacific\nMonthly & Area averaged",
       subtitle = "low-pass filtered(grey) & linear residuals(red)"
  )
time_index=1:N-1
pacific_data<-Trop.Pac_tmp.anom # dates temp_anom res.pac lin.trnd anom.smth res.smth
C= mean(pacific_data$anom.smth) # 0.08
A=sd(pacific_data$anom.smth) # 0.32
sd(pacific_data$res.smth) # 0.18
temp_anoma.smth=pacific_data$anom.smth
res.smth=pacific_data$res.smth
FFT.res= tibble(idx= 1:N-1,spc=fft(res.smth),amp= Mod(spc))
idx.mx=FFT.res%>% arrange(desc(amp))%>% pull(idx)
idx.mx[1:8]
idx.mx2=idx.mx[1:4]
idx.mx[9:17]
# decadal
idx.deca=FFT.res%>% arrange(desc(amp))%>% pull(idx)%>%{.[9:16]}

FFT.mxharm=FFT.res%>% mutate( spc= ifelse(idx %in% idx.mx,spc,0i))
FFT.mxharm2=FFT.res%>% mutate( spc= ifelse(idx %in% idx.mx2,spc,0i))
FFT.deca=FFT.res%>% mutate( spc= ifelse(idx %in% idx.deca,spc,0i))

library(gsignal)
tibble(dates=time(ersst),y.max= Re(ifft(FFT.mxharm$spc)))%>%
  ggplot(aes(x=dates,y=y.max))+geom_line()
tibble(dates=time(ersst),y.max= Re(ifft(FFT.mxharm2$spc)))%>%
  ggplot(aes(x=dates,y=y.max))+geom_line()
tibble(dates=time(ersst),y.max= Re(ifft(FFT.deca$spc)))%>%
  ggplot(aes(x=dates,y=y.max))+geom_line()
# nls optim 1
FFT.res[3,]
phi1=Arg(FFT.res[[3,2]])# 0.3380
Amp=FFT.res[[3,3]]/N # 0.05
fit_trop_pac <- nls(
  sgnl ~ Amp* sin((2*pi/T) * t + phi) + C+ d*t,
  data = list(sgnl = Trop.Pac_tmp.anom$res.smth, t = (1:N) - 1),
  start = list(T=N/2, phi = phi1, C = 0,d=3e-7),
  control = nls.control(maxiter = 1000)
)
AIC(fit_trop_pac) # -1307.056
coef1=coef(fit_trop_pac)
T.fit=coef1[[1]]
phi.fit=coef1[[2]]
t=(1:N)-1
trop.pac_period=coef(fit_trop_pac)[[1]]
Y.fit=tibble(dates=time(ersst),
             y1.fit=Amp*sin((2*pi/T.fit)*t)+coef1[[3]]+coef1[[4]]*t)
Y.fit%>%ggplot(aes(x=dates,y=y1.fit))+geom_line()
# second group 10:16
# start
FFT.res[15,]
Amp2=FFT.res[[15,3]]
phi2=Arg(FFT.res[[15,2]])
T2= 132