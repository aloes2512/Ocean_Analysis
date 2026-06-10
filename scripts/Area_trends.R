# Area trends tanh
# Define the tanh model for nls
# A: Amplitude of the transition
# k: Steepness of the warming/cooling trend
# x0: The center point of the transition (e.g., a specific year/month)
# b: Vertical offset (baseline)
#1. format data
library(terra)
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
idx.mx= idx.mx[1:8]
idx.mx2=idx.mx[1:4]
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


#==========
#===========
library(minpack.lm)
my_data=Trop.Pac_tmp.anom%>%mutate(time_index=1:N-1)
my_data%>%ggplot(aes(x=dates,y=res.pac))+geom_line()
res.pac=Trop.Pac_tmp.anom$res.pac
fit_tanh <- nlsLM(res.pac~ A * tanh(k * (time_index - x0)) + b,
                  data = my_data,
                  start = list(A = 0.3, k = 0.001, x0 = 1000, b = 0.08))
# Generate predictions
my_data$pred_tanh.res <- predict(fit_tanh)
my_data$pred_tanh<-NULL
# Plot the residuals and the tanh fit

ggplot(my_data, aes(x = dates)) +
  geom_line(aes(y = res.pac), color = "grey70", alpha = 0.6) + # Raw residuals
  geom_line(aes(y = pred_tanh.res), color = "firebrick", linewidth = 1) + # Tanh trend
  labs(title = "Equatorial Pacific: \nTanh Trend on Linear Residuals",
       y = "Residual Anomaly",
       x = "Year") +
  theme_minimal()
# trend predicted:
colnames(my_data)
my_data=my_data%>%mutate(pred.tanh=pred_tanh.res+lin.trnd)
my_data%>%ggplot(aes(x=dates,y=pred.tanh))+geom_line()
coef(fit_tanh)
round(coef(fit_tanh),2)
