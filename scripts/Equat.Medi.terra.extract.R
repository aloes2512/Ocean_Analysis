# Anomalies_detrend_harmonic
library(tidyverse)
# data from terra and formated
library(terra)
library(lubridate)
weighted_sst_stack<-rast("data/weighted_sst_stack.tif")
# 2. Convert to relative weights (0 to 1 scale)
# This divides every cell by the area of an equatorial cell
active_weights <- rast("data/active_weights.tif")
# Now your denominator is perfectly synced to the available data
sum_of_ocean_weights <- global(active_weights, fun="sum", na.rm=TRUE)
ocean.anom.sum=unlist(global(weighted_sst_stack, fun="sum", na.rm=TRUE))
ocean_anom_noaa <- ocean.anom.sum / as.numeric(sum_of_ocean_weights)

#1. The Bounding Box / Extent Generator
make_region_extent <- function(xmin = -180, xmax = 180, ymin = -90, ymax = 90) {
  terra::ext(xmin, xmax, ymin, ymax)
}
equatorial_ext <- make_region_extent(-180, 180, -5, 5)
#2. The Cropping and Spatial Weighting Function
extract_weighted_region <- function(raster_brick, region_ext, mask_poly = NULL) {
  cropped_brick <- terra::crop(raster_brick, region_ext)
}

equatorial_stack=extract_weighted_region(weighted_sst_stack,equatorial_ext )
equatorial_weigths=extract_weighted_region(active_weights,equatorial_ext)
# sum to get ts
equator.anom.sum=unlist(global(equatorial_stack, fun="sum", na.rm=TRUE))
sum_equator_weights <- global(equatorial_weigths, fun="sum", na.rm=TRUE)
equator_anom<- equator.anom.sum / as.numeric(sum_equator_weights)
Equator_mean.anom=tibble(dt.mnth=time(weighted_sst_stack),
                         equator_anomaly=equator_anom-mean(equator_anom))
Equator_mean.anom%>%ggplot(aes(x=dt.mnth,y=equator_anomaly))+geom_line()
colnames(Equator_mean.anom)
Equator_mean.anom=rename(Equator_mean.anom,"Anomaly"=equator_anomaly)
# Gibraltar -5° W 36 °N
# Venedig    12°E 45° N
#Beirut      34° E 35.5° N
#Al-Agheila  19° E 30° N
Medi_ext=make_region_extent(-5, 34, 30, 36)
Medi_stack=extract_weighted_region(weighted_sst_stack,Medi_ext )
Medi_weigths=extract_weighted_region(active_weights,Medi_ext)
Medi.anom.sum=unlist(global(Medi_stack, fun="sum", na.rm=TRUE))
sum_Medi_weights <- global(Medi_weigths, fun="sum", na.rm=TRUE)
Medi_anom<- Medi.anom.sum / as.numeric(sum_Medi_weights)
Medi_mean.anom=tibble(dt.mnth=time(weighted_sst_stack),
                         Medi_anomaly=Medi_anom-mean(Medi_anom))
Medi_mean.anom%>%ggplot(aes(x=dt.mnth,y=Medi_anomaly))+geom_line()
Medi_mean.anom=rename(Medi_mean.anom,"Anomaly"=Medi_anomaly)

#==================
library(itsmr)

extend.preind=function(area.ts){
  M.ssn=c("season",6,"season",12)
  N=NROW(area.ts)
  N_base=floor(N/2)
  ts=area.ts%>%mutate(Anomaly=Resid(Anomaly,M.ssn)) #eliminate season
  preind.lmt= ts$dt.mnth[N_base]
  # extend circular
  circular_indices <- ((0:(N - 1)) %% N_base) + 1 #N_base= 1058; circular-length = N
  Preind.ext=tibble(dt.mnth=ts$dt.mnth,
                    ext.preind=ts$Anomaly[circular_indices])%>%
                    mutate(anoma.hr1_6=hr(ext.preind,N/1:6))
  return(Preind.ext)
  }
Medi.preind=extend.preind(area.ts =Medi_mean.anom )

Medi.preind%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=ext.preind),col="grey")+
  geom_line(aes(y=anoma.hr1_6),col=2)
