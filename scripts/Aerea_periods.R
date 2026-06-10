library(terra)
ersst= terra::rast("~/projects/Ocean_analysis/data/gistemp1200_GHCNv4_ERSSTv5.nc")
# Or regional (e.g., tropical Pacific polygon) to be save make coordinates first
coords <- rbind(c(-179, 5), c(-120, 5), c(120, -5), c(-179, -5), c(-179, 5))
poly_trop_pac <- vect(coords, type="polygon", crs = crs(ersst))
dates=time(ersst) # rangge "1880-01-15" "2026-02-15"
# multible polygons
# Example: 3 polygons (from matrices)
coords1 <- rbind(c(-179, 5), c(-120, 5), c(120, -5), c(-179, -5), c(-179, 5))
coords2 <- rbind(c(-150, 0), c(-100, 0), c(-100, -5), c(-150, -5), c(-150, 0))
coords3 <- rbind(c(50, 2), c(100, 2), c(100, -3), c(50, -3), c(50, 2))
# Combine into one SpatVector
polys <- vect(list(coords1, coords2, coords3),
              type = "polygons",
              crs = crs(ersst))

# Save all at once
writeVector(polys, "multipolygons.gpkg", overwrite=TRUE)
# get it from disc
polys<-vect("multipolygons.gpkg")
# Or more detailed multi-polygon for full tropical Pacific (e.g., 20°S–20°N, 110°E–90°W)
poly_trop_pac2 <- vect("MULTIPOLYGON(((110 20, 180 20, 180 -20, 110 -20, 110 20)),
                               ((-180 20, -90 20, -90 -20, -180 -20, -180 20)))",
                       crs = crs(ersst))

ersst_trop.pac=crop(ersst, poly_trop_pac2)
temp_anom.troppac=terra::global(ersst_trop.pac,"mean",na.rm=T)$mean
Trop.Pac_tmpanom=tibble(dates=time(ersst),temp_anom=temp_anom.troppac)
Trop.Pac_tmpanom%>%ggplot(aes(x=dates,y=temp_anom))+
  geom_line()+ labs(x= "month", y="tempanomaly K",
                    title =" Tropical Pacific",subtitle = "Monthly & Area averaged")
library(itsmr)
colnames(Trop.Pac_tmpanom) #dates temp_anom
M.ssn=c("season",12)
Trop.Pac_trends=Trop.Pac_tmpanom%>%
  mutate(yr.anom=Resid(temp_anom,M.ssn),
         yr.smth=smooth.fft(yr.anom,0.1),
         trd1=trend(yr.anom,1),
         trd2=trend(yr.anom,2),
         trd3=trend(yr.anom,3),
         trd4=trend(yr.anom,4))

Trop.Pac_trends%>%ggplot(aes(x=dates,y=yr.smth))+geom_line()
Trop.Pac_trends%>%pivot_longer(cols=-dates,
                               names_to = "trend_poly",
                               values_to = "anomalies")%>%
                      subset(trend_poly %in% c("trd1","trd2","trd3","trd4"))%>%
                      ggplot(aes(x=dates,y=anomalies,col=trend_poly))+
                      geom_line()+labs(x="",y="K",title = "Tropical Pacific Polynomial Trend",
                                       subtitle = "Temperature Anomalies")
# max slope
max.slope=Trop.Pac_trends$trd4%>%diff()%>%max() #0.002472736
which.max(diff(Trop.Pac_trends$trd4)) # 1753
which.min(Trop.Pac_trends$trd4) # 359
T.est= (which.max(diff(Trop.Pac_trends$trd4)) -
          which.min(Trop.Pac_trends$trd4))*4 # 5576 month 465 years
T.est.2=(which.max(Trop.Pac_trends$trd4)-
           which.min(Trop.Pac_trends$trd4))*4 # 5580 month 465 years
# amplitude model sin estimate sd *sqrt(2)
amp.trop.pac=Trop.Pac_trends$yr.anom%>%sd()%>% {.*sqrt(2)} # 0.5623393
# estimate from diff
Trop.Pac_trends%>%dplyr::select(dates,trd4)%>%
  mutate(trdiff=c(NA,diff(trd4)))%>%
  ggplot(aes(x=dates,y=trdiff))+geom_line()
#  ======
which.max(Trop.Pac_trends$trd4) #1754
which.min(Trop.Pac_trends$trd4) #359
Trop.Pac_trends[359,] # 1909-11-15
t0=1754-359
Y.manual=tibble(dates=Trop.Pac_trends$dates,
             dtnm= 1:length(dates)-1,
             y= amp.trop.pac*sin((2*pi/T.est)*(dtnm-1752-50)))
Y.manual%>%ggplot(aes(x=dates,y=y))+geom_line(col=2)+
  geom_line(data = Trop.Pac_trends,aes(x=dates,y=trd2))
min1=which.min(Y.manual$y) #409
min2=which.min(Trop.Pac_trends$trd2) #410
const=(Trop.Pac_trends[[min2,6]]-
  Y.manual[[min1,3]])
Y.manual=Y.manual%>%mutate(y=y+const)
Y.manual%>%ggplot(aes(x=dates,y=y))+geom_line(col=2)+
  geom_line(data = Trop.Pac_trends,aes(x=dates,y=trd2))
