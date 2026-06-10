# from Area_trends.R
ersst_trop.pac=crop(ersst, poly_trop_pac)
temp_anom.trop.pac=terra::global(ersst_trop.pac,"mean",na.rm=T)$mean
Trop.Pac_tmp.anom=tibble(dates=time(ersst),temp_anom=temp_anom.trop.pac)
mean(Trop.Pac_tmp.anom$temp_anom) # 0.08
Trop.Pac=Trop.Pac_tmp.anom%>% mutate(trd1=trend(temp_anom,1),
                                     trd2=trend(temp_anom,2),
                                     trd3=trend(temp_anom,3),
                                     trd4= trend(temp_anom,4))
Trop.Pac%>%
  pivot_longer(cols=-c(dates,temp_anom),
               names_to = "poly.trend",
               values_to ="values" )%>%
  ggplot(aes(x=dates,y=values,col=poly.trend))+geom_line()
