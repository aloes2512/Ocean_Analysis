Ocean_solar_anom<-readRDS("data/Ocean_solar_anom")
Ocean_solar_anom$Solar_Retained_Residuals <- Ocean_solar_anom$Anomaly - Ocean_solar_anom$Clean_Baseline_Trend

Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Solar_Retained_Residuals),col="grey")+
  geom_line(aes(y=Clean_Baseline_Trend),col=2)+
  labs(title = "Decomposition Ocean Anomalies",
       subtitle = "phaselocked to solarpower\n trend fitted to solar power")

# total harmonics
Ocean_decomp_anomalies<-Ocean_solar_anom%>%
  mutate(Anomaly.total=Solar_Retained_Residuals+har.trd)%>%
  dplyr::select(dt.mnth,Anomaly.total,trd3)
Ocean_decomp_anomalies%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Anomaly.total))+
  geom_smooth(aes(y=Anomaly.total))+
  labs(x="",title="Solar Locked Anomalies")
#-----
anoma.total<-Ocean_decomp_anomalies$Anomaly.total
dt.mnth=Ocean_decomp_anomalies$dt.mnth
last(dt.mnth) # 2026.25 Apr 2026
extend_signal <- function(x, n_back = 3000) {
  N <- length(x)
  N_ext=N+n_back
  t_ext= seq(from=-(n_back-1),N)

  idx_ext=((t_ext - 1) %% N) + 1
  x_ext <- x[idx_ext]
  return(x_ext)

}
sign_ext=extend_signal(anoma.total)
N=length(sign_ext) # 5116
Anoma.ext=tibble(
  dt.ext=seq(1600,last(dt.mnth),by=1/12),
  sign_ext=extend_signal(anoma.total)) %>%
  mutate(sign.long.hr=hr(sign_ext,N/1:20))
Anoma.ext%>%
  ggplot(aes(x=dt.ext))+
  geom_line(aes(y=sign_ext),col="grey")+
  geom_line(aes(y=sign.long.hr),col=2)
# residuals only
res.anomalies=Ocean_solar_anom$Solar_Retained_Residuals
Res.anoma.ext=tibble(
  dt.ext=seq(1600,last(dt.mnth),by=1/12),
  ext.anoma=extend_signal(res.anomalies))%>%
    mutate(long.ext.hr=hr(ext.anoma,N/3:15))
Res.anoma.ext%>%ggplot(aes(x=dt.ext))+
  #geom_line(aes(y=ext.anoma),col="grey")+
  geom_line(aes(y=long.ext.hr),col=2)
