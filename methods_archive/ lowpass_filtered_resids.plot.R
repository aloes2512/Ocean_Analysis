Ocean_solar_anom%>%
  mutate(smth.res.11=smooth.fft(Solar_Retained_Residuals,f=0.065),
         smth.res.3.4=smooth.fft(Solar_Retained_Residuals,f=0.02))%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Solar_Retained_Residuals),col="grey")+
  geom_line(aes(y=smth.res.11),col=2)+
  geom_line(aes(y=smth.res.3.4),col=4)+
  labs(x="",title = "Resids of sol.phaselocked-trend",
       subtitle="low-pass filtered  \ncutoff 11.2 years (red) 3.4 years (blue)")
ggsave("figs/Resids_phaselocked-trd.lowpass.png")
##