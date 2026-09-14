Ocean_solar_anom<-readRDS("data/Ocean_solar_anom.rds")

# ext backwards Solar_locked_anomaly.R
library(gsignal)
extend_backwards <- function(x, n_back = 3000) {
  N <- length(x)
  N_ext=N+n_back
  t_ext= seq(from=-(n_back-1),N)

  idx_ext=((t_ext - 1) %% N) + 1

  x[idx_ext]
}

Long_solar.periods<-Ocean_solar_anom%>%
  dplyr::select(dt.mnth,Anomaly,trend.solar,trd5)
# add trd.harm to anomaly to get long.sum
#apply smth.res, trd.period, long.sum
N =length(Long_solar.periods$dt.mnth) #2116
Long_solar.periods=Long_solar.periods%>% mutate(trd.diff=trd5-trend.solar,
                                    long.Anomaly=hr(Anomaly,N/2:20),
                                    all.anomaly=trd.diff+long.Anomaly)
Long_solar.periods%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=long.Anomaly,col="long.Anoma"))+
  geom_line(aes(y=trend.solar,col="trend.solar" ))+
  geom_line(aes(y=trd5,col="poly5"))+
  geom_line(aes(y=trd.diff))+
  geom_line(aes(y=all.anomaly,col="all.anomaly"),linetype = 2)
n_back=3000
#trd.period=trd5-trend.solar==> trd.diff
#long.sum=trd.period+smth.res
# smth.res=smooth.fft(Anomaly,f=0.02)==>long.Aomaly
# all.anomaly=trd.diff+long.Anomaly
Long_dominant=Long_solar.periods%>%
  dplyr::select(dt.mnth,long.Anomaly,trd.diff,all.anomaly)
Long_dominant$dt.mnth[1] # 1850
# extend backwards all.anomaly ( base anomalies == Anomaly;+ trd.diff)
my_dominant=Long_dominant$all.anomaly
Ext.dominant= tibble(dates_ext= seq(1604,by=1/12,length.out=N+n_back),
                     anoma.ext=extend_backwards(my_dominant,n_back = 3000))

#============
Ext.dominant %>%
  ggplot(aes(x = dates_ext, y = anoma.ext)) +
  geom_line() +

  # 1. Background Rectangles (Dates must match column type)
  annotate("rect",
           xmin = 1645, xmax = 1715+11/12,
           ymin = -Inf, ymax = Inf, fill = "blue", alpha = 0.15) +
  annotate("rect",
           xmin = 1730, xmax = 1750+11/12,
           ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.15) +
  # NEW: Added Dalton Minimum Shading
  annotate("rect",
           xmin = 1790, xmax = 1830+11/12,
           ymin = -Inf, ymax = Inf, fill = "darkgreen", alpha = 0.12) +

  # 2. Reference Line for the end of LIA
  geom_vline(xintercept = 1850,
             linetype = "dashed", color = "darkred") +

  # 3. Text Labels
  annotate("text",
           x = 1680,
           y = 0,                       # Set to the middle of your y-axis scale
           label = "Maunder Minimum",
           angle = 90,
           vjust = 0.5,                 # Centers the text on the 'y' coordinate
           size = 3.0,
           color = "blue4",
           fontface = "bold") +

  annotate("text",
           x = 1740,
           y = 0,
           label = "18th C Warmth",
           angle = 90,
           vjust = 0.5,
           size = 3.0,
           color = "orange4",
           fontface = "bold") +

  # NEW: Added Dalton Minimum Text
  annotate("text",
           x = 1810,                    # Centered between 1790 and 1830
           y = 0,
           label = "Dalton Minimum",
           angle = 90,
           vjust = 0.5,
           size = 3.0,
           color = "darkgreen",
           fontface = "bold") +

  theme_minimal() +
  labs(title = "Ext. sum of Dominant Harmonics",
       subtitle = "harmonic part of phaselocked trend\n smoothed retained resids",
       x = "Year", y = "Dominant Harmonic")
