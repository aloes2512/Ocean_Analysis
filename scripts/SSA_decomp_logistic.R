library(Rssa)
library(tidyverse)
library(itsmr)
Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
dt.mnth=Ocean_solar_anomaly$dt.mnth
#Ocean_solar_anomaly<-Ocean_solar_anomaly%>%mutate(low.=smooth.fft(anomaly,0.1))

#Ocean_solar_anomaly<-Ocean_solar_anomaly[-c(1:4),]
N<-NROW(Ocean_solar_anomaly) # 2120
trd_solar<-Ocean_solar_anomaly$trend.solar
Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trd_solar))
  #geom_line(aes(y=low.))
L <- 1058
s <- ssa(trd_solar, L = L)
# 1. Compute wcor
w <- wcor(s, groups = 1:20)
# 2. Extract matrix and convert to long-format tibble
w_mat <- as.matrix(w)

w_df <- expand.grid(F1 = 1:20, F2 = 1:20) %>%
  mutate(value = abs(as.vector(w_mat)))

# 3. Plot with clear axes
ggplot(w_df, aes(x = F1, y = F2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:20, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:20, expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "W-Correlation Matrix",
    x = "Eigenvector Index",
    y = "Eigenvector Index",
    fill = "W-Cor"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8),
    panel.grid = element_blank()
  )
plot(s,type="values")
plot(s,type="paired",idx=1:10)


# 1. Extract the reconstructed trend component vector
rc <- reconstruct(s, groups = list(Trend = 1:4,
                                   trd=1:2,
                                   Har2=5:6,
                                   Har3=7:8),
                                   Har4=9:10)
rc.full<- reconstruct(s,groups = list(all=1:40))

tibble(dt.mnth=dt.mnth,
  full.trd=unlist(rc.full$all),
  )%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=full.trd),lwd=1.5)+
  geom_line(aes(y=trd_solar),col=2)+
  labs(subtitle = "compare solar trend(black)\nssa-reconstructed(red) ")
Trend_vec <- unlist(rc$Trend) #comp 1:4
trend_vec<- unlist(rc$trd)
har_vec2 <- unlist(rc$Har2) # comp 5:6

#trend_vec<-trend_vec[-c(2113:2116)]
# 2. Build a data frame (replace time_vec with your actual time or index vector)
df_trend <- tibble(
  Time = Ocean_solar_anomaly$dt.mnth, # or your date/year axis: Ocean_solar_anomaly$time
  Trend = Trend_vec,
  trend= trend_vec,
  Har2=unlist(rc$Har2),
  Original = Ocean_solar_anomaly$trend.solar
)

plt.ssa.trend<-df_trend%>%pivot_longer(cols=-Time,names_to = "Compts",values_to = "anom")%>%
  ggplot(aes(x=Time,y=anom,col=Compts))+geom_line()+
  labs(x=NULL,y="Anomaly / Trend Units",
       title = "Solar-Locked Decomposed",
       subtitle= "ssa components Trend (trend); Hale(Har2) ")

Ocean_anoma.decomp=df_trend%>%dplyr::select("dt.mnth"=Time,
                                            "trend.ssa"=Trend,
                                            "har.ssa"=Har2,
                                            "trend.solar"=Original)


#===========
# fit to generalized logistic
library(tidyverse)
library(minpack.lm)
# 1.1 Fit standard 4-parameter logistic model to the 1-2 SSA trend
fit_log <- nlsLM(
  anomaly ~ L / (1 + exp(-k * (dt.mnth - t0))) + b,
  data = Ocean_solar_anomaly,
  start = list(L = 1.0, k = 0.03, t0 = 1990, b = -0.2)
)
round(coef(fit_log),3) #  L       k    t0         b
                        #  1.086 0.047 2015.991   -0.126
t0=coef(fit_log)[["t0"]] # 2015.991
#which.max(diff(predict(fit_log))) # 1992
TP=Ocean_solar_anomaly$dt.mnth[which.max(diff(predict(fit_log)))] # 2015.991
diff(predict(fit_log))[1992] # 0.001061991
Trend.upper.lim=coef(fit_log)["L"]+coef(fit_log)["b"] # 0.9594697
Trend.lower.lim=coef(fit_log)["b"] # -0.126381
Trend.range=coef(fit_log)["L"] # Trend.range
logs.fun=function(t,coefs){
  L=coefs["L"]
  k=coefs["k"]
  t0=coefs["t0"]
  b=coefs["b"]
  y= L / (1 + exp(-k * (t - t0))) + b
}
Logstc.fit<-tibble(dt.mnth=dt.mnth,
          y=predict(fit_log),
          anomaly=Ocean_solar_anomaly$anomaly)

Logstc.fit%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anomaly),col="grey")+
  geom_line(aes(y=y),col=2)+
  geom_hline(yintercept = Trend.lower.lim,linetype = 3)+
  geom_hline(yintercept = Trend.upper.lim,linetype = 3)+
  geom_vline(xintercept = TP,linetype = 4)+
  labs(x=NULL,title = "Growth Curve fitted to Anomaly",
       subtitle = "asymptots= horizontal dots;\n turning point= vertical line ")
# 2. Extract fitted industrial trend and compute residual

df_trend <- df_trend %>%
  mutate(
    industrial_logistic = predict(fit_log),
    har_trend    = Trend - industrial_logistic,
    anomaly=Ocean_solar_anomaly$anomaly,
    natural_residual = anomaly - industrial_logistic
  )
df_trend%>%ggplot(aes(x=Time))+
  geom_line(aes(y=predict(fit_log),col="indus"))+
  geom_line(aes(y=Trend,col="Trend"))+
  geom_line(aes(y=har_trend,col="har.trd"))
# 3. Visualize the separation
plt.preind.sep=ggplot(df_trend, aes(x = Time)) +
  geom_line(aes(y = Trend, color = "SSA Trend"), lwd = 1.1) +
  geom_line(aes(y = industrial_logistic, color = "Logistic "), lwd = 1, linetype = 5) +
  geom_line(aes(y = natural_residual, color = "Residual "), lwd = 0.6) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray50") +
  geom_vline(xintercept = t0,linetype = 2)+
  geom_hline(yintercept = Trend.lower.lim,linetype = 3)+
  geom_hline(yintercept = Trend.upper.lim,linetype = 3)+

  scale_color_manual(values = c(
    "SSA Trend" = "firebrick",
    "Logistic " = "black",
    "Residual " = "navy"
  )) +
  labs(
    title = "Solar Locked Trend Decomposed",subtitle="Industrial(Logistic) / Natural(Residual) ",
    x = "",
    y = "Anomaly",
    color = "Comp:"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
plt.preind.sep


ggsave("figs/Separ.preind.png")

