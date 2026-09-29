#Trend methods
library(tidyverse)
library(itsmr)
source("scripts/Solar_locked.trend.R")
loc.lst=ls()
#1. sol locked trend
Ocean_trends<-Ocean_solar_anomaly%>%
  dplyr::select(dt.mnth,anomaly,trend.solar)
aicc_trd.solar<- -5010.735
rm(list=loc.lst)
N=NROW(Ocean_trends) #2120
#2. poly
anomaly=Ocean_trends$anomaly
dt.mnth=Ocean_trends$dt.mnth
poly.mdl=lm(anomaly ~ stats::poly(dt.mnth,7),data=Ocean_trends)
AIC(poly.mdl)#-4486.657
poly.mdl.=lm(anomaly ~ stats::poly(dt.mnth,27),data =Ocean_trends )
res_poly=residuals(poly.mdl.)
n<-length(res_poly) #2120
k<- 27+1
rss_poly=sum(res_poly^2) # 10.50885
aicc_poly<- n * log(rss_poly / n) + 2 * k + n * (log(2 * pi) + 1) + (2 * k * (k + 1)) / (n - k - 1)
aicc_poly # -5177.665
AIC(poly.mdl.)# -5176.442
Ocean_trends$trd7<-predict(poly.mdl)
Ocean_trends$trd27<-predict(poly.mdl.)
Ocean_trends%>%dplyr::select(-trd7)%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.solar,col="trd.sol"))+
  geom_line(aes(y=trd27,col="trd.ply"))
#3.low pass filtered trend periods> 108 month
Ocean_trends<-Ocean_trends%>%mutate(trd.low=smooth.fft(anomaly,f=0.0185))
Ocean_trends%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trd.low))
rss.low=sum((Ocean_trends$anomaly-Ocean_trends$trd.low)^2)
rss.low # 11.85
n<-length(Ocean_trends$trd.low) #2120
m= round(0.0185*n/2) #20
p = 2*m + 1 #{(sine + cosine for each bin, plus the mean/intercept)}
k = p +1 # +1 for variance
aicc_low=n * log(rss.low / n) + 2 * k + n * (log(2 * pi) + 1) + (2 * k * (k + 1)) / (n - k - 1)
# -4893.677
Ocean_trends%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trd27,col="poly"))+
  geom_line(aes(y=trd.low,col="trd.low"))+
  geom_line(aes(y=trend.solar,col="trd.sol"))
#4. trd ssa
L <- 1060
anomaly<-Ocean_trends$anomaly
s <- ssa(anomaly, L = L)
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
plot(s,type="paired",idx=1:10)
# #1. Extract the reconstructed trend component vector
rc <- reconstruct(s, groups = list(Trend = 1:4,
                                   trd=1:2,
                                   Har2=5:6,
                                   Har3=8:9))
rc.full<- reconstruct(s,groups = list(all=1:40))
# #2 plott
tibble(dt.mnth=dt.mnth,
       full.trd=unlist(rc.full$all),
)%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=full.trd),col="grey",lwd=1)+
  geom_line(aes(y=rc$Trend),col=2)+
  geom_line(aes(y=rc$Har2),col=3)+
  geom_line(aes(y=rc$Har3))
#4. loess

loes.mdl=loess(anomaly ~ dt.mnth,span = 0.25)
Ocean_trends$trd.loes<-predict(loes.mdl)
# combine all
Ocean_trends%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trd27,col="poly.trd"))+
  geom_line(aes(y=trd.low,col="low.trd"))+
  geom_line(aes(y=trend.solar,col="sol.trd"))+
  geom_line(aes(y=trd.loes,col="loes"))+
  geom_line(aes(y=rc$Trend,col="ssa.trd"),lwd=1.2,linetype= 3)+
  labs(x=NULL,y="Anomaly Trend",title= "Trend Extract.-Methods",
       subtitle = "SSA,Sol.-lcked,Low-pass,Polynomial")

