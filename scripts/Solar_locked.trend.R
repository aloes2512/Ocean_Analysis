df=readRDS("data/NOAA.OCEAN.ANOMALIES.rds")$data #2121 obs
library(tidyverse)
# median not used here
names(df)# dt.mnth sst.mean... anoma.mean

N=NROW(df) #2121
# eliminate annual and semi annual (at equator)
library(itsmr)
M.ssn=c("season",12,"season",6)
# combine solar_monthly with ocean anomaly
solar_mnthly=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))
Ocean_anomaly<-df%>% mutate(anomaly=Resid(anoma.mean,M.ssn))
# combine solar and ocean data dim 2120 4
Ocean_solar_anomaly=Ocean_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,anomaly,"Solar_Variation"=SI)
# anchor  at time-series solar power
library(splines)
# find zero crossings of Solar Irradiation SI
find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}
# nodes <- find_nodes(solar_mnthly$SI)
# time coordinates of zero-crossings
zero_crossing_times <- solar_mnthly$dt.mnth[find_nodes(solar_mnthly$SI)]
length(zero_crossing_times) #19
# 2. Fit a continuous spline locked @ exact nodes
# degree = 1 creates  sawtooth a continuous line that bends at the zero-crossings
# spline: degree = 3 creates a smooth, continuous curve through the zero-crossings
trd.solar <- lm(anomaly ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                data = Ocean_solar_anomaly)
Lckd.trd<-Ocean_solar_anomaly%>%mutate(
                 trd.sol=predict(trd.solar),
                 solar.res=anomaly-trd.sol)

plt.Baseline.trd<-Lckd.trd%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=solar.res),col="grey")+
  geom_line(aes(y=trd.sol),col=2)+
  labs(x=NULL,y="anomaly",title="Baseline Trend",
       subtitle="locked to solar irradiation\nresiduals(grey) mapping natural harmoniccs")
ggsave("figs/Baseline_trend.png")
# select absolute minimum;start and end 1950 trough; turning point

# Set standardized column names
col_names <- c("dt.mnth", "trd.sol")

# 1. Absolute minimum
min.pt <- Lckd.trd %>%
  arrange(trd.sol) %>%
  slice(1) %>%
  {.[, c(1, 4)]} # 1912 -0.278

# 2. Start 1950 Trough (Target value: ~0.0138)
strt.trgh <- Lckd.trd %>%
  subset(dt.mnth < 1950 & dt.mnth > 1940) %>%  # local window around 1944 peak/trough
  arrange(desc(trd.sol)) %>%
  slice(1) %>%
  {.[, c(1, 4)]} #1944.  0.0134

# end trough


# Save target height explicitly (or set target_val <- 0.0138 directly)
# 1. Start height of the trough (target value at t_start)
target_val <- strt.trgh[[2]] # 0.01378636
t_start    <- strt.trgh[[1]] # 1944.5


# 2. Extract rising segment post-1950
rising_segment <- Lckd.trd %>%
  subset(dt.mnth >= 1950 & dt.mnth <= 1985)

# 3. Exact interpolated time (t_end) where trd.sol reaches target_val
t_end <- approx(
  x = rising_segment$trd.sol,
  y = rising_segment$dt.mnth,
  xout = target_val
)
trd.end <-t_end$x # 0.01378
trgh_end<-tibble(x=unlist(t_end$y),
                 y=unlist(t_end$x)) # 1977.205
plt.Baseline.trd+
  geom_vline(xintercept = t_start,linetype = 2)+
  geom_point(data=strt.trgh,aes(x=dt.mnth,y=trd.sol))+
  geom_vline(xintercept = trgh_end$x,linetype = 2)+
  geom_point(data=strt.trgh,aes(x=dt.mnth,y=trd.sol))+
  geom_point(data=trgh_end,aes(x=x,y=y))

# 4. Exact trough duration in  years
trough_duration <- t_end$y - t_start #32.70477 years

# 3. End 1950 Trough (First point >= 1950 with trd.sol >= target_val)
end.trgh <- c(dt.mnth=t_end$y,trd.sol=t_end$x)

# 4. Turning Point
TP<-Ocean_solar_anomaly%>%mutate(
  trd.sol=predict(trd.solar))%>%
  subset(dt.mnth>=1975)%>%
  dplyr::select(dt.mnth,trd.sol,anomaly)%>%
  mutate(slp.trd=c(NA,diff(trd.sol)))%>%
  na.omit%>%arrange(slp.trd)%>%slice(1)%>%
  dplyr::select(dt.mnth,trd.sol)
plt.Baseline.trd+
  geom_vline(xintercept = t_start,linetype = 2)+
  geom_point(data=strt.trgh,aes(x=dt.mnth,y=trd.sol))+
  geom_vline(xintercept = trgh_end$x,linetype = 2)+
  geom_point(data=strt.trgh,aes(x=dt.mnth,y=trd.sol))+
  geom_point(data=trgh_end,aes(x=x,y=y))+
  geom_point(data=TP,aes(x=dt.mnth,y=trd.sol),shape= 1)

# Combine into one tibble
pts <- bind_rows(
  "Minimum"           = min.pt     %>% setNames(col_names),
  "Start 1950 Trough" = strt.trgh  %>% setNames(col_names),
  "End 1950 Trough"   = end.trgh   %>% setNames(col_names),
  "Turning Point"     = TP         %>% setNames(col_names),
  .id = "point_type"
)

library(ggrepel)
dt.mnth<-Lckd.trd$dt.mnth
plt.lckd.trd<-ggplot(Lckd.trd, aes(x = dt.mnth, y = trd.sol)) +
  # Background solar residual line
  geom_line(aes(y = solar.res), color = "gray70", alpha = 0.6) +

  # Fitted solar-locked trend line
  geom_line(color = "#F8766D", linewidth = 1) +

  # Highlight key points
  geom_point(data = pts, aes(x = dt.mnth, y = trd.sol), size = 2, color = "black") +

  # Add labels with repulsion to prevent line/point overlap
  geom_text_repel(
    data = pts,
    aes(x = dt.mnth, y = trd.sol, label = point_type),
    size = 2.5,

    box.padding = 0.5,
    point.padding = 0.3
  ) +

  labs(
    title = "Solar-locked Trend\n & Residual",
    subtitle = "extracted with spline\nfit to solar power",
    x = NULL,
    y = "anomaly"
  ) +
  theme_minimal()
ggsave("figs/SolLckd_trend.png")
# calculate aic
res=residuals(trd.solar)
n=length(res) #2120
rss=sum(res^2)
p=coef(trd.solar)%>%length() #23
k=p+1
aicc_solar <- n * log(rss / n) + 2 * k + n * (log(2 * pi) + 1) + (2 * k * (k + 1)) / (n - k - 1)
Ocean_solar_anomaly$trend.solar<-predict(trd.solar)
aicc_solar # -5010.75
#==============

#--------

#===========
# periods of trend.solar estimate reference secular

secul.mdl<-loess(anomaly ~ dt.mnth,span = 0.25,data =Ocean_solar_anomaly )
Lckd.trd<-Lckd.trd%>%mutate(secular=predict(secul.mdl))
Lckd.trd<-Lckd.trd%>%
  mutate(trend.res=trd.sol-secular)
# difference sol locked trend secular (smoothed) analyzed
my_signal=Lckd.trd$trd.sol-Lckd.trd$secular
FFT.har.trd=tibble(idx=1:N-1,
                   spc=fft(my_signal),
                   amp=Mod(spc))
idx.mx=FFT.har.trd%>%arrange(desc(amp))%>%head(20)%>%pull(idx)

FFT.har.trd%>%subset(idx<30)%>%
  ggplot(aes(x=idx))+
  geom_point(aes(y=amp))+
  labs(x="rel.frequency",
       title="Fourier Coef. Harmonic Trend")
id=idx.mx[1:14]+1 # 9  8 10 11 13  6  7
library(gsignal)
# 7 max harmonics
spc<-FFT.har.trd$spc%>%as.vector()
sel.mx<-which(Mod(spc)>5) # length 14
spc.mx<-rep(0+0i,N)
spc.mx[sel.mx]<-spc[sel.mx]
# reconstruct max 7 periods
har.part_trend=tibble(dt.mnth=Ocean_solar_anomaly$dt.mnth,
       y=Re(ifft(spc.mx)))%>%
  ggplot(aes(x=dt.mnth,y=y))+
  geom_line()
library(pracma)

y <- Re(ifft(spc.mx))
x <- Ocean_solar_anomaly$dt.mnth

# 1. Find local maxima (peaks)
p_max <- pracma::findpeaks(y)
max_x <- x[p_max[, 2]]
period.avrg=diff(max_x)%>%mean()%>%round() # 22
max_y <- p_max[, 1]          # Original positive values

# 2. Find local minima (troughs) by passing -y
p_min <- pracma::findpeaks(-y)
min_x <- x[p_min[, 2]]
min_y <- y[p_min[, 2]]       # Extract true original negative values using indices

# 3. Combine both
all_x <- c(max_x, min_x)
all_y <- c(max_y, min_y)
PKs<-tibble(x.loc=all_x,
            y.heigt=all_y)
har.part_trend=har.part_trend+
  geom_point(data=PKs,aes(x=x.loc,y=y.heigt),col=2)+
  labs(x=NULL,y="anomaly",title = "Fourier reconstr. Dominant Periods")


# diff=trd.solar-max.har.res)
plt.trd.res<-Lckd.trd%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trd.sol,col="trd.solar"))+
  geom_line(aes(y=trend.res,col="trd.res"))+
  geom_line(aes(y=secular,col="secular"))+
  labs(title="Solar-locked-trend",
       subtitle = "low-pass-filtered extracted periods")
plt.trd.res
#======== Rssa
## without selected external reference
library(Rssa)
my_trend<-Lckd.trd$trd.sol
N<-length(my_trend) #2120
L= round(N/3)# to reduce leakage and spectral resolution
s <- ssa(my_trend, L = L)
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
# reconstruct
# 1. Extract the reconstructed trend component vector
rc <- reconstruct(s, groups = list(Trend = 1:2,

                                   Har2=5:6,
                                   Har3=7:9))
rc.full<- reconstruct(s,groups = list(all=1:40))
#=======
# calculate periods
# Extract frequency, period, and damping factor for the harmonic pair (e.g., components 1 & 2)
p <- parestimate(s, groups = list(c(7:9)), method = "esprit")
# p(7:9) 242 month

plt.ssa.trd<-Lckd.trd %>% mutate(Trend=unlist(rc$Trend),
                    Har2=unlist(rc$Har2),
                    Har3=unlist(rc$Har3))%>%
  dplyr::select(-c(anomaly,
                   Solar_Variation,
                   solar.res,
                   secular,
                   trend.res))%>%
  pivot_longer(cols=-dt.mnth,
               names_to="compnts",
               values_to = "anomaly")%>%
  ggplot(aes(x=dt.mnth,y=anomaly,col=compnts))+
  geom_line()+
  labs(x=NULL,title="SSA filtered Trend",
       subtitle="compared with: dominant period(Har2)\n and solar locked trend")
ggsave("figs/SSA-Trend.png")

plt.comp.ssa_lckd<-Lckd.trd%>%
  mutate(full.trd=unlist(rc.full$all),
         Trend_vec = unlist(rc$Trend) )%>% #comp 1:2
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=full.trd),lwd=1.2)+
  geom_line(aes(y=trd.sol),col=2)+
  geom_line(aes(y=Trend_vec),col=3,lwd=1.5)+
  labs(subtitle = "compare: solar-locked-trend(black)\nssa-reconstructed(red)\nSSA-Trend_vector(comp 1:4;green) ")
ggsave("figs/trend_compare.png")
trend_vec<- unlist(rc$Trend)
har_vec2 <- unlist(rc$Har2) # comp 5:6
rc.anoma<-reconstruct(s, groups = list(Trend = 1:2,
                                       harmn=3:40))

Lckd.trd%>%
  mutate(full.trd=unlist(rc.full$all),
         Trend_vec = unlist(rc$Trend) )%>% #comp 1:4
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trd.sol),col=2,lwd= 1.1)+
  geom_line(aes(y=Trend_vec),col=3,lwd=1.2)+
  labs(x=NULL,y= "anomaly/trend",
       title = "Solar Locked Trend ",subtitle= "Trend: SSA components 1:4")


plt.harmnics=Lckd.trd%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.res),col=2)+
  geom_line(aes(y=unlist(rc.anoma$harmn)),col=4)+
  labs(x=NULL,y="Residuals",subtitle = "Difference: Solar-Locked- / Secular-Trend (red)\n SSA Harmonics (blue) ")
plt.harmnics
library(pracma)

y <- unlist(rc.anoma$harmn)
x <- dt.mnth

# 1. Find local maxima (peaks)
p_max <- pracma::findpeaks(y)
max_x <- x[p_max[, 2]]
period.avrg=diff(max_x)%>%mean()%>%round() # 23
max_y <- p_max[, 1]          # Original positive values

# 2. Find local minima (troughs) by passing -y
p_min <- findpeaks(-y)
min_x <- x[p_min[, 2]]
min_y <- y[p_min[, 2]]       # Extract true original negative values using indices

# 3. Combine both
all_x <- c(max_x, min_x)
all_y <- c(max_y, min_y)
PKs<-tibble(x.loc=all_x,
            y.heigt=all_y)

Lckd.trd%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=unlist(rc.anoma$harmn)),col=4)+
  geom_point(data=PKs,aes(x=x.loc,y=y.heigt),size = 2)+
  labs(x=NULL,y="anomaly",title="SSA reconstructed Harmonic",
       subtitle = paste0("quasi harmonic, estimated period: ",period.avrg," yrs"))
ggsave("figs/SSA_dominant_harmoni.png")

