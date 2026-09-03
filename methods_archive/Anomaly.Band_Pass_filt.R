library(tidyverse)
library(minpack.lm)
library(Rssa)
library(itsmr)
Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")[-c(1:4),]
N<-NROW(Ocean_solar_anomaly) # 2112
L <- 132
trd_solar<-Ocean_solar_anomaly$trend.solar
s_sol <- ssa(trd_solar, L = L)
rc <- reconstruct(s_sol, groups = list(Trend = 1:2))
trend_vec <- as.numeric(rc$Trend)


df_trend <- tibble(
  Time = Ocean_solar_anomaly$dt.mnth, # or your date/year axis: Ocean_solar_anomaly$time
  Trend = trend_vec,
  Original = trd_solar
)

df_values <- data.frame(
  Component = 1:length(s_sol$sigma),
  SingularValue = s_sol$sigma,
  VarianceExplained = (s_sol$sigma^2) / sum(s_sol$sigma^2) * 100
)


ggplot(df_values[1:20, ], aes(x = Component, y = SingularValue)) +
  geom_point(color = "steelblue", size = 2) +
  geom_line(color = "steelblue") +
  scale_y_log10() +
  scale_x_continuous(breaks = 1:20) +
  labs(
    title = "SSA Singular Values (Log Scale)",
    x = "Component Index",
    y = "Singular Value (Sigma)"
  ) +
  theme_minimal(base_size = 11)

# plot cor matrix
w <- wcor(s_sol, groups = 1:10)
# 2. Extract matrix and convert to long-format tibble
w_mat <- as.matrix(w)

w_df <- expand.grid(F1 = 1:10, F2 = 1:10) %>%
  mutate(value = abs(as.vector(w_mat)))

# 3. Plot with clear axes
ggplot(w_df, aes(x = F1, y = F2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:10, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:10, expand = c(0, 0)) +
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
# reconstruct
res <- reconstruct(s_sol, groups = list(
  Secular_Trend = c(1),
  Multidecadal  = c(2, 3),
  Decadal       = c(4, 5)
))

res_df=tibble(dt.mnth=Ocean_solar_anomaly$dt.mnth,
              trd_solar=Ocean_solar_anomaly$trend.solar,
              Secular.trd=res$Secular_Trend,
              Multidecadal=res$Multidecadal,
              Decadal=res$Decadal)
time<-Ocean_solar_anomaly$dt.mnth
plt.ssa.decomp=res_df%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Secular.trd),col = "firebrick", lwd = 1.5)+
  geom_line(aes(y=trd_solar),col="grey",lwd=1)+
  geom_line(aes(y=Decadal+Multidecadal),col="navy")+
  labs(x="",y="anomaly",
       title="SSA decomposed Solar Trend",
       subtitle= "decomposed original @ harmonics")
print(plt.ssa.decomp)
# 1. Fit standard 4-parameter logistic model to the 1-2 SSA trend
L <- 64
s <- ssa(trd_solar, L = L)
rc <- reconstruct(s, groups = list(Trend = 1:2))
trend_vec <- as.numeric(rc$Trend)
df_trend <- tibble(
  Time = Ocean_solar_anomaly$dt.mnth, # or your date/year axis: Ocean_solar_anomaly$time
  Trend = trend_vec,
  Original = trd_solar
)

fit_log <- nlsLM(
  Trend ~ L / (1 + exp(-k * (Time - t0))) + b,
  data = df_trend,
  start = list(L = 1.0, k = 0.03, t0 = 1990, b = -0.2)
)
round(coef(fit_log),3) #  L       k    t0       b
                      # 1.082 0.047 2015.870 -0.126
t0=coef(fit_log)[["t0"]] # 2015.87
Trend.max=coef(fit_log)["L"]+coef(fit_log)["b"] # 0.956
Trend.min=coef(fit_log)["b"] # -0.1263
Trend.range=coef(fit_log)["L"] # 1.082458
# natural harmonics == pre industrial
M.ssn<-c("season",12,"season",6)
Ocean_solar_anomaly<-Ocean_solar_anomaly%>%
  mutate(anomaly=Resid(Anomaly,M.ssn),Secular.trd=predict(fit_log))
my_anomaly.resd<-Ocean_solar_anomaly%>%
  dplyr::select(dt.mnth,Secular.trd,anomaly,trend.solar)%>%
  mutate(Total.resd=anomaly-trend.solar,
         total.low_pass=hr(Total.resd,N/1:20))
my_anomaly.resd%>%ggplot(aes(x=dt.mnth,y=Total.resd))+
  geom_line(col="grey")+
  geom_smooth(method="loess",span=0.1)+
  geom_line(aes(y=total.low_pass),col=2)+
  labs(x="",title = "Pre-Industrial Harmonics",
       subtitle = "Low-Pass filtered: loess (blue),hr (red) ")

# separate by bandpass using Schwabe ~ 11.3  135.6 month
# Schwabe range 10.8 to 11.4 years
h.2schwabe<- c(10,12)*12 # 120: 144
h.hale<- c(19.8,20.5)*12 # 237.6 246.0
h.gleis<- c(80,90)*12 # 960 1056
h.sum<- c(0.1,21)*12
my_anomaly.resd<-Ocean_solar_anomaly%>%
  dplyr::select(dt.mnth,Secular.trd,anomaly,trend.solar)%>%
  mutate(Total.resd=anomaly-trend.solar,
         schwabe.low_pass=hr(Total.resd,h.2schwabe),
         hale.low_pass=hr(Total.resd,h.hale),
         gleis.low_pass=hr(Total.resd,h.gleis),
         sum.low_pas=hr(Total.resd,h.sum))


my_anomaly.resd%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Total.resd),col="grey")+
  geom_line(aes(y=schwabe.low_pass,col="schwb"))+
  geom_line(aes(y=hale.low_pass,col="hale"))+
  geom_line(aes(y=gleis.low_pass,col="gleis"))+
  geom_line(aes(y=sum.low_pas,col="sum"))+
  labs(x="",title = "Total Residuals (pre-industrial)",
       subtitle = "harmonic parts (Schwabe, Hale)")
# bandpass harmonics only
my_anomaly.resd%>%ggplot(aes(x=dt.mnth))+
  #geom_line(aes(y=Total.resd),col="grey")+
  geom_line(aes(y=schwabe.low_pass,col="schwb"))+
  geom_line(aes(y=hale.low_pass,col="hale"))+
  #geom_line(aes(y=gleis.low_pass,col="gleis"))+
  geom_line(aes(y=sum.low_pas,col="sum"),lwd= 1.2)+
  labs(x="",title = "Total Residuals (pre-industrial)",
       subtitle = "harmonic parts (Schwabe, Hale)")
# is the schwabe plot the consequence of an
#interaction of two frequencies
library(lomb)

# Compute high-density periodogram strictly across the Schwabe window (10 to 12 years)
# Output periods evaluated at micro-steps without zero-padding
lsp_schwb <- lsp(
  my_anomaly.resd$Total.resd,
  times = my_anomaly.resd$dt.mnth,
  from = 10, to = 12,           # Target search window in years
  ofac = 20,                    # Oversampling factor (20x denser than standard FFT)
  type = "period"
)
str(lsp_schwb)
# Extract peak periods T1 and T2 directly from the oversampled spectrum
library(gsignal)
length(lsp_schwb$power) # 58
PEAK_df=tibble(period=seq(10,12,length.out=58),
               power=lsp_schwb$power)
period=seq(10,12,length.out=58)
Peaks=findpeaks(lsp_schwb$power,DoubleSided = T)
PEAK_df%>%ggplot(aes(x=period,y=power))+geom_line()+
  geom_vline(xintercept = period[17],linetype = 2)+
  geom_vline(xintercept = period[45],linetype = 2)
period[Peaks$loc%>%as.numeric()]#10.56140 11.08772 11.54386 11.96491
# fit linear
# Fixed sideband frequencies derived from periodogram peaks
T1 <- 10.56
T2 <- 11.54

# Direct linear regression model for Schwabe residual dynamics
fit_schwabe <- lm(
  Total.resd ~ sin(2*pi/T1 * dt.mnth) + cos(2*pi/T1 * dt.mnth) +
    sin(2*pi/T2 * dt.mnth) + cos(2*pi/T2 * dt.mnth),
  data = my_anomaly.resd
)

# Extract reconstructed wave
my_anomaly.resd$schwabe_reconstructed <- predict(fit_schwabe)
my_anomaly.resd%>%ggplot(aes(x=dt.mnth,y=schwabe_reconstructed))+
  geom_line()

