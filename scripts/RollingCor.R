Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
N_train<-NROW(Ocean_solar_anomaly)/2
df_train<-Ocean_solar_anomaly[1:N_train,]
df_full <- Ocean_solar_anomaly
gam_model <- gam(Anomaly ~ s(dt.mnth, k = 180, bs = "cc"), data = Ocean_solar_anomaly)

gam_train <- gam(Anomaly ~ s(dt.mnth, k = 90, bs = "cc"), data = df_train)

df_full$ext.preind<-predict(gam_train, newdata = df_full)
df_full<-df_full%>%mutate(trend.ind=Anomaly-ext.preind)
DF=df_full%>%dplyr::select(dt.mnth,Anomaly,ext.preind,trend.ind)

# Rolling_corr
library(slider)
# Your two parallel vectors
vector_obs <- DF$Anomaly
vector_res <- DF$ext.preind
dt.mnth<-DF$dt.mnth
# Calculate the rolling correlation over a 132-month (11-year) window
# 'slide_index' passes a slice of the indices to your function
rolling_corr <- slide_index_dbl(
  .x = 1:length(vector_obs),
  .i = 1:length(vector_obs),
  .f = function(idx) {
    # If the window doesn't have a full 11 years of data, return NA
    if(length(idx) < 132) return(NA)

    # Calculate simple vector correlation on the current slice of indices
    cor(vector_obs[idx], vector_res[idx])
  },
  .before = 131,   # Looks 131 steps backward + current step = 132 months
  .complete = TRUE # Automatically returns NA for incomplete windows at the start
)

RollingCor=tibble(Time=dt.mnth,rolling_corr=as.numeric(rolling_corr))
plt.cor=RollingCor%>%ggplot(aes(x=Time,y=rolling_corr))+
  geom_line(col="blue")+
  labs(y="11-Year Rolling Correlation",
       title="Correlation of harmonic fitted observations\n  and logistic detrended residuals")
print(plt.cor)
rolling_corr.1=as.numeric(rolling_corr)[1:2000]
corlim<-which(rolling_corr.1>0.5)%>%last()
dt.mnth[corlim]# 1945
