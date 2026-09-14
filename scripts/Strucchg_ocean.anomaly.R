# use strucchange with Ocean Anomaly
Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
library(tidyverse)
library(itsmr)
library(strucchange)
# Eliminate seasons
M.ssn=c("season",12,"season",6)
Ocean_solar_anomaly<-Ocean_solar_anomaly%>%mutate(anomaly=Resid(Anomaly,M.ssn))
anomaly_df<-Ocean_solar_anomaly%>%dplyr::select(dt.mnth,anomaly)
# select h= 0.05 == 9 yrs and 0.1
mosum_anomaly <- efp(anomaly ~ 1, data = anomaly_df, type = "OLS-MOSUM", h = 0.05)
mosum_anomaly1 <- efp(anomaly ~ 1, data = anomaly_df, type = "OLS-MOSUM", h = 0.1)
str(mosum_anomaly)
str(mosum_anomaly1)
tms<-seq(0.0246,0.975,length.out=2012)
tms1<-seq(0.0496,0.95,lenght.out=1906)
# plot MOSUM
tibble(times=tms*2116/12+1850,
       values=mosum_anomaly$process)%>%
  ggplot(aes(x=times,y=values))+geom_line()
plot(mosum_anomaly, main = "MOSUM anomaly change process")
------
  tibble(times=tms1*2116/12+1850,
         values=mosum_anomaly1$process)%>%
  ggplot(aes(x=times,y=values))+geom_line()

#----
# 2. Second Pass: Pinpoint the exact breakpoint date
bp_anomaly <- breakpoints(anomaly~ 1, data = anomaly_df, breaks = 1)
bp_anomaly1 <- breakpoints(anomaly~ 1, data = anomaly_df, h=0.2)

# breakdate:0.7310964 percent
k=bp_anomaly$breakpoints # 1547
k1=bp_anomaly1$breakpoints
time_seq <- seq(1850, by = 1/12, length.out = 2116)
# Extract the exact breakpoint date
breakdate <- time_seq[k]
breakdate1 <-time_seq[k1]
# confidence interval
ci_bp <- confint(bp_anomaly)

# Observation indices for 2.5%, point estimate, 97.5%
k_ci <- ci_bp$confint
conf=time_seq[as.numeric(k_ci)]
# Map all three bounds directly to your time sequence
breakdate_ci <- time_seq[k_ci]
# Extract 2.5%, breakpoint, and 97.5% point estimates
conf <- time_seq[as.numeric(k_ci)]

# Draw point estimate and bounds separately
anomaly_df %>%   #subset(dt.mnth>1975&dt.mnth<1980)%>%

  ggplot(aes(x = dt.mnth, y = anomaly)) +
  geom_line(color = "black") +
  # Main breakpoint (solid red line)
  geom_vline(xintercept = conf[2], color = "firebrick", linewidth = 1) +
  # Confidence bounds (dashed blue lines)
  geom_vline(xintercept = conf[1], color = "navy", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = conf[3], color = "navy", linetype = "dashed", linewidth = 1.2)
#select and extend pre-ind
N<-NROW(anomaly_df) # 2116
k_full<-floor(N/12) # 176 years
gam_full  <- gam(anomaly ~ s(dt.mnth, k = k_full,  bs = "cc"), data = anomaly_df)

anomaly_train<-anomaly_df%>%subset(dt.mnth<1978.833)
N_pre.in=1978-1850 #128
k_train= 128

gam_train  <-gam(anomaly ~ s(dt.mnth, k = k_train,  bs = "cc"), data = anomaly_train)
train.seq<-anomaly_train$anomaly
pre.extended <- predict(gam_train, newdata = anomaly_df)
summary(pre.extended)
anomaly_df<-anomaly_df%>%mutate(pre.extended = predict(gam_train, newdata = anomaly_df),
                    anoma.filt=predict(gam_full),
                    indust=anoma.filt-pre.extended)
anomaly_df%>%ggplot(aes(x=dt.mnth)) +
  geom_line(aes(y=pre.extended),col=4)+
  geom_line(aes(y=anoma.filt),col=3)+
  geom_line(aes(y=indust),col=2)
# low pass periods > N/20 ~ 8.8 yrs
anomaly_df<-anomaly_df%>%mutate(pre.ind=hr(pre.extended,N/1:20),
                    obs.filt=hr(anoma.filt,N/1:20),
                    ind.filt=hr(indust,N/1:20))
anomaly_df%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=pre.ind),col=4)+
  geom_line(aes(y=obs.filt),col=3)+
  geom_line(aes(y=ind.filt),col=2)+
  labs(x="",y="ocean anomaly",
       title=" Ocean Industr Trends",
       subtitle="observed(green) minus \npre-industrial(blue) ")
