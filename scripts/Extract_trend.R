# Trend separation methods
# data derived from "Solar_locked_anomaly.R"
library(tidyverse)
library(itsmr)
library(mgcv)
# 1. split by period
Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
N=NROW(Ocean_solar_anomaly)
Ocean_anomaly<-Ocean_solar_anomaly%>%
  dplyr::select(dt.mnth,Anomaly)
Ocean_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Anomaly))+
  labs(x="",title= "Global Ocean Temperature Anomalies",
       subtitle = "anomalies defined as SST -mean(SST)",
       caption = "NOAA monthly interpolated SST")
# select periods > 105 month ~ 8.8 years
Ocean_anomaly<-Ocean_anomaly%>%
              mutate(lng.perds.9=hr(Anomaly,N/1:20),
                     res.9=Anomaly-lng.perds.9,
                     lng.perds.12=hr(Anomaly,N/1:12))

Ocean_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=lng.perds.9,col="perds>9"))+
  geom_line(aes(y=res.9,col="frqu>0.01"))+
  labs(x="",title = "Splitted Anomaly",
       subtitle = "periods > 9 years & frqu > 0.01 [1/mnth]")
# 2. Solar phase locked & secular
plt.trds_sollck.lowp<-Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.solar,col="sol.lckd"))+
  geom_line(data=Ocean_anomaly,aes(x=dt.mnth,y=lng.perds.12,col="lowP.filt"))+
  labs(x="",title = "Anomaly Trends",
       subtitle= "low-pass filter/ solar-phase locked splines")
# 2.1 Decompose solar locked trend by secular trend
Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
# plt trd.solar
Ocean_solar_anomaly<-Ocean_solar_anomaly%>%rename("trd.solar"=trend.solar)
plt.trend.sol<-Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trd.solar,col="trd.solar"))
# 2.2 secular trend.solar span=1
Ocean_solar_anomaly<-Ocean_solar_anomaly%>%mutate(secular.trd=predict(
  loess(trd.solar ~ dt.mnth, data = Ocean_solar_anomaly, span = 1)
))
plt.trd.secul=plt.trend.sol+
  geom_line(data=Ocean_solar_anomaly,
            aes(y=secular.trd,col="secul.trd"))
# 2.3.periodic_harmonics (trd.solar - secular.trd)
Ocean_solar_anomaly<-Ocean_solar_anomaly%>%
  mutate(periodic.trd=trd.solar-secular.trd)
plt.trd.decomp=plt.trd.secul+
  geom_line(data=Ocean_solar_anomaly ,
            aes(y=periodic.trd,col="per.trd"))+
  labs(x="",y="trends [K]",title = "Anomaly Long Periods",
       subtitle="per.trd=solar.locked-secular")

# 3. GAM filtered
Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")

N<-NROW(Ocean_solar_anomaly)
#  3.1
gam_model <- gam(Anomaly ~ s(dt.mnth, k = 180, bs = "cc"), data = Ocean_solar_anomaly)
Ocean_solar_anomaly$gamtrd180<-predict(gam_model)
Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=gamtrd180,col="gam180trd"))+
  geom_line(aes(y=trend.solar,col="trend.solar"))+
  labs(x="",
       title = "Trend-Comparison: GAM spline vs\n Solar-Power Phase-locked Trend")
# 3.2 train with half obs ("Solar_anomalies_extended.R")
N_train=floor(N/2) # 1058
df_train<-Ocean_solar_anomaly[1:N_train,]
df_full <- Ocean_solar_anomaly
# 2. GAM fit  cyclic  Splines (bs = "cc")

# k = 90 fits the  1905er trough perfect : cyclical smooth (bs = "cc"),
# fitting with Penalized Maximum Likelihood Estimation 90 parameters say beta 1 t0 beta 90
gam_train <- gam(Anomaly ~ s(dt.mnth, k = 90, bs = "cc"), data = df_train)
df_full$ext.preind<-predict(gam_train, newdata = df_full)
df_full<-df_full%>%mutate(trend.ind=Anomaly-ext.preind)
df_full%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=ext.preind,col="ext.preind"))+
  geom_line(aes(y=Anomaly,col="Anomaly"))+
  geom_line(aes(y=trend.ind,col="trend.ind"))+
  geom_vline(xintercept = df_full$dt.mnth[N_train],
             linetype = 2,col=2)+
  ggtitle("Splitted Anomaly: Train / Extd.-train\nPreind defined: before 1938")
# 2.2 Low pass N/20 limited Anomaly & Split Ind./ Preind
DF=df_full%>%dplyr::select(dt.mnth,Anomaly,ext.preind,trend.ind)
DF=DF%>%mutate(Anoma.9=hr(Anomaly,N/1:20),
               Preind.9=hr(ext.preind,N/1:20),
               Ind.trnd.9=hr(trend.ind,N/1:20))
DF%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Anoma.9,col="Anomaly"))+
  geom_line(aes(y=Preind.9,col="Preind."))+
  geom_line(aes(y=Ind.trnd.9,col="Indust."))+
  labs(x="",title = "Decomposed Anomaly",
       subtitle= "Preindustrial & Deviation= Difference")
# poly is masked by gsignal::poly and needs to be specified
# --------
df_full<-df_full%>%
  mutate(poly5=predict(lm(df_full$trend.ind~ stats::poly(df_full$dt.mnth,5))))
TP5=which.max(diff(df_full$poly5)) #1766 1765
plt.GAM_diff.ext=df_full%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.ind,col="trend.ind"))+
  geom_line(aes(y=poly5,col="poly5"))+
  geom_vline(xintercept = 1850+TP5/12,linetype = 2)+
  labs(x="",title= "GAM filtered: Anomaly \n Pre-WWII_extended",
       subtitle = "trend.ind =Anomaly-ext.preind")

df_full$dt.mnth[1766] # 1997.083
M.ctr<-1997-1850 # 147
df.ctr<-df_full%>%mutate(t.ctr=dt.mnth-1850)%>%
  dplyr::select(t.ctr,trend.ind)

# fit logistic
fit_logistic <- nls(
  trend.ind ~ A0+K*(1+exp(-B*(t.ctr-M.ctr)))^(-1/nu),
  data = df.ctr,
  start=list(A0=-0.004,K=0.8,B=0.12,nu=1.6),
  control = nls.control(maxiter = 200)
)
coefs5=coef(fit_logistic)
# 2nd iteration
A0=coefs5["A0"] #-0.005944196
K=coefs5["K"]   # 0.6627517
B=coefs5["B"]   # 0.109869
nu=coefs5["nu"] # 1.544103
#2nd iteration
fit_M_stability <- nls(
  trend.ind ~ A0+K*(1+exp(-B*(t.ctr-M.ctr)))^(-1/nu),
  data = df.ctr,
  start=list(M.ctr= 146),
  control = nls.control(maxiter = 200)
)
coef(fit_M_stability) # 146.5237
M_final=1850+146.5 # 1996.5
df_full$y.logistic=predict(fit_M_stability)
df_full%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.ind,col="trend.ind"),alpha=0.5)+
  geom_line(aes(y=poly5,col="poly5"))+
  geom_line(aes(y=y.logistic,col="y.logistic"))+
  labs(x="",title= "Anomaly-Trend Industrial-times",
       subtitle="fitted with polynomial 5th°;gerealized logistic")
Ocean_solar_anomaly$logistic.trd<-df_full$y.logistic
# 3.2 extended preind
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
# new train
df_newtrain<-Ocean_solar_anomaly[1:corlim,]
df_full <- Ocean_solar_anomaly
#adapt k to shorter ts
k.new=floor(180/2116*corlim)
gam_newtrain <- gam(Anomaly ~ s(dt.mnth, k = k.new, bs = "cc"), data = df_newtrain)
df_full$new.preind<-predict(gam_train, newdata = df_full)
df_full<-df_full%>%mutate(trend.ind=Anomaly-new.preind)
df_full%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=new.preind,col="new.preind"))+
  geom_line(aes(y=Anomaly,col="Anomaly"))+
  geom_line(aes(y=trend.ind,col="trend.ind"))+
  geom_vline(xintercept = df_full$dt.mnth[corlim],
             linetype = 2,col=2)+
  ggtitle("Splitted Anomaly: NewTrain\nnew.preind defined: before 1945")
# visualize
DF=df_full%>%dplyr::select(dt.mnth,Anomaly,new.preind,trend.ind)
DF=DF%>%mutate(Anoma.9=hr(Anomaly,N/1:20),
               Preind.9=hr(new.preind,N/1:20),
               Ind.trnd.9=hr(trend.ind,N/1:20))
DF%>%ggplot(aes(x=dt.mnth))+
  #geom_line(aes(y=Anoma.9,col="Anomaly"))+
  geom_line(aes(y=Preind.9,col="Preind."))+
  geom_line(aes(y=Ind.trnd.9,col="Indust."))+
  geom_vline(xintercept = df_full$dt.mnth[corlim],col=2,linetype = 2)+
  labs(x="",title = "Decomposed Anomaly",
       subtitle= "Preindustrial & Deviation= Industrial")



# compare the 4 methods
#Ocean_solar_anomaly<-Ocean_solar_anomaly%>%rename("trend.solar"=Baseline_Trend)
plt.comp4=Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trd5,col="trd5"))+
  geom_line(aes(y=trend.solar,col="trend.solar"))+
  geom_line(aes(y=gamtrd180,col="gamtrd180"))+
  geom_line(data=Ocean_anomaly,aes(x=dt.mnth,y=lng.perds.12,col="lowP.filt"))+
  labs(title = "Compare 4 Methods: Extract Trend")
Ocean_solar_anomaly$low.pass<-Ocean_anomaly$lng.perds.12
saveRDS(Ocean_solar_anomaly,"data/Ocean_solar_anomaly.rds")
# model industrial trend.ind+mean(df_full$new.preind)
Ocean_solar_anomaly<-Ocean_solar_anomaly%>%
  mutate(ind.model=logistic.trd+mean(df_full$new.preind),
         gam.lowpas=hr(gamtrd180,N/1:8))
plt.trd=Ocean_solar_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=ind.model,col="Logis.Mdl"),linewidth = 1.2)+
  #geom_line(aes(y=low.pass),col=2)+
  geom_line(aes(y=gam.lowpas,col="low.ps_obs"))+
  geom_line(aes(y=low.pass-ind.model,col="difference"),linewidth = 1.2)+
  geom_line(aes(y=trd7,col="poly.trd7"))+
  labs(x="",title = "Industrial Times Logistic Model\nlong periods residual harmonic")


