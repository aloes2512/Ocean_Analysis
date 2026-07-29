# Base knowledge on Solarforcing see Youtube
#browseURL("https://www.youtube.com/watch?v=BhNfOZ6YLbc")
library(tidyverse)
DF<-readRDS("data/NOAA.OCEAN.ANOMALIES.rds")$data
colnames(DF)
global.ts<-DF%>%dplyr::select(dt.mnth,"Anomaly"=anoma.mean)
N_global=NROW(global.ts) # 2116
dt.mnth=global.ts$dt.mnth
N_base=which(dt.mnth==1958) #1297
N_prewar=which(dt.mnth==1938) #1057
dt.mnth[1058] # 1938.083
library(itsmr)
#eliminate season
M.ssn=c("season",6,"season",12)
global.ts<-global.ts%>%mutate(anomaly=Resid(Anomaly,M.ssn))
# checked mean. e-19
# limit calculated with rollingcor
dt.preind.limit=1958.0# global.ts$dt.mnth[floor(N_global/2)] # 1938.083
dt.prewar.limit<-1938
global.ts  <- global.ts%>%mutate(anom.hr=hr(anomaly,N_global/1:20), # 9 yrs
                                 res.hr_anom=anomaly-anom.hr)

library(itsmr)
M.trd=c("trend",1)
preind_global.ts<- subset(global.ts,dt.mnth<=dt.preind.limit)%>%
  mutate(anomaly=Resid(Anomaly,M.trd))
prewar_global.ts<-subset(global.ts,dt.mnth<=dt.prewar.limit)%>%
  mutate(anomaly=Resid(Anomaly,M.trd))
prewar_global.ts%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anomaly))
global.ts%>%ggplot(aes(x=dt.mnth,y=anom.hr))+
  geom_line(aes(y=res.hr_anom,col="res.hr_anom"))+
  geom_line(aes(col="hr1_20"))+
  labs(x="",title="Low.pass Filtrd Anomalies(red)",
       subtitle = "residuals = diff(hr1_20,anomalies)")
#-----------
preind_global.ts%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anomaly,col="preind.anomaly"))+
  geom_line(aes(y=hr(anomaly,N_global/1:8),col="hr1_6"))+
  labs(title = "Normalized Preind.anomalies\n anomalies low-pass-filtered")
#--------------
prewar_global.ts%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anomaly,col="prewar.anomaly"))+
  geom_line(aes(y=hr(anomaly,N_global/1:8),col="anoma.hr1_6"))+
  labs(title = "Normalized Prewar.anomalies\n anomalies low-pass-filtered")

# extend prewar.anomaly 1850 to N_global=2026; N_prewar = 1057
circular_indices.1938 <- ((0:(N_global - 1)) %% N_prewar) + 1 #N_base= 1297; circular-length 2116
Prewar.ext=tibble(dt.mnth=global.ts$dt.mnth,
       ext.prewar=prewar_global.ts$anomaly[circular_indices.1938])%>%
       mutate(ext.prewar=ext.prewar-mean(ext.prewar),anoma.hr1_6=hr(ext.prewar,N_global/1:6))
Prewar.ext$ext.prewar%>%mean()#e-18
Prewar.ext%>%
  ggplot(aes(x=dt.mnth))+geom_line(aes(y=ext.prewar))+
  geom_line(aes(y=anoma.hr1_6,col="anoma.hr1_6"))+
  labs(x="",title="Pre WWII Anomaly Extended",
       subtitle = "reference for modern anomaly change")
#--------------------
# extend preind.anomaly 1850 to N_global=2026
circular_indices <- ((0:(N_global - 1)) %% N_base) + 1 #N_base= 1297; circular-length 2116
Preind.ext=tibble(dt.mnth=global.ts$dt.mnth,
                  ext.preind=preind_global.ts$anomaly[circular_indices])%>%
  mutate(ext.preind=ext.preind-mean(ext.preind),anoma.hr1_6=hr(ext.preind,N_global/1:6))
Preind.ext$ext.preind%>%mean()#e-18
Preind.ext%>%
  ggplot(aes(x=dt.mnth))+geom_line(aes(y=ext.preind))+
  geom_line(aes(y=anoma.hr1_6,col="anoma.hr1_6"))+
  labs(x="",title="Preindustrial Anomaly Extended",
       subtitle = "reference for modern anomaly change")

#===================
# backcast: extended to 1600
ende=global.ts$dt.mnth%>%last() # 2026.25
#extended timescale
dts_historic=seq(1600,ende,by=1/12)
N_hist=length(dts_historic) # 5116
circular_historic<-((1:N_hist-1))%%N_global+1 # range N_global 1850 to 2026.25

# center historic also to mean zero
mean(Preind.ext$ext.preind[circular_historic]) # --0.1180949
Anomaly.historic=tibble(dts_historic=seq(1600,ende,by=1/12),
                        anoma.hist=Preind.ext$ext.preind[circular_historic],
                        anomaly.historic=anoma.hist-mean(anoma.hist))
Anomaly.historic%>%ggplot(aes(x=dts_historic))+
  geom_line(aes(y=anomaly.historic,col="ext.preind"))+
  geom_line(aes(y=hr(anomaly.historic,N_global/1:6),col="anoma.hr1_6"))+
  labs(x="",title = "backcast of preind anomalies to 1600")#
#-------------
# alternative filter
global.ts_1.4  <- global.ts%>%mutate(anom.hr=hr(anomaly,N_global/1:4),
                                 res.hr_anom=anomaly-anom.hr)

global.ts_1.4%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anom.hr))+
  geom_line(aes(y=res.hr_anom))
Anomaly.historic.hr=tibble(x=dts_historic,
                        anoma.hist=Preind.ext$ext.preind[circular_historic],
                        anomaly.historic=anoma.hist-mean(anoma.hist),
                        longhr=hr(anomaly.historic,N_hist/5))

Anomaly.historic.hr%>%ggplot(aes(x=dts_historic))+
  geom_line(aes(y=anomaly.historic,col="ext.preind"))+
  geom_line(aes(y=longhr,col="Gleisberg"))+
  labs(x="",title = "Preind. Anomalies+ Gleisberg",
       subtitle = "anomalies 1850::1938 circular extended &\nmax.like fit : 1023.2 month period")
# Gemini :~1700: The trough aligns perfectly with the deep climax of the Maunder Minimum.
#         ~1840: The next major trough catches the Dalton Minimum beautifully.
#         In solar physics, 1749/1750 is a legendary boundary line.
          #When Rudolf Wolf back-calculated the historical sunspot numbers
          #to establish the standard numbering system we use today,
          #he chose 1749 as the start of Solar Cycle 0.
# superposition of several harmonics
Anomaly.historic.hr18=tibble(x=dts_historic,
                           anoma.hist=Preind.ext$ext.preind[circular_historic],
                           anomaly.historic=anoma.hist-mean(anoma.hist),
                           longhr=hr(anomaly.historic,N_hist/1:18))
Anomaly.historic.hr18%>%ggplot(aes(x=dts_historic))+
  geom_line(aes(y=anomaly.historic,col="ext.preind"))+
  geom_line(aes(y=longhr,col="18 harmonics"))
# include wolf
Anomaly.historic.hr37=tibble(x=dts_historic,
                             anoma.hist=Preind.ext$ext.preind[circular_historic],
                             anomaly.historic=anoma.hist-mean(anoma.hist),
                             longhr=hr(anomaly.historic,N_hist/1:37),
                             gleisberg=hr(anomaly.historic,N_hist/5))
Anomaly.historic.hr37%>%ggplot(aes(x=dts_historic))+
  geom_line(aes(y=anomaly.historic,col="ext.preind"))+
  geom_line(aes(y=longhr,col="37 harmonics"))+
  geom_line(aes(y=gleisberg,col="gleisberg"),linewidth =1.3)+
  labs(title = " Historic Anomalies Reconstructed",
       subtitle = "from pre-industrial anomalies 1850 to 1958")


#============
# calculate trend as difference

## difference observed - pre-idustrial
glob.data<-global.ts%>%left_join(Preind.ext,by="dt.mnth")%>%
  dplyr::select(dt.mnth,Anomaly,ext.preind)
# seasonality to remove from ext.preind
M=c("season",12,"season",6)
glob.data=glob.data%>%mutate(ext.preind=Resid(ext.preind,M),
                             difference=Anomaly-ext.preind,
                             trd5=trend(difference,5))

# turning point of trd5 =  max diff
#TP=glob.data$dt.mnth[which.max(glob.data$trd5%>%diff())] #1767
# set TP
TP<-1996
glob.data%>%ggplot(aes(x=dt.mnth))+geom_line(aes(y=trd5))+
  geom_vline(xintercept=TP,linetype= 2)
# signal to fit with Richards curve
t0.ctr=TP-1850 # 176.1667

# difference observed - pre idustrial
# formula
t.ctr=global.ts$dt.mnth-1850
M.ctr=TP-1850 # 146.333
df.ctr=tibble(t.ctr=t.ctr,
              diff_anomaly=glob.data$difference)
"Y(t.ctr)= A0+A1*t.ctr+K*(1+exp(-B*(t.ctr-M.ctr)))^(-1/nu)"
#      A0         K           B           nu
#-0.003889104   0.635245515  0.120588599  1.592539744
start=list(A0=-0.004,K=1.2,B=0.12,nu=1.6)
fit_logistic <- nls(
  diff_anomaly  ~ A0+A1*t.ctr+K*(1+exp(-B*(t.ctr-M.ctr)))^(-1/nu),
  data = df.ctr,
  start=list(A0=-0.004,A1=-0.01,K=1.2,B=0.12,nu=1.6)
)
coefs=coef(fit_logistic )
diff_anomaly.fit=coefs["A0"]+coefs["A1"]*t.ctr+coefs["K"]*(1+exp(-coefs["B"]*(t.ctr-M.ctr)))^(-1/coefs["nu"])
fit_1st=tibble(tme=t.ctr,y.fitted=predict(fit_logistic))
fit_1st%>%
  ggplot(aes(x=tme,y=y.fitted))+geom_line()+labs(title = "approx nls fit\nwith fixed M.ctr")
# 2nd. iteration
# Fit optimizing ONLY M.ctr

A0_opt=coefs["A0"]
A1_opt=coefs["A1"]
K_opt=coefs["K"]
B_opt=coefs["B"]
nu_opt=coefs["nu"]
fit_M_stability <- nls(
  diff_anomaly ~ A0_opt + A1_opt * t.ctr + K_opt * (1 + exp(-B_opt * (t.ctr - M.ctr)))^(-1 / nu_opt),
  data = df.ctr,
  start = list(M.ctr = 143), # Use your exact polynomial value here
  control = nls.control(maxiter = 200)
)
M.ctr=coef(fit_M_stability) # 145.9989
summary(fit_M_stability)
#-----
# 1. Extract the high-frequency modern deviations
df.ctr$richards_fit <- predict(fit_M_stability) # Or your chosen model matrix
# check stability of M.ctr
df.ctr$modern_residuals <- df.ctr$diff_anomaly - df.ctr$richards_fit
df.ctr%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=richards_fit))+
  geom_vline(xintercept=M.ctr+1850,linetype= 2)+
  labs(x="",title = "Difference Observed/\nPre-Industrial",subtitle="Generalized Logistic Fitted ")
#
Pre.industrials=mutate(global.ts,residuals_final=Anomaly-predict(fit_M_stability))
Pre.industrials%>%ggplot(aes(x=dt.mnth))+geom_line(aes(y=residuals_final),col="grey")
Pre.industrials=tibble(dt.mnth=global.ts$dt.mnth,
                  logistic.fit=predict(fit_M_stability),
                  anoma.preind=global.ts$Anomaly-logistic.fit,
                  anoma.hr1_6=hr(anoma.preind,N_global/1:6),
                  anoma.hr1_20=hr(anoma.preind,N_global/1:20),
                  mdl0=logistic.fit+anoma.hr1_6,
                  mdl20=logistic.fit+anoma.hr1_20)
#saveRDS(Preind.ext,"data/Preind.ext_logistic.trd.rds")

colnames(Pre.industrials)
preind.plot=Pre.industrials%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anoma.preind),col="grey")
preind.plot+  geom_line(aes(y=anoma.hr1_6,col="anoma.hr1_6"))+
  geom_line(aes(y=anoma.hr1_20,col="anoma.hr1_20"))+
  geom_line(aes(y=logistic.fit,col="logistic.fit"),linewidth = 1.5)

Pre.industrials%>%
  ggplot(aes(x=dt.mnth)) +
  geom_line(aes(y=logistic.fit,col="logistic.fit"),linewidth = 1.2)+
  geom_line(aes(y=mdl0,col="model_hr1_6"))+
  geom_line(aes(y=mdl20,col="model_hr1_20"))+
  labs(x="",title = "Ocean Anomaly Model",subtitle = "Trend Logistic Fit +\npre-industrial harmonics")
ggsave("figs/Preind_logistic.trd.png")
global_mdl=global.ts%>%
  mutate(model_logist.hr1_20=Pre.industrials$mdl20,
         logistic.fit=Pre.industrials$logistic.fit)
global_mdl%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Anomaly),col="grey")+
  geom_line(aes(y=model_logist.hr1_20,col="20 preind.har"))+
  geom_line(aes(y=logistic.fit,col="trd.logistic"))+
  labs(x="",title = "Model: Global Observed Anomalies",
       subtitle = "trd:Generalized Logistic Fit;\npre-industr: 20 harmonics model")
ggsave("figs/global.logist.20har.model.tiff")
# end of preindustrial?
library(itsmr)
observed_anomalies<-global.ts$Anomaly
harmonic_model_observed<-hr(observed_anomalies,N_global/1:20)
# 1. Calculate the pure Richards residuals (Observed - Pure S-Curve)
# (Assuming 'richards_fit' contains only the monotonic generalized logistic trend)
#from Anomalies_detrend_harmonic.R: fit_M_stability
richards_res <- observed_anomalies - predict(fit_M_stability)

# 2. Fit the 20 longest/most dominant periodic components using hr
# This will find the optimal continuous frequencies and fit amplitudes globally
harmonic_model_residuals <- hr(richards_res, N_global/1:20)


library(tidyverse)
library(slider)

# Your two parallel vectors
vector_obs <- harmonic_model_observed
vector_res <- harmonic_model_residuals

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
# Time> 1957 first min @ 1958.0
##browseURL("https://gml.noaa.gov/aggi/aggi.html")
#browseURL("https://de.wikipedia.org/wiki/Klimasensitivit%C3%A4t")
# "The Winter is coming" Atlantic Meridional Overturning Circulation (AMOC)
path="~/Desktop/Klima_Energiewende/Berichte_Veröffentlichungen/"
reports<-list.files(path)
reports[79] # "Winter Is Coming_ .......




