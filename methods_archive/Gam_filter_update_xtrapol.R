library(tidyverse)
library(mgcv)
df<-readRDS("data/NOAA.OCEAN.ANOMALIES.rds")$data
global.ts<-df%>%dplyr::select(dt.mnth,"Anomaly"=anoma.mean)
N=NROW(global.ts) # 2116
dt.mnth=global.ts$dt.mnth
# define N/2 as limit pre-industrial modern times
dt.preind.limit=global.ts$dt.mnth[floor(N/2)] # 1938.083

k_full<-180 # number of base functions fit cyclic with splines
gam_full  <- gam(Anomaly ~ s(dt.mnth, k = k_full,  bs = "cc"), data = global.ts)
require(broom)
GAM_fit=gam_full%>% augment()%>%
  dplyr::select(dt.mnth,Anomaly,"GAM180_fit"=.fitted,"glob_res"=.resid)
#-----
GAM_fit%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=glob_res,col="glob_res"))+
  geom_line(aes(y=GAM180_fit,col="GAM180_fit"))
GAM_fit|> # periods > 12 mnth
            ggplot(aes(x=dt.mnth)) +
            geom_line(aes(y=GAM180_fit),col=2)
# find nodes
find_nodes <- function(x) {
  # Find where the sign of the signal changes from negative to positive, or vice versa
  crossings <- which(diff(sign(x)) != 0) + 1
  return(crossings)
}

find_nodes(GAM_fit$GAM180_fit) #`10   25   47  180  199  325  351 1082 1113 1123 1154 1231 1238 1289
                              #1313 1333 1340 1426 1444 1470 1486 1522`
timespan=dt.mnth[c(10,1082)]# 1850.750 1940.083

# 1. Create a clean dataframe with your decimal dates and filtered data
# Assuming 'dt.mnth' is your decimal date and 'glob_res.smth' is your filtered data
GAM_fit=GAM_fit%>%mutate(time=dt.mnth-1850)%>%
  rename("signal"=GAM180_fit)
colnames(GAM_fit)
df<-GAM_fit%>%dplyr::select(time,signal)

# 2. Subset the data strictly to your justified preindustrial window (indices 16 to 1065)
df_preindustrial <- df[10:1082, ] ## 1850.750 1940.083
df_preindustrial%>%ggplot(aes(x=time))+geom_line(aes(y=signal))
# ======
library(gsignal)

# 1. Extract your perfectly bounded preindustrial signal (indices 16 to 1065)
# This segment has a length of N = 1050 months
preind_signal <- df_preindustrial$signal
N_base <- length(preind_signal) # 1073

# 2. Transform the baseline segment into the frequency domain
X_freq <- fft(preind_signal)

# 3. Determine how many months you need to fill total (from index 1 to the end of your data)
N_total <- length(dt.mnth) # 2116

# 4. Create the circular extension using a modulo index operation
# This seamlessly tiles the frequency-domain properties forward
circular_indices <- ((0:(N_total - 1)) %% N_base) + 1
extended_signal <- preind_signal[circular_indices]

# 5. Integrate back into your main dataframe
df <- data.frame(
  dt.mnth=dt.mnth,
  time = dt.mnth-1850,
  signal = df$signal,
  extended_baseline = extended_signal
)
df%>%ggplot(aes(x=dt.mnth,y=extended_baseline))+
  geom_line(aes(y=extended_baseline))+
  labs(x="",title="Extended Pre-industrial Anomalies")
# 6. Isolate the modern anomaly
df<-df%>%mutate(diff_anomaly=signal-extended_baseline)
df%>%ggplot(aes(x=time,y=extended_baseline))+
  geom_line(aes(y=extended_baseline,col="ext.preind"))+
  geom_line(aes(y=diff_anomaly,col="diff.anoma"))+
  labs(x="",title="Preindustrial Anomalies and\n Difference to Observed Anomalies",
       subtitle="data filtered with GAM 180 base functions")
#===============
# Fit with generalized logistic (look Wikipedia for definition)
## estimate M (t0) from 5th order polynomial fit
# 1. Fit the 5th-order polynomial over your full timeline
# Using I() ensures R treats the powers as raw mathematical exponents
poly_mod <- lm(diff_anomaly ~ time + I(time^2) + I(time^3) + I(time^4) + I(time^5), data = df)
poly_mod%>%augment()%>%dplyr::select(time,"poly5"=.fitted)%>%
  ggplot(aes(x=time+1850))+geom_line(aes(y=poly5))
# Extract the coefficients cleanly
coefs <- coef(poly_mod)
# 2. Predict onto a high-resolution grid to find the exact inflection point
# The inflection point of the modern surge is where the first derivative (slope)
# reaches its local maximum post-1950.
M_estimated <- dt.mnth[which.max(predict(poly_mod)%>%diff())]
print(paste("Estimated Inflection Point (M):", round(M_estimated,3)))
t.ctr=dt.mnth-1850
M.ctr=M_estimated-1850 # 146.333
df.ctr=tibble(t.ctr=dt.mnth-1850,
              diff_anomaly=df$diff_anomaly)
#"Y(t.ctr)= A0+A1*t.ctr+K*(1+exp(-B*(t.ctr-M.ctr)))^(-1/nu)"
#      A0         K           B           nu
#-0.003889104   0.635245515  0.120588599  1.592539744
fit_logistic <- nls(
  diff_anomaly  ~ A0+A1*t.ctr+K*(1+exp(-B*(t.ctr-M.ctr)))^(-1/nu),
  data = df.ctr,
  start=list(A0=-0.004,A1=-0.01,K=1.2,B=0.12,nu=1.6)
)
fitted.coef=coef(fit_logistic)
Y=tibble(dt.mnth=dt.mnth,
         y.fit=predict(fit_logistic ))
Y%>%ggplot(aes(x=dt.mnth,y=y.fit))+geom_line()+
  labs(x="",title = "Difference of Observed Anomalies and \n extended preidustrial anomalies",
       subtitle = "fitted with generalized logistic function")
# check stability of M.ctr
# Lock your optimized coefficients as constants
A0_opt <- fitted.coef["A0"] #0.0171524311
A1_opt <- fitted.coef["A1"] #-0.0004182251
K_opt  <- fitted.coef["K"]  #0.7089882090
B_opt  <- fitted.coef["B"] #0.1059179351
nu_opt <- fitted.coef["nu"] #1.5623638822

# Fit optimizing ONLY M.ctr
fit_M_stability <- nls(
  diff_anomaly ~ A0_opt + A1_opt * t.ctr + K_opt * (1 + exp(-B_opt * (t.ctr - M.ctr)))^(-1 / nu_opt),
  data = df.ctr,
  start = list(M.ctr = 143), # Use your exact polynomial value here
  control = nls.control(maxiter = 200)
)
coef(fit_M_stability) # 146.027 was estimated
library(zoo)
as.yearmon(coef(fit_M_stability)+1850)# Jan 1996 (Apr)
summary(fit_M_stability)
# 1. Extract the high-frequency modern deviations
df.ctr$richards_fit <- predict(fit_M_stability) # Or your chosen model matrix
df.ctr$modern_residuals <- df.ctr$diff_anomaly - df.ctr$richards_fit

# 2. Add these residuals to your extended preindustrial baseline
# This creates a "total corrected natural anomaly" series
df.ctr$ext.preind<-extended_signal
df.ctr$total_natural_baseline <- df.ctr$ext.preind + df.ctr$modern_residuals
df.ctr%>%ggplot(aes(x=dt.mnth))+geom_line(aes(y=total_natural_baseline))
# cut out a sequence staring and ending with zero
# 1. Run your existing function on the modern residuals
nodes <- find_nodes(df.ctr$total_natural_baseline)

# 2. Extract the first and last detected node indices
idx_start <- nodes[1] # 9
idx_end   <- nodes[length(nodes)] # 2101

# 3. Clip the residuals to this perfectly periodic window
clean_res_tile <- df.ctr$total_natural_baseline[idx_start:idx_end]
N_tile <- length(clean_res_tile) # 2093

# 4. Set up the deep past timeline (1600 to 1849)
hist_time <- seq(1600, 1849 + (11/12), by = 1/12)
N_hist <- length(hist_time)

# 5. Tile both the base signal and the clipped residuals backward
backward_base_idx <- ((-N_hist:-1) %% length(preind_signal)) + 1
backward_res_idx  <- ((-N_hist:-1) %% N_tile) + 1

#ALTERNATIVE :
# 1. Run the standard FFT on your node-to-node periodic modern residuals
N_tile <- length(clean_res_tile)
fft_res <- fft(clean_res_tile)

# 2. Extract amplitudes and phases for each frequency bin
# We scale by N_tile to normalize the coefficients back to physical temperature values
amplitude <- Mod(fft_res) / N_tile
phase     <- Arg(fft_res)

# 3. Create your historical time steps (measured relative to the original sequence length)
# Since your original data step size is 1 month (1/12 of a year):
hist_time <- seq(1600, 1849 + (11/12), by = 1/12)

# Create an index vector 't_hist' that counts backwards from the start of your tile
# If your modern tile starts at month index 'idx_start', the historical steps are:
t_hist <- seq(from = 1 - length(hist_time), to = 0, by = 1)

# 4. Reconstruct the signal by summing the Fourier components over the new timeline
# We only need the first half of the FFT vector due to symmetry (Nyquist limit)
half_N <- floor(N_tile / 2) # 1046

# Vectorized reconstruction loop
y_hist <- numeric(length(t_hist))
for (k in 1:half_N) {
  # Frequency for bin k
  omega <- 2 * pi * (k - 1) / N_tile

  # Accumulate the wave components (multiplying by 2 accounts for the negative frequencies)
  if (k == 1) {
    y_hist <- y_hist + amplitude[k] # DC component (mean value)
  } else {
    y_hist <- y_hist + 2 * amplitude[k] * cos(omega * t_hist + phase[k])
  }
}

# 5. Build your clean historical tibble
Y.back <- tibble(
  dts = hist_time,
  y_hist = y_hist
)

# Plot to verify
Y.back %>% ggplot(aes(x = dts)) + geom_line(aes(y = y_hist))
