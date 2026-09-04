Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
N=NROW(Ocean_solar_anomaly) #2116
library(itsmr)
# 1. Fit seasonal residual
M.ssn <- c("season", 12, "season", 6)
df_full <- Ocean_solar_anomaly %>%
  mutate(anomaly = Resid(Anomaly, M.ssn)) %>%
  dplyr::select(dt.mnth, anomaly)
N <- nrow(df_full) # 2116
df_train <- df_full[1:floor(N/2), ] # Data pre-1938

# 2. Fit pre-industrial baseline on train data only
# Using a lower k to avoid overfitting noise in the training set
gam_train <- mgcv::gam(anomaly ~ s(dt.mnth, k = 30, bs = "cc"), data = df_train)
# 3. Calculate difference relative to the pre-industrial baseline
df_results <- df_full %>%
  mutate(
    extrapol_pre1938 = predict(gam_train, newdata = df_full),
    # Difference from the projected pre-industrial trend
    industrial_effect = anomaly - extrapol_pre1938 ,
    industr.smth= hr(industrial_effect,N/1:20)
  )
df_results%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=extrapol_pre1938,col="pre1938"))+
  geom_line(aes(y=industr.smth,col="ind.smth"))+
  labs(x="",subtitle= "anomaly split:pre 1938 @ industrial")
# 4. Fit 4-parameter logistic model to the industrial effect
fit_log <- nlsLM(
  industrial_effect ~ b + L / (1 + exp(-k * (dt.mnth - t0))),
  data = df_results,
  start = list(
    b = -0.15,  # Lower asymptote (allows negative shift)
    L = 0.80,   # Total vertical span / scale
    k = 0.05,   # Growth rate
    t0 = 1988   # Inflection point
  )
)
# Extract fitted values & parameters
df_results$logist_fit <- predict(fit_log)
params <- round(coef(fit_log),3)
# b        L        k       t0
#-0.003    0.632    0.106 1988.715
# 5. Plotting
t0=params[["t0"]]
lower_asymptote <- params["b"] # -0.003
upper_asymptote <- params["b"] + params["L"] # 0.629
df_results %>%
  ggplot(aes(x = dt.mnth)) +
  geom_line(aes(y = industrial_effect, color = "Observed"), alpha = 0.5) +
  geom_line(aes(y = logist_fit, color = "Logistic"), linewidth = 1.2) +
  geom_hline(yintercept = lower_asymptote, linetype = 2, color = "blue") +
  geom_hline(yintercept = upper_asymptote, linetype = 2, color = "red") +
  geom_vline(xintercept = params["t0"], linetype = 3) +
  scale_color_manual(values = c("Observed " = "grey50", "Logistic" = "cyan4")) +
  labs(
    x = "Year",
    y = "Anomaly Effect",
    title = "Industrial Model ",
    subtitle = paste0(" Turning-point (t0): ", round(params["t0"], 1))
    )
ggsave("figs/Industr.mdl.png")
#--------
colnames(df_results)
df_results <- df_results%>%
  mutate(
    signl = hr(anomaly,N/1:20),
    poly7 = trend(anomaly,7),
    difference = signl - extrapol_pre1938)

df_results%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=extrapol_pre1938,col="extrapol"))+
  geom_line(aes(y=difference,col="differnc"))+
  labs(x="",title = "Anomaly Pre 1938 @ Difference",
       subtitle="method:  extrapolated GAM ")
df_results%>%ggplot(aes(x=dt.mnth,y=signl))+
  geom_line(col="grey")+
  geom_line(aes(y=poly7))+
  geom_line(aes(y=logist_fit))+
  geom_vline(xintercept = 1997)+
  geom_vline(xintercept = t0)
#
# fit richards
## step 1 take  t0 to from logistic fit, and optim  with t0 fixed
t0_fixed= params[["t0"]] # estimated from logistic
B_fixed= params[["k"]]
fit_anomaly <- nls(
  anomaly ~ L + Delta / (1 + exp(-B_fixed * (dt.mnth - t0_fixed)))^(1 / nu),
  data =df_results,
  start = list(L = 0, Delta = 0.9, nu = 1)
)
coef1=coef(fit_anomaly)
## step2 t0 variable
L=coef1["L"] # -0.1247615
Delta=coef1["Delta"] # 0.52288
nu=coef1["nu"] # 1.035196
fit_anomaly_2<- nls(
  anomaly ~ L + Delta / (1 + exp(-B_fixed * (dt.mnth - t0)))^(1 / nu),
  data = df_results,
  start = list(t0=t0_fixed) #  1988
)
coef(fit_anomaly_2) # 1990.126
Rich_final=tibble(dt.mnth=df_full$dt.mnth,
  richrds=predict(fit_anomaly_2),
  sgnl=df_results$anomaly)
ceilng=L+Delta# 0.409
Rich_final%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=sgnl),col="grey")+
  geom_line(aes(y=richrds,col="Richard"))+
  geom_vline(xintercept = coef(fit_anomaly_2),
             linetype = 2)+
  geom_hline(yintercept = ceilng,linetype= 3)+
  labs(x="",y="trends",title="Anomalies fitted to Richards")
ceilng=L+Delta# ceilng 1.8
#------------
# 1. Fit initial logistic model
fit_log <- nlsLM(
  industrial_effect ~ b + L / (1 + exp(-k * (dt.mnth - t0))),
  data = df_results,
  start = list(b = -0.11, L = 0.85, k = 0.05, t0 = 1988)
)

log_coef <- coef(fit_log) #
B_fixed  <- as.numeric(log_coef["k"]) # Fixed forcing rate coefficient

# 2. Fit Richards model with FIXED B, leaving t0 free for ocean dynamics
fit_richards_fixedB <- nlsLM(
  industrial_effect ~ L + Delta / (1 + exp(-B_fixed * (dt.mnth - t0)))^(1 / nu),
  data = df_results,
  start = list(
    L     = as.numeric(log_coef["b"]),
    Delta = as.numeric(log_coef["L"]),
    t0    = as.numeric(log_coef["t0"]), # Free to shift due to climate modes
    nu    = 1.0
  ),
  lower = c(L = -0.5, Delta = 0.5, t0 = 1950, nu = 0.01), # Prevents asymptote collapse
  control = nls.lm.control(maxiter = 500)
)
# Extract parameters
p_B <- coef(fit_richards_fixedB) #
round(p_B,3)
#   L      Delta    t0       nu
# 0.001    0.611 1980.529    0.571
df_results$richardsB<-predict(fit_richards_fixedB )
# 3. Compare timing shift
t0_log     <- log_coef["t0"]
t0_richard <- p_B[["t0"]]
shift_years <- t0_richard - t0_log # -8.186014
#-----
# or gam smoothed
# Fit a stable GAM smooth over the raw industrial effect
gam_industrial <- mgcv::gam(industrial_effect ~ s(dt.mnth, k = 60, bs = "tp"), data = df_results)
df_results$industrial_gam <- predict(gam_industrial)
# or gam smoothed
# Fit a stable GAM smooth over the raw industrial effect
gam_industrial <- mgcv::gam(industrial_effect ~ s(dt.mnth, k = 60, bs = "tp"), data = df_results)
df_results$industrial_gam <- predict(gam_industrial)

ggplot(df_results, aes(x = dt.mnth)) +
  geom_line(aes(y = industrial_gam, color = "GAM Trend"), linewidth = 0.8) +
  geom_line(aes(y = richardsB, color = "Richards Model"), linewidth = 1.1) +
  geom_vline(xintercept = p_B["t0"], linetype = 2) +
  scale_color_manual(values = c("GAM Trend" = "grey40", "Richards Model" = "firebrick")) +
  labs(
    x = NULL, y = "Anomaly Trend",
    title = "Industrial Anomaly: GAM Trend vs. Richards Model"
  ) +
  theme_minimal()

ggplot(df_results, aes(x = dt.mnth)) +
  geom_line(aes(y = industrial_gam, color = "GAM Trd"), lwd = 0.6) +
  geom_line(aes(y = richardsB, color = "Richards"), lwd = 1.1) +
  geom_vline(xintercept = p_B["t0"], linetype = 2,col=3) +
  scale_color_manual(values = c("GAM Trd" = "grey40", "Richards" = "firebrick")) +
  labs(
    x = NULL, y = "Anomaly Trend",
    title = "Industrial Anomaly Model",
    subtitle="GAM Trend vs. Richards Model"
  ) +
  theme_minimal()

