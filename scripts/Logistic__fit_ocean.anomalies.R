# this file was corrected as mean had to be changed
Ocean_anoma=readRDS("ocean.averaged_monthly.temperatures.rds")
# compare Ocean_anoma with corrected
M.ssn=c("season",12)
Ocean_anoma=Ocean_anoma%>%mutate(yr.anom=Resid(ocean_anom,M.ssn),
                                 yr.smth=smooth.fft(yr.anom,0.1),
                                 res3=Resid(yr.smth,3),
                                 trd3=trend(yr.smth,3))
colnames(Ocean_anoma)

Ocean_anoma%>%ggplot(aes(x=dates,y=yr.smth))+geom_line()
dates=Ocean_anoma$dates
dates[1] # 1880-01-15
dates[1081] # 1970-01-15
Ocean_anoma[[1081,4]] # 0.0000318
as.numeric(max(Ocean_anoma$yr.smth)) # 0.003464049
which.max(Ocean_anoma$yr.smth) #1745
lngth=which.max(Ocean_anoma$yr.smth)-1081 # 664
L=((Ocean_anoma[which.max(Ocean_anoma$yr.smth),4]-Ocean_anoma[[1081,4]]))%>%as.numeric()
slp.strt=L/lngth %>%as.numeric() #5.169047e-06
x0= dates[1200+240]%>% as.numeric() # 10940
my_data= Ocean_anoma$yr.smth
# Define a simple logistic function
# y = L / (1 + exp(-k*(x-x0)))
# L = max value, k = steepness, x0 = midpoint year
#==========
# 1. Prepare a clean data frame take the offset into account
# Ensure dates and my_data are exactly the same length (1754)
# 1. Prepare Data
df_fit <- data.frame(
  y = as.numeric(Ocean_anoma$yr.smth),
  x = as.numeric(dates)
)

# 2. Estimate Starting Parameters
# Let's assume the curve starts at -2e-04 and rises to +2e-04
start_C <- -0.002          # The 'floor'
start_L <-  0.005          # The total 'rise' (bottom to top)
start_x0 <- 10940            # Midpoint (approx year 2000)
start_k  <- 0.001            # Daily growth rate

# 3. Fit the Symmetric Logistic
fit_logistic <- nls(
  y ~ C + (L / (1 + exp(-k * (x - x0)))),
  data = df_fit,
  start = list(C = start_C, L = start_L, k = start_k, x0 = start_x0),
  control = nls.control(maxiter = 1000, warnOnly = TRUE)
)

# 4. Extract Trend
logistic_trend <- predict(fit_logistic)

# 4. Get the trend
logistic_trend <- predict(fit_logistic)
AIC(fit_logistic) # -22429.81
length(logistic_trend) # 1754
tibble(dates=dates,logistic_trnd=predict(fit_logistic))%>%
  ggplot(aes(x=dates,y=logistic_trend))+geom_line()
res.logistic=Ocean_anoma$yr.smth-logistic_trend
res.poly3=Ocean_anoma$yr.smth-Ocean_anoma$trd3
mean(sum(res.logistic^2)) # 0.000285369
mean(sum(res.poly3^2)) # 0.0002643852
#================
#fit_physicist
# 1. Ensure x is centered to avoid large-number exponent issues
# Using years is best: 1600.0, 1600.08, etc.
df_full <- data.frame(
  y = as.numeric(Ocean_anoma$yr.smth),
  x = as.numeric(format(dates, "%Y")) + as.numeric(format(dates, "%j"))/366
)

# 2. The Model
# we use cbind() to tell plinear which terms have linear coefficients
# Term 1: 1 (intercept 'b')
# Term 2: x (slope 'm')
# Term 3: tanh(k*(x-x0)) (amplitude 'A')
fit_physicist <- nls(
  y ~ cbind(1, x, tanh(k * (x - x0))),
  data = df_full,
  start = list(k = 0.04, x0 = 1990), # Adjusted x0 based on your 'sharp rise' observation
  algorithm = "plinear",
  control = nls.control(maxiter = 1000, minFactor = 1/1024)
)
fit_physicist
AIC(fit_physicist) # -22800.16
# 3. Extract the trend and residuals
df_full$trend <- predict(fit_physicist)
df_full$resids <- df_full$y - df_full$trend
