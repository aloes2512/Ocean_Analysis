library(MASS)
library(tidyverse)

# Preparation: Round to integer for NegBinom requirements
SP_intpol <- SP_intpol %>%
  mutate(SP_round = round(SP.intpl))

# Fit global Negative Binomial (Estimates theta automatically)
# Here we fit it against time to see the secular trend
fit_nb <- glm.nb(SP_round ~ nmdt, data = SP_intpol)

# To compare to your 'k=60' GAM, we need to give glm.nb more flexibility.
# We can use a natural spline basis to mimic the GAM wiggles:
library(splines)
fit_nb_spline <- glm.nb(SP_round ~ ns(nmdt, df = 60), data = SP_intpol)

# Check the estimated theta (overdispersion)
# If this is very different from 1, your GAM family=negbin(1) might be mis-specified
theta_est <- fit_nb_spline$theta # 2.203
# another anchor code
library(MASS)
library(tidyverse)

# 1. Periods extracted from your GAM-FFT (converted to months)
gam_anchors_mo <- c(10.96, 10.42, 9.92, 11.57, 8.01, 12.25) * 12

# 2. Fit the 'Rigid' NegBinom GLM
# This model estimates the 'best' theta (overdispersion) automatically
# unlike the family = negbin(1) in your GAM.
harmonic_terms <- map_chr(gam_anchors_mo, ~{
  paste0("sin(2*pi*nmdt/", .x, ") + cos(2*pi*nmdt/", .x, ")")
}) %>% paste(collapse = " + ")

nb_rigid_formula <- as.formula(paste("SP_round ~ nmdt +", harmonic_terms))

fit_nb_rigid <- glm.nb(nb_rigid_formula, data = SP_intpol)

# 3. Compare the Theta
# If fit_nb_rigid$theta is much larger than 1, your GAM was under-dispersed
message("Estimated Theta: ", round(fit_nb_rigid$theta, 3))
# =========
# Comparison
library(tidyverse)
library(MASS)
library(mgcv)

# 1. Define your anchors from the GAM-FFT (Top 6 Schwabe + Long cycles)
# Assuming these are in years
solar_anchors_yr <- c(10.96, 10.42, 9.92, 11.57, 8.01, 12.25, 105, 210)
solar_anchors_mo <- solar_anchors_yr * 12

# 2. Build the Comparison Models
# Model A: Your existing GAM
gam_fit <- gam(SP.intpl ~ s(nmdt, k = 60),
               family = negbin(1),
               data = SP_intpol)

# Model B: Rigid GLM.NB using the anchors
# We use 'round' because glm.nb requires integer counts
harmonic_terms <- map_chr(solar_anchors_mo, ~{
  paste0("sin(2*pi*nmdt/", .x, ") + cos(2*pi*nmdt/", .x, ")")
}) %>% paste(collapse = " + ")

nb_formula <- as.formula(paste("round(SP.intpl) ~ nmdt +", harmonic_terms))
nb_fit <- glm.nb(nb_formula, data = SP_intpol)

# 3. Compare the "Solar Fingerprint"
comparison_metrics <- tibble(
  model = c("GAM (Local Spline)", "GLM.NB (Fixed Anchors)"),
  aic = c(AIC(gam_fit), AIC(nb_fit)),
  theta = c(1, nb_fit$theta), # GAM was forced to 1
  logLik = c(logLik(gam_fit), logLik(nb_fit))
)

print(comparison_metrics)
#==========
plot(gam_fit)
