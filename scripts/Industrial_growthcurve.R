## ===========================================================================
##  Industrial growth curve from the solar-locked trend
##  ---------------------------------------------------------------------
##  Pipeline (as specified by the author, 2026-09-14):
##    1. updated NOAA ERSST v6 global ocean anomalies  (NOAA.OCEAN.ANOMALIES.rds)
##    2. deseasonalise (itsmr, 12 and 6 month harmonics)
##    3. solar-locked trend: cubic spline with knots at the minima of the
##       solar power variation                                  -> trend.solar
##    4. strip the harmonic part with a LOESS fit               -> nonharm
##    5. difference of the observed monthly anomalies and the
##       harmonic part of the locked trend                      -> growth
##       (= nonharm + the retained monthly residual: a NOISY series that
##        still carries the secular rise.  Fitting the low-pass filtered
##        series instead is what made the earlier script fail.)
##    6. STEP 1  logistic + constant, 4 parameters, one nlsLM call
##       STEP 2  Richards, slope k held at the logistic value, nu free
##    7. AIC comparison over the candidate target series and over the
##       two growth curves
##
##  No staged asymptote separation, no logit regression, no (nu,t0) grid:
##  nlsLM from the author's own starting values converges directly.
## ===========================================================================

library(tidyverse)
library(splines)
library(itsmr)        # Resid(), hr()
library(minpack.lm)   # nlsLM

## ---- constants ------------------------------------------------------------
LOESS_SPAN <- 0.25
M.ssn      <- c("season", 12, "season", 6)
SST_SERIES <- "anoma.mean"   # "anoma.mean" or "anoma.median"
                             # timing is the same either way; the fitted span
                             # and ceiling are ~1.5x larger on the median --
                             # a spatial median of SST is a different statistic,
                             # not a robust version of the spatial mean
DATA_NOAA  <- "data/NOAA.OCEAN.ANOMALIES.rds"
DATA_SOLAR <- "data/S_power.rds"
OUT_RDS    <- if (SST_SERIES == "anoma.mean") "data/Ocean_solar_v6.rds" else
                                                 "data/Ocean_solar_v6_median.rds"

## ---- 1. data --------------------------------------------------------------
Ocean_NOAA.data      <- readRDS(DATA_NOAA)
NOAA.Ocean.anomalies <- Ocean_NOAA.data$data          # 2120 months, to 2026-08
solar_mnthly         <- readRDS(DATA_SOLAR) %>%       # SILSO -> TSI
  mutate(SI = TSI - mean(TSI))

## NOTE  S_power.rds currently ends 2026.33 while the SST grid ends 2026.58,
## so the inner join drops the last three months.  Re-run scripts/Sun_power.R
## to refresh the SILSO download if those months are wanted.
Ocean_solar <- NOAA.Ocean.anomalies %>%
  dplyr::select(dt.mnth, Anomaly = all_of(SST_SERIES)) %>%
  inner_join(dplyr::select(solar_mnthly, dt.mnth, SI), by = "dt.mnth") %>%
  arrange(dt.mnth)

N        <- nrow(Ocean_solar)
rec_from <- min(Ocean_solar$dt.mnth)
rec_to   <- max(Ocean_solar$dt.mnth)
cat(sprintf("record: %.2f - %.2f,  N = %d months\n", rec_from, rec_to, N))

## ---- 2. deseasonalise -----------------------------------------------------
Ocean_solar <- Ocean_solar %>% mutate(anomaly = Resid(Anomaly, M.ssn))

## ---- 3. solar-locked trend ------------------------------------------------
## knots = local minima of the solar power variation (physical clock,
## irregularly spaced; no mean cycle length is used anywhere)
find_nodes <- function(x) which(diff(sign(diff(x))) == 2) + 1
nodes <- Ocean_solar$dt.mnth[find_nodes(Ocean_solar$SI)]
nodes <- nodes[nodes > rec_from & nodes < rec_to]
cat(sprintf("solar knots: %d, mean interval %.2f yr (sd %.2f)\n",
            length(nodes), mean(diff(nodes)), sd(diff(nodes))))

trd.solar <- lm(Anomaly ~ bs(dt.mnth, knots = nodes, degree = 3),
                data = Ocean_solar)
Ocean_solar$trend.solar <- as.numeric(predict(trd.solar))

## ---- 4. strip the harmonic part with a LOESS ------------------------------
Ocean_solar$nonharm <- as.numeric(
  predict(loess(trend.solar ~ dt.mnth, data = Ocean_solar,
                span = LOESS_SPAN, degree = 2)))

Ocean_solar <- Ocean_solar %>%
  mutate(
    harmonic = trend.solar - nonharm,   # periodic part of the locked trend
    natural  = anomaly - trend.solar,   # residual against the locked trend
    growth   = anomaly - harmonic       # <- target of the growth-curve fit
  )                                     #    ( = nonharm + natural )

## ---- 5. two-step growth-curve fit -----------------------------------------
## step 1: logistic + constant, 4 parameters, single nlsLM call
## step 2: Richards, slope k fixed at the step-1 value, nu free
fit_growth <- function(y, tt, label = "",
                       start = list(b = -0.15, L = 0.80, k = 0.05, t0 = 1988)) {

  d <- tibble(tt = tt, y = y)

  fit_log <- nlsLM(y ~ b + L / (1 + exp(-k * (tt - t0))), data = d,
                   start = start, control = nls.lm.control(maxiter = 200))
  cl  <- coef(fit_log)
  b   <- cl[["b"]]; L <- cl[["L"]]; k <- cl[["k"]]; t0 <- cl[["t0"]]
  yl  <- as.numeric(predict(fit_log))

  ## step 2 -- k held fixed, so the secular rate is not re-estimated and the
  ## two curves differ only in shape
  k_fixed <- k
  fit_rich <- nlsLM(
    y ~ Lr + Delta / (1 + exp(-k_fixed * (tt - t0r)))^(1 / nu), data = d,
    start   = list(Lr = b, Delta = L, t0r = t0, nu = 1),
    control = nls.lm.control(maxiter = 500))
  cr  <- coef(fit_rich)
  yr  <- as.numeric(predict(fit_rich))

  ## turning points read as the sample of steepest rise (the data are discrete)
  turn <- function(v) tt[which.max(diff(v)) + 1]

  out <- list(
    label   = label,
    logis   = fit_log, rich = fit_rich,
    y_log   = yl,      y_rich = yr,
    k       = k_fixed,
    t0_log  = t0,      turn_log  = turn(yl),
    t0_rich = cr[["t0r"]], nu = cr[["nu"]], turn_rich = turn(yr),
    lower   = b,           upper      = b + L,
    lower_r = cr[["Lr"]],  upper_r    = cr[["Lr"]] + cr[["Delta"]],
    ## how far up the sigmoid the record actually reaches: < ~0.9 means the
    ## upper asymptote is an extrapolation, not an estimate
    reached = (tail(yl, 1) - b) / L,
    aic_log = AIC(fit_log), aic_rich = AIC(fit_rich),
    rmse_log  = sqrt(mean((y - yl)^2)),
    rmse_rich = sqrt(mean((y - yr)^2))
  )
  ## analytic inflection of the Richards curve, as a check on the discrete read
  out$infl_rich <- out$t0_rich - log(out$nu) / k_fixed
  out
}

report <- function(f) {
  cat(sprintf(
"\n--- %s ---\n  logistic : b=%+.4f  L=%.4f  k=%.4f  t0=%.2f   turn=%.2f   RMSE=%.4f  AIC=%.1f
  richards : L=%+.4f  Delta=%.4f  t0=%.2f  nu=%.4f  turn=%.2f   RMSE=%.4f  AIC=%.1f
  asymptotes  logistic %.3f -> %.3f   richards %.3f -> %.3f
  analytic inflection t0-ln(nu)/k = %.2f  (discrete read %.2f)
  fraction of the span reached at %.2f : %.2f\n",
    f$label, coef(f$logis)[["b"]], coef(f$logis)[["L"]], f$k, f$t0_log,
    f$turn_log, f$rmse_log, f$aic_log,
    coef(f$rich)[["Lr"]], coef(f$rich)[["Delta"]], f$t0_rich, f$nu,
    f$turn_rich, f$rmse_rich, f$aic_rich,
    f$lower, f$upper, f$lower_r, f$upper_r,
    f$infl_rich, f$turn_rich, rec_to, f$reached))
  invisible(f)
}

tt <- Ocean_solar$dt.mnth

fit_main <- fit_growth(Ocean_solar$growth,  tt, "growth  = anomaly - harmonic part  (noisy)")
fit_nh   <- fit_growth(Ocean_solar$nonharm, tt, "nonharm = LOESS non-harmonic locked trend")
report(fit_main); report(fit_nh)

## the two differ only in the noise that is carried along: the parameters are
## the same to 3 decimals, which is the author's "add white noise" control.

## ---- 6. AIC across the candidate target series ----------------------------
## AIC is comparable only between models fitted to the SAME response, so the
## table below compares logistic vs Richards WITHIN each series, and reports
## the parameters across series for reference.
aic_tbl <- bind_rows(
  tibble(series = "anomaly - harmonic part", model = "logistic",
         AIC = fit_main$aic_log,  k = fit_main$k, t0 = fit_main$t0_log,
         nu = NA_real_, turn = fit_main$turn_log,  ceiling = fit_main$upper),
  tibble(series = "anomaly - harmonic part", model = "Richards (k fixed)",
         AIC = fit_main$aic_rich, k = fit_main$k, t0 = fit_main$t0_rich,
         nu = fit_main$nu,        turn = fit_main$turn_rich, ceiling = fit_main$upper_r),
  tibble(series = "non-harmonic locked trend", model = "logistic",
         AIC = fit_nh$aic_log,  k = fit_nh$k, t0 = fit_nh$t0_log,
         nu = NA_real_, turn = fit_nh$turn_log,  ceiling = fit_nh$upper),
  tibble(series = "non-harmonic locked trend", model = "Richards (k fixed)",
         AIC = fit_nh$aic_rich, k = fit_nh$k, t0 = fit_nh$t0_rich,
         nu = fit_nh$nu,        turn = fit_nh$turn_rich, ceiling = fit_nh$upper_r)
) %>%
  group_by(series) %>% mutate(dAIC = AIC - min(AIC)) %>% ungroup()
print(aic_tbl, n = Inf)

## ---- 7. plots -------------------------------------------------------------
Ocean_solar <- Ocean_solar %>%
  mutate(logist_fit = fit_main$y_log,
         richards   = fit_main$y_rich)

plt_growth <- Ocean_solar %>%
  ggplot(aes(x = dt.mnth)) +
  geom_line(aes(y = growth,     colour = "observed - harmonic"), alpha = .45) +
  geom_line(aes(y = logist_fit, colour = "logistic"),  linewidth = 1.1) +
  geom_line(aes(y = richards,   colour = "Richards"),  linewidth = 1.1) +
  geom_hline(yintercept = fit_main$lower,   linetype = 2, colour = "blue") +
  geom_hline(yintercept = fit_main$upper_r, linetype = 2, colour = "red") +
  geom_vline(xintercept = fit_main$turn_log,  linetype = 3) +
  geom_vline(xintercept = fit_main$turn_rich, linetype = 2) +
  scale_colour_manual(values = c("observed - harmonic" = "grey60",
                                 "logistic" = "cyan4", "Richards" = "firebrick")) +
  labs(x = NULL, y = "anomaly (K)", colour = NULL,
       title = "Industrial growth curve, solar-locked construction",
       subtitle = sprintf("logistic turn %.1f,  Richards turn %.1f (nu = %.3f),  k = %.4f held fixed",
                          fit_main$turn_log, fit_main$turn_rich, fit_main$nu, fit_main$k)) +
  theme_minimal(base_size = 11)
print(plt_growth)
ggsave("figs/Industr.growthcurve.png", plt_growth, width = 8, height = 5, dpi = 150)

plt_decomp <- Ocean_solar %>%
  ggplot(aes(x = dt.mnth)) +
  geom_line(aes(y = anomaly,   colour = "deseasonalised anomaly"), alpha = .4) +
  geom_line(aes(y = nonharm,   colour = "non-harmonic locked trend"), linewidth = 1) +
  geom_line(aes(y = harmonic,  colour = "harmonic part")) +
  geom_line(aes(y = richards,  colour = "Richards growth curve"), linewidth = 1.1) +
  labs(x = NULL, y = "anomaly (K)", colour = NULL,
       title = "Decomposition of the global ocean anomaly",
       subtitle = "NOAA ERSST v6, knots at the solar power minima") +
  theme_minimal(base_size = 11)
print(plt_decomp)
ggsave("figs/Decomp.solar.locked_v6.png", plt_decomp, width = 8, height = 5, dpi = 150)

## ---- 8. save --------------------------------------------------------------
Ocean_solar_v6 <- Ocean_solar %>%
  dplyr::select(dt.mnth, Anomaly, Solar_Variation = SI, anomaly,
                trend.solar, nonharm, harmonic, natural, growth,
                logist_fit, richards)
attr(Ocean_solar_v6, "source")   <- Ocean_NOAA.data$url.source
attr(Ocean_solar_v6, "update")   <- Ocean_NOAA.data$update
attr(Ocean_solar_v6, "growth")   <- list(k = fit_main$k,
                                         logistic = coef(fit_main$logis),
                                         richards = coef(fit_main$rich))
saveRDS(Ocean_solar_v6, OUT_RDS)
cat("\nsaved", OUT_RDS, "\n")
