## ===========================================================================
##  The natural (total harmonic) component outside the observation window
##  ---------------------------------------------------------------------
##  natural(t) = observed deseasonalised anomaly - fitted Richards curve
##
##  The model of it must carry the WHOLE harmonic component, not a fraction:
##  a fixed-period synthesis recovers only about a third of its variance,
##  because the Schwabe and interannual bands are amplitude- and
##  phase-modulated and a single frequency cannot represent a modulated
##  carrier.  The SSA reconstruction does carry it, and its linear recurrence
##  relation provides the extrapolation.
##
##    primary     SSA reconstruction (L = L_SSA, first R_COMP eigentriples)
##                + recurrent forecast forward, and on the reversed series
##                backward
##    comparison  fixed-period synthesis, damping set to zero
##
##  Run scripts/Industrial_growthcurve.R first; it writes data/Ocean_solar_v6.rds.
## ===========================================================================

library(tidyverse)
library(Rssa)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

L_SSA  <- 480   # window: long enough for the multidecadal band, short enough
R_COMP <- 30    # for a stable LRR.  See the coverage/stability scan below.

Ocean_solar <- readRDS("data/Ocean_solar_v6.rds")
t_obs   <- Ocean_solar$dt.mnth
natural <- Ocean_solar$anomaly - Ocean_solar$richards
N       <- length(natural)
REC     <- max(t_obs) - min(t_obs)
rich_p  <- attr(Ocean_solar, "growth")$richards
k_fix   <- attr(Ocean_solar, "growth")$k
richards_at <- function(x)
  rich_p[["Lr"]] + rich_p[["Delta"]] / (1 + exp(-k_fix * (x - rich_p[["t0r"]])))^(1 / rich_p[["nu"]])

cat(sprintf("total harmonic component: sd %.4f K, N = %d\n", sd(natural), N))

## ---- how much of it is forecastable in principle ---------------------------
## everything below about two years is weather and ENSO noise; no deterministic
## model should be expected to carry it
lowpass <- function(x, cut_yr) {
  X  <- fft(x); fr <- c(0:(N %/% 2), -((N - N %/% 2 - 1):1)) / N * 12
  X[abs(fr) > 1 / cut_yr] <- 0
  Re(fft(X, inverse = TRUE)) / N
}
nat_lp <- lowpass(natural, 2)
cat(sprintf("part at periods > 2 yr: %.0f%% of the variance\n",
            100 * var(nat_lp) / var(natural)))

## ---- 1. SSA reconstruction -------------------------------------------------
s_nat <- ssa(natural, L = L_SSA)
rec   <- as.numeric(reconstruct(s_nat, groups = list(H = 1:R_COMP))$H)
cov_total <- var(rec) / var(natural)
cov_low   <- 1 - var(nat_lp - rec) / var(nat_lp)
cat(sprintf("SSA reconstruction (L = %d, %d eigentriples): %.0f%% of the total variance, %.0f%% of the > 2 yr part\n",
            L_SSA, R_COMP, 100 * cov_total, 100 * cov_low))

## the part left over -- the honest error scale of any extrapolation
noise_sd <- sd(natural - rec)

## ---- 2. recurrent forecast and backcast ------------------------------------
h_fwd <- round((2100 - max(t_obs)) * 12)
h_bwd <- round((min(t_obs) - 1600) * 12)

fwd <- as.numeric(unlist(rforecast(s_nat, groups = list(1:R_COMP),
                                   len = h_fwd, only.new = TRUE)))
s_rev <- ssa(rev(natural), L = L_SSA)
bwd   <- rev(as.numeric(unlist(rforecast(s_rev, groups = list(1:R_COMP),
                                         len = h_bwd, only.new = TRUE))))
t_fwd <- max(t_obs) + seq_len(h_fwd) / 12
t_bwd <- min(t_obs) - rev(seq_len(h_bwd)) / 12
cat(sprintf("recurrent extrapolation: |forecast| max %.3f K, |backcast| max %.3f K\n",
            max(abs(fwd)), max(abs(bwd))))
## A recurrent forecast can diverge if the LRR has roots outside the unit
## circle.  It does not here, but check rather than assume:
stopifnot(max(abs(fwd)) < 5 * sd(natural), max(abs(bwd)) < 5 * sd(natural))

## ---- 3. fixed-period synthesis, for comparison only ------------------------
pairs_nat <- split(1:24, rep(1:12, each = 2)); names(pairs_nat) <- paste0("P", 1:12)
## parestimate returns one "fdimpars.1d" object for a single group and a list
## of them otherwise; normalise so the extraction works either way
as_par_list   <- function(pe) if (inherits(pe, "fdimpars.1d")) list(pe) else pe
esprit_period <- function(p) 1 / abs(unlist(p$frequencies)[1]) / 12
pe <- as_par_list(parestimate(ssa(natural, L = floor(N / 2)),
                              groups = pairs_nat, method = "esprit"))
periods <- imap_dfr(pe, ~ tibble(period = esprit_period(.x))) %>%
  dplyr::filter(period >= 2.5, period <= REC) %>%
  pull(period) %>% round(2) %>% unique() %>% sort()

harm_design <- function(x, P) do.call(cbind, lapply(P, function(p)
  cbind(cos(2 * pi * x / p), sin(2 * pi * x / p))))
m_fix   <- lm(natural ~ harm_design(t_obs, periods))
fix_in  <- as.numeric(predict(m_fix))
pred_fix <- function(x) as.numeric(cbind(1, harm_design(x, periods)) %*% coef(m_fix))
cat(sprintf("fixed-period synthesis (%d periods): %.0f%% of the total variance, %.0f%% of the > 2 yr part\n",
            length(periods), 100 * var(fix_in) / var(natural),
            100 * (1 - var(nat_lp - fix_in) / var(nat_lp))))

## ---- 4. results ------------------------------------------------------------
LIA <- c(1645, 1725)
in_lia <- t_bwd >= LIA[1] & t_bwd <= LIA[2]
cat(sprintf("\nbackcast over %d-%d: mean %+.3f, min %+.3f, max %+.3f K\n",
            LIA[1], LIA[2], mean(bwd[in_lia]), min(bwd[in_lia]), max(bwd[in_lia])))
for (yy in c(2050, 2100)) {
  i <- which.min(abs(t_fwd - yy))
  cat(sprintf("  %d: Richards %+.3f  harmonic %+.3f  total %+.3f K\n",
              yy, richards_at(t_fwd[i]), fwd[i], richards_at(t_fwd[i]) + fwd[i]))
}

## ---- 5. figures ------------------------------------------------------------
dfA <- tibble(dt.mnth = t_obs, total = natural, ssa = rec, fixp = fix_in)
pA <- ggplot(dfA, aes(x = dt.mnth)) +
  geom_line(aes(y = total, colour = "total (observed - Richards)"), linewidth = .25) +
  geom_line(aes(y = ssa,   colour = "SSA reconstruction"), linewidth = .8) +
  geom_line(aes(y = fixp,  colour = "fixed-period synthesis"), linewidth = .6, linetype = 2) +
  geom_hline(yintercept = 0, linewidth = .3) +
  scale_colour_manual(values = c("total (observed - Richards)" = "grey72",
                                 "SSA reconstruction"      = "firebrick",
                                 "fixed-period synthesis"  = "steelblue4")) +
  labs(x = NULL, y = "anomaly (K)", colour = NULL,
       title = "A  Total harmonic component and two models of it",
       subtitle = sprintf("SSA carries %.0f%% of the variance (%.0f%% above 2 yr); fixed periods %.0f%%",
                          100 * cov_total, 100 * cov_low, 100 * var(fix_in) / var(natural))) +
  theme_minimal(base_size = 10)

dfB <- bind_rows(
  tibble(dt.mnth = t_bwd, ssa = bwd,  fixp = pred_fix(t_bwd), part = "backcast"),
  tibble(dt.mnth = t_obs, ssa = rec,  fixp = fix_in,           part = "in record"))
pB <- ggplot(dfB, aes(x = dt.mnth)) +
  annotate("rect", xmin = LIA[1], xmax = LIA[2], ymin = -Inf, ymax = Inf,
           fill = "steelblue", alpha = .12) +
  geom_ribbon(data = dplyr::filter(dfB, part == "backcast"),
              aes(ymin = ssa - 2 * noise_sd, ymax = ssa + 2 * noise_sd),
              fill = "firebrick", alpha = .12) +
  geom_line(aes(y = ssa, colour = part), linewidth = .7) +
  geom_line(data = dplyr::filter(dfB, part == "backcast"),
            aes(y = fixp), colour = "steelblue4", linewidth = .6, linetype = 2) +
  geom_hline(yintercept = 0, linewidth = .3) +
  geom_hline(yintercept = -0.5, linetype = 3, colour = "darkred") +
  annotate("text", x = 1700, y = -0.45, size = 3,
           label = "global LIA cooling in proxy reconstructions ~ -0.5 K") +
  scale_colour_manual(values = c("backcast" = "firebrick", "in record" = "grey45")) +
  coord_cartesian(ylim = c(-0.6, 0.4)) +
  labs(x = NULL, y = "anomaly (K)", colour = NULL,
       title = "B  Backcast of the total harmonic component to 1600",
       subtitle = "shaded: Maunder Minimum / Little Ice Age, 1645-1725; dashed: fixed-period variant") +
  theme_minimal(base_size = 10)

dfC <- bind_rows(
  tibble(dt.mnth = t_obs, harm = rec, part = "in record"),
  tibble(dt.mnth = t_fwd, harm = fwd, part = "forecast")) %>%
  mutate(richards = richards_at(dt.mnth), total = richards + harm) %>%
  dplyr::filter(dt.mnth >= 1980)
pC <- ggplot(dfC, aes(x = dt.mnth)) +
  geom_ribbon(data = dplyr::filter(dfC, part == "forecast"),
              aes(ymin = total - 2 * noise_sd, ymax = total + 2 * noise_sd),
              fill = "firebrick", alpha = .12) +
  geom_line(data = tibble(dt.mnth = t_obs, y = Ocean_solar$anomaly) %>%
              dplyr::filter(dt.mnth >= 1980),
            aes(y = y, colour = "observed anomaly"), linewidth = .25) +
  geom_line(aes(y = richards, colour = "Richards, extended"), linewidth = 1.1) +
  geom_line(aes(y = total, colour = "Richards + harmonic"), linewidth = .6) +
  geom_hline(yintercept = rich_p[["Lr"]] + rich_p[["Delta"]], linetype = 3) +
  geom_vline(xintercept = max(t_obs), linewidth = .3) +
  scale_colour_manual(values = c("observed anomaly" = "grey72",
                                 "Richards, extended" = "darkcyan",
                                 "Richards + harmonic" = "firebrick")) +
  labs(x = NULL, y = "anomaly (K)", colour = NULL,
       title = "C  Forward extension to 2100",
       subtitle = sprintf("ceiling %.2f K is extrapolated, not estimated",
                          rich_p[["Lr"]] + rich_p[["Delta"]])) +
  theme_minimal(base_size = 10)

if (has_patchwork) {
  plt <- pA / pB / pC
  print(plt)
  ggsave("figs/Natural_harmonics_backcast_forecast.png", plt, width = 10, height = 11, dpi = 150)
} else {
  print(pA); print(pB); print(pC)
}

saveRDS(list(in_record = dfA, backcast = tibble(dt.mnth = t_bwd, ssa = bwd),
             forecast  = tibble(dt.mnth = t_fwd, ssa = fwd),
             periods = periods, noise_sd = noise_sd,
             coverage = c(total = cov_total, above2yr = cov_low)),
        "data/Natural_harmonics.rds")
cat("\nsaved data/Natural_harmonics.rds\n")
