# A2_run.R : main analysis, parallel over records
base <- "./"   # run with the working directory set to the coral_test folder
wd <- base
source(paste0(wd, "A2_lib.R"))
dat <- readRDS(paste0(wd, "out/records.rds")); R <- dat$R
hale <- build_hale(dat$sun)
saveRDS(hale, paste0(wd, "out/hale.rds"))
cat("Hale cycle minima:", round(hale$minima, 1), "\n")
NSURR <- 2000; NSTAT <- 1000; NDATE <- 200

process <- function(r, hale, seed) {
  set.seed(seed)
  df <- r$df; yrs <- min(df$year):max(df$year)
  x <- approx(df$year, df$x, xout = yrs)$y; ngap <- length(yrs) - nrow(df)
  N <- length(x); xs <- pretreat(x, yrs)
  # --- spectrum
  pg <- pgram(xs); band <- pg$f >= 1 / 25 & pg$f <= 1 / 15
  jb <- which(band)[which.max(pg$p[band])]; pk <- 1 / pg$f[jb]
  fit <- ar1_fit(x - mean(x)); fit2 <- ar1_fit(xs)
  nf <- length(pg$f); SP <- matrix(NA, NSURR, nf); smax <- numeric(NSURR); smax2 <- numeric(NSURR)
  SS <- vector("list", NSTAT)
  for (i in 1:NSURR) {
    s <- ar1_sim(N, fit); ss <- pretreat(s, yrs); q <- pgram(ss)$p; SP[i, ] <- q; smax[i] <- max(q[band])
    if (i <= NSTAT) SS[[i]] <- stat_all(ss, yrs, hale)
    s2 <- ar1_sim(N, fit2); s2 <- (s2 - mean(s2)) / sd(s2); smax2[i] <- max(pgram(s2)$p[band])
  }
  q95 <- apply(SP, 2, quantile, 0.95); qmean <- colMeans(SP)
  spec <- list(f = pg$f, p = pg$p, q95 = q95, qmean = qmean, peakP = pk, peakpow = pg$p[jb],
               pct_point = mean(SP[, jb] < pg$p[jb]), p_band = (1 + sum(smax >= max(pg$p[band]))) / (1 + NSURR),
               p_band_filtfit = (1 + sum(smax2 >= max(pg$p[band]))) / (1 + NSURR), phi = fit$phi, phi_filt = fit2$phi)
  # --- stationarity / phase
  ob <- stat_all(xs, yrs, hale)
  pv <- function(o, s) if (is.na(o)) NA else (1 + sum(s >= o, na.rm = TRUE)) / (1 + sum(!is.na(s)))
  sR <- function(nm) sapply(SS, function(z) z[[nm]])
  pvals <- c(RN = pv(ob$RN, sR("RN")), RS = pv(ob$RS, sR("RS")), RH = pv(ob$RH, sR("RH")), maxr = pv(ob$maxr, sR("maxr")))
  thr <- c(RN = quantile(sR("RN"), .95, na.rm = TRUE), RS = quantile(sR("RS"), .95, na.rm = TRUE), RH = quantile(sR("RH"), .95, na.rm = TRUE))
  # surrogate window-period spread (sd of best period) for comparison
  spS <- sapply(SS, function(z) sd(z$bestP)); spO <- sd(ob$bestP)
  pstab <- (1 + sum(spS <= spO)) / (1 + NSTAT)            # small sd = stable
  # --- dating error
  dstat <- NULL
  if (NDATE > 0) {
    dstat <- lapply(1:NDATE, function(i) {
      w <- cumsum(c(0, rnorm(N - 1, 0, 0.1)))[N:1][N:1]   # placeholder, overwritten below
      w <- rev(cumsum(c(0, rnorm(N - 1, 0, 0.1))))        # error 0 at youngest (last) year, grows toward past
      tn <- yrs + w; xx <- approx(tn, x, xout = yrs)$y
      xx[is.na(xx)] <- x[is.na(xx)]
      xs2 <- pretreat(xx, yrs); z <- stat_all(xs2, yrs, hale)
      c(RN = z$RN, RS = z$RS, RH = z$RH, phN = z$fullN[["ph"]], phS = z$fullS[["ph"]], phH = z$fullH[["ph"]],
        ampN = z$fullN[["amp"]], ampS = z$fullS[["amp"]], bpsd = z$bpsd, maxr = z$maxr)
    })
    dstat <- do.call(rbind, dstat)
  }
  list(name = r$name, N = N, ngap = ngap, yrs = range(yrs), x = x, xs = xs, spec = spec, obs = ob, pvals = pvals, thr = thr,
       pstab = pstab, spO = spO, spSmed = median(spS), dstat = dstat,
       surr_full = t(sapply(SS, function(z) c(N = z$fullN[["ph"]], S = z$fullS[["ph"]], H = z$fullH[["ph"]]))))
}

library(parallel)
cl <- makeCluster(min(7, detectCores() - 1))
wd0 <- normalizePath(".")
clusterExport(cl, c("hale", "NSURR", "NSTAT", "NDATE", "process", "wd0"))
clusterEvalQ(cl, { source(file.path(wd0, "A2_lib.R")); NULL })
seeds <- 2026 + seq_along(R)
t0 <- Sys.time()
res <- clusterMap(cl, function(r, s) process(r, hale, s), R, seeds, SIMPLIFY = FALSE, .scheduling = "dynamic")
stopCluster(cl)
cat("elapsed", format(Sys.time() - t0), "\n")
names(res) <- names(R)
saveRDS(res, paste0(wd, "out/results.rds"))
