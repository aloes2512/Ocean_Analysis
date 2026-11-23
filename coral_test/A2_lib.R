# A2_lib.R : functions
base <- "./"   # run with the working directory set to the coral_test folder
P_NODAL <- 18.613; T0_NODAL <- 2006.0; P_SIN <- 22.0; T0_SIN <- 2000.0
PGRID <- seq(15, 25, by = 0.1)

pretreat <- function(x, years) {            # LOESS ~100 yr span removed, standardised
  N <- length(x); sp <- min(1, 100 / N)
  fit <- loess(x ~ years, span = sp, degree = 2, family = "gaussian")
  r <- x - fitted(fit); (r - mean(r)) / sd(r)
}
bandpass <- function(x, lo = 15, hi = 25, ramp = 0.2) {   # FFT mask, cosine ramps in frequency
  n <- length(x); nf <- 2^ceiling(log2(2 * n)); xp <- c(x - mean(x), rep(0, nf - n))
  f <- (0:(nf - 1)); f <- ifelse(f <= nf / 2, f, f - nf) / nf; af <- abs(f)
  f1 <- 1 / hi; f2 <- 1 / lo; w <- rep(0, nf)
  w[af >= f1 & af <= f2] <- 1
  lw <- af >= f1 * (1 - ramp) & af < f1; w[lw] <- 0.5 * (1 - cos(pi * (af[lw] - f1 * (1 - ramp)) / (f1 * ramp)))
  uw <- af > f2 & af <= f2 * (1 + ramp); w[uw] <- 0.5 * (1 + cos(pi * (af[uw] - f2) / (f2 * ramp)))
  Re(fft(fft(xp) * w, inverse = TRUE) / nf)[1:n]
}
hilbert_q <- function(x) {                  # quadrature component
  n <- length(x); nf <- 2^ceiling(log2(2 * n)); X <- fft(c(x, rep(0, nf - n)))
  h <- rep(0, nf); h[1] <- 1; h[2:(nf / 2)] <- 2; h[nf / 2 + 1] <- 1
  Im(fft(X * h, inverse = TRUE) / nf)[1:n]
}
pgram <- function(x, pad = 4) {
  n <- length(x); nf <- n * pad; w <- spec.taper(rep(1, n), p = 0.1)
  X <- fft(c((x - mean(x)) * w, rep(0, nf - n)))
  k <- 1:(nf / 2); list(f = k / nf, p = (Mod(X[k + 1])^2) / sum(w^2))
}
ar1_fit <- function(x) { a <- ar.yw(x, aic = FALSE, order.max = 1); list(phi = min(0.99, max(-0.5, a$ar[1])), s2 = a$var.pred) }
ar1_sim <- function(n, fit) { e <- rnorm(n + 100, 0, sqrt(fit$s2)); as.numeric(stats::filter(e, fit$phi, method = "recursive"))[101:(n + 100)] }

# LS at fixed period P with reference t0: returns amplitude and phase (peak offset in rad)
ls_phase <- function(t, y, P, t0) {
  th <- 2 * pi * (t - t0) / P; cc <- cos(th); ss <- sin(th)
  M <- matrix(c(sum(cc * cc), sum(cc * ss), sum(cc * ss), sum(ss * ss)), 2)
  co <- solve(M, c(sum(y * cc), sum(y * ss)))
  c(amp = sqrt(sum(co^2)), ph = atan2(co[2], co[1]))
}
# LS against Hale series h and quadrature hq
ls_hale <- function(y, h, hq) {
  M <- matrix(c(sum(h * h), sum(h * hq), sum(h * hq), sum(hq * hq)), 2)
  co <- solve(M, c(sum(y * h), sum(y * hq))); c(amp = sqrt(sum(co^2)), ph = atan2(co[2], co[1]))
}
# best period (grid) in window
best_period <- function(t, y) {
  t <- t - mean(t); y <- y - mean(y)
  th <- outer(t, 2 * pi / PGRID); C <- cos(th); S <- sin(th)
  Scc <- colSums(C * C); Scs <- colSums(C * S); Sss <- colSums(S * S)
  Syc <- colSums(y * C); Sys <- colSums(y * S); det <- Scc * Sss - Scs^2
  a <- (Syc * Sss - Sys * Scs) / det; b <- (Sys * Scc - Syc * Scs) / det
  amp <- sqrt(a^2 + b^2); j <- which.max(amp); c(P = PGRID[j], amp = amp[j])
}
circR <- function(ph) Mod(mean(exp(1i * ph)))
circmean <- function(ph) Arg(mean(exp(1i * ph)))
rayleigh_p <- function(ph) { n <- length(ph); R <- n * circR(ph); min(1, exp(sqrt(1 + 4 * n + 4 * (n^2 - R^2)) - (1 + 2 * n))) }

# Hale signed sunspot series (annual, band-passed + quadrature), from SILSO monthly
build_hale <- function(s) {
  s <- s[s$sn >= 0, ]
  sm <- stats::filter(s$sn, rep(1 / 13, 13), sides = 2)
  tt <- s$t; cand <- which(!is.na(sm))
  # cycle minima: lowest smoothed value in each window of +-5 yr, kept if local minimum there
  mins <- c()
  for (i in cand) { w <- which(abs(tt - tt[i]) <= 5 & !is.na(sm)); if (sm[i] == min(sm[w])) mins <- c(mins, i) }
  mt <- tt[mins]; ms <- sm[mins]
  # merge near-duplicates (<8 yr apart): keep lowest
  keep <- rep(TRUE, length(mt)); for (i in 2:length(mt)) if (mt[i] - mt[i - 1] < 8) { if (ms[i] < ms[i - 1]) keep[i - 1] <- FALSE else keep[i] <- FALSE }
  mt <- mt[keep]
  yrs <- 1749:2025; sn <- tapply(s$sn, floor(s$t), mean)[as.character(yrs)]
  # cycle index by time: cycle k begins at minimum mt[k]; sign alternates; cycle starting at 1755.2 (cycle 1) = +1
  idx <- findInterval(yrs + 0.5, mt)              # 0 = before first min (cycle 0)
  c1 <- which.min(abs(mt - 1755.2)); sign <- ifelse((idx - c1) %% 2 == 0, 1, -1)
  hs <- sn * sign
  hb <- bandpass(hs - mean(hs), 15, 25)
  list(years = yrs, signed = hs, h = hb, hq = hilbert_q(hb), minima = mt)
}

# all statistics of a pretreated series xs on years yrs
stat_all <- function(xs, yrs, hale, win = 100, step = 10) {
  y <- bandpass(xs); N <- length(y); out <- list()
  st <- seq(1, N - win + 1, by = step)
  W <- lapply(st, function(a) a:(a + win - 1))
  phN <- ampN <- phS <- ampS <- phH <- ampH <- bp <- ba <- cen <- numeric(0)
  for (w in W) {
    t <- yrs[w]; yy <- y[w]
    a <- ls_phase(t, yy, P_NODAL, T0_NODAL); phN <- c(phN, a[2]); ampN <- c(ampN, a[1])
    a <- ls_phase(t, yy, P_SIN, T0_SIN); phS <- c(phS, a[2]); ampS <- c(ampS, a[1])
    k <- match(t, hale$years)
    if (!anyNA(k)) { a <- ls_hale(yy, hale$h[k], hale$hq[k]); phH <- c(phH, a[2]); ampH <- c(ampH, a[1]) } else { phH <- c(phH, NA); ampH <- c(ampH, NA) }
    b <- best_period(t, yy); bp <- c(bp, b[1]); ba <- c(ba, b[2]); cen <- c(cen, mean(t))
  }
  rb <- function(p) if (sum(!is.na(p)) >= 3) circR(p[!is.na(p)]) else NA
  # full-record phases
  fN <- ls_phase(yrs, y, P_NODAL, T0_NODAL); fS <- ls_phase(yrs, y, P_SIN, T0_SIN)
  k <- match(yrs, hale$years); ok <- !is.na(k)
  fH <- if (sum(ok) > 60) ls_hale(y[ok], hale$h[k[ok]], hale$hq[k[ok]]) else c(amp = NA, ph = NA)
  # lagged correlation with Hale series (record lags sunspots by L)
  r <- rep(NA, 6)
  if (sum(ok) > 60) for (L in 0:5) { hh <- hale$h; ii <- which(ok); jj <- ii[ii + L <= N]
    kk <- k[jj]; tgt <- jj + L; sel <- !is.na(kk) & tgt <= N
    r[L + 1] <- cor(y[tgt[sel]], hale$h[kk[sel]]) }
  list(cen = cen, phN = phN, phS = phS, phH = phH, ampN = ampN, ampS = ampS, ampH = ampH, bestP = bp, bestA = ba,
       RN = rb(phN), RS = rb(phS), RH = rb(phH), fullN = fN, fullS = fS, fullH = fH, r = r,
       maxr = if (all(is.na(r))) NA else max(abs(r), na.rm = TRUE), bpsd = sd(y), y = y)
}
