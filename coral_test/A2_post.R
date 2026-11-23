base <- "./"   # run with the working directory set to the coral_test folder
wd <- base; source(paste0(wd, "A2_lib.R"))
res <- readRDS(paste0(wd, "out/results.rds")); dat <- readRDS(paste0(wd, "out/records.rds")); R <- dat$R; hale <- readRDS(paste0(wd, "out/hale.rds"))
set.seed(2026)
nm <- names(res)
deg <- function(x) round(x * 180 / pi)
# ---- sign of IDPO relative to Fiji_AB (phase reference only)
yI <- res$FijiTonga_IDPO$obs$y; yF <- res$Fiji_AB_d18O$obs$y
ya <- res$FijiTonga_IDPO$yrs; yb <- res$Fiji_AB_d18O$yrs
cm <- intersect(ya[1]:ya[2], yb[1]:yb[2])
cidpo <- cor(yI[match(cm, ya[1]:ya[2])], yF[match(cm, yb[1]:yb[2])])
idpo_flip <- cidpo < 0; cat("IDPO vs Fiji_AB corr of band-passed:", round(cidpo,2), "
")
# ---- table per record
rows <- lapply(nm, function(n) {
  z <- res[[n]]; s <- z$spec; o <- z$obs; N <- z$N
  rl <- r <- o$r; jr <- if (all(is.na(r))) NA else which.max(abs(r)) - 1
  data.frame(record = n, site = R[[n]]$site, indep = R[[n]]$indep, start = z$yrs[1], end = z$yrs[2], N = N,
    peak_P = round(s$peakP, 1), pct_AR1_point = round(100 * s$pct_point, 1), p_band = round(s$p_band, 3),
    p_band_filtfit = round(s$p_band_filtfit, 3), phi = round(s$phi, 2),
    res_186_20 = N >= 270, res_20_22 = N >= 220, res_186_22 = N >= 120,
    winP_median = round(median(o$bestP), 1), winP_min = min(o$bestP), winP_max = max(o$bestP), winP_sd = round(sd(o$bestP), 2),
    winP_sd_surr_median = round(z$spSmed, 2), p_stable = round(z$pstab, 3),
    nwin = length(o$bestP), amp_win_med = round(median(o$bestA), 3),
    RN = round(o$RN, 2), p_RN = round(z$pvals["RN"], 3), RS = round(o$RS, 2), p_RS = round(z$pvals["RS"], 3),
    RH = round(o$RH, 2), p_RH = round(z$pvals["RH"], 3), maxr = round(o$maxr, 2), lag_maxr = jr, p_maxr = round(z$pvals["maxr"], 3),
    stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows); write.csv(tab, paste0(wd, "out/table_records.csv"), row.names = FALSE)
print(tab[, c("record", "start", "end", "peak_P", "pct_AR1_point", "p_band", "p_band_filtfit", "winP_median", "winP_sd", "winP_sd_surr_median", "p_stable")], row.names = FALSE)
print(tab[, c("record", "RN", "p_RN", "RS", "p_RS", "RH", "p_RH", "maxr", "lag_maxr", "p_maxr")], row.names = FALSE)
# multiple testing note: count of p<0.05 among independent
ind <- tab[tab$indep, ]; cat("\nIndependent-site records: n =", nrow(ind), "\n")
cat("p<0.05 (band):", sum(ind$p_band < .05), " (point pct>=95):", sum(ind$pct_AR1_point >= 95), " RN:", sum(ind$p_RN < .05), " RS:", sum(ind$p_RS < .05), " RH:", sum(ind$p_RH < .05), " maxr:", sum(ind$p_maxr < .05, na.rm = TRUE), "\n")
cat("p<0.05 across all", nrow(tab), "records: band", sum(tab$p_band < .05), " RN", sum(tab$p_RN < .05), " RS", sum(tab$p_RS < .05), " RH", sum(tab$p_RH < .05, na.rm = TRUE), "\n")
# ---- cross-record phases
indn <- nm[tab$indep]
fp <- function(n, which) { o <- res[[n]]$obs; ph <- switch(which, N = o$fullN[["ph"]], S = o$fullS[["ph"]], H = o$fullH[["ph"]]); ph }
cross <- data.frame(record = nm, phN = sapply(nm, fp, "N"), phS = sapply(nm, fp, "S"), phH = sapply(nm, fp, "H"), indep = tab$indep)
cross$phN[cross$record == "FijiTonga_IDPO"] <- cross$phN[cross$record == "FijiTonga_IDPO"] + if (idpo_flip) pi else 0
cross$phS[cross$record == "FijiTonga_IDPO"] <- cross$phS[cross$record == "FijiTonga_IDPO"] + if (idpo_flip) pi else 0
cross$phH[cross$record == "FijiTonga_IDPO"] <- cross$phH[cross$record == "FijiTonga_IDPO"] + if (idpo_flip) pi else 0
cross$deg_N <- deg(cross$phN); cross$deg_S <- deg(cross$phS); cross$deg_H <- deg(cross$phH)
write.csv(cross[, c("record", "indep", "deg_N", "deg_S", "deg_H")], paste0(wd, "out/table_fullphases.csv"), row.names = FALSE)
print(cross[, c("record", "indep", "deg_N", "deg_S", "deg_H")], row.names = FALSE)
ry <- function(sel, col) { p <- cross[[col]][sel]; p <- p[!is.na(p)]; c(n = length(p), Rbar = round(circR(p), 2), mean_deg = deg(circmean(p)), p = round(rayleigh_p(p), 3)) }
ci <- cross$indep
cat("\nRayleigh, independent sites (full record):\n"); print(rbind(nodal = ry(ci, "phN"), sine22 = ry(ci, "phS"), hale = ry(ci, "phH")))
pac <- cross$record %in% c("HonTre_SST", "Fiji_AB_d18O", "Amedee_SrCa", "Rarotonga_SrCa", "Ryukyu_SST")
cat("Rayleigh, Pacific independent only:\n"); print(rbind(nodal = ry(pac, "phN"), sine22 = ry(pac, "phS"), hale = ry(pac, "phH")))
cat("Rayleigh, all coral/proxy records (dependent!):\n"); allp <- !(cross$record %in% c("CET", "PAGES2k_GMST", "ERSSTv6_global")); print(rbind(nodal = ry(allp, "phN"), sine22 = ry(allp, "phS"), hale = ry(allp, "phH")))
# common-period phases 1734-1996, independent sites
cy <- 1750:1996; cp <- t(sapply(indn, function(n) { z <- res[[n]]; yy <- z$yrs[1]:z$yrs[2]; y <- z$obs$y[match(cy, yy)]
  k <- match(cy, hale$years); c(N = ls_phase(cy, y, P_NODAL, T0_NODAL)[["ph"]], S = ls_phase(cy, y, P_SIN, T0_SIN)[["ph"]], H = ls_hale(y, hale$h[k], hale$hq[k])[["ph"]]) }))
cat("\nCommon period 1750-1996 phases (deg):\n"); print(round(cp * 180 / pi))
cat("Rayleigh common period:\n"); print(rbind(nodal = c(Rbar = round(circR(cp[, "N"]), 2), p = round(rayleigh_p(cp[, "N"]), 3)), sine22 = c(round(circR(cp[, "S"]), 2), round(rayleigh_p(cp[, "S"]), 3)), hale = c(round(circR(cp[, "H"]), 2), round(rayleigh_p(cp[, "H"]), 3))))
# surrogate null for cross-record Rbar (independent AR1 surrogates per site, from stored surrogate full phases)
ns <- nrow(res[[1]]$surr_full); cat("surrogate-based cross-site Rbar null (independent sites, 1000 draws):\n")
for (cl in c("N", "S", "H")) { nullR <- replicate(5000, circR(sapply(indn, function(n) { m <- res[[n]]$surr_full[, cl]; m <- m[!is.na(m)]; sample(m, 1) })))
  obsR <- circR(cross[[paste0("ph", cl)]][ci]); cat(cl, "obs Rbar", round(obsR, 2), " p(surrogate) =", round((1 + sum(nullR >= obsR)) / 5001, 3), "\n") }
# ---- dating error
dr <- lapply(nm, function(n) { z <- res[[n]]; d <- z$dstat; o <- z$obs
  data.frame(record = n, RN_obs = round(o$RN, 2), RN_dated_med = round(median(d[, "RN"]), 2), RN_thr95 = round(z$thr["RN.95%"], 2), frac_RN_above_thr = round(mean(d[, "RN"] > z$thr["RN.95%"]), 2),
    RS_obs = round(o$RS, 2), RS_dated_med = round(median(d[, "RS"]), 2), RS_thr95 = round(z$thr["RS.95%"], 2), frac_RS_above = round(mean(d[, "RS"] > z$thr["RS.95%"]), 2),
    RH_obs = round(o$RH, 2), RH_dated_med = round(median(d[, "RH"], na.rm = TRUE), 2), RH_thr95 = round(z$thr["RH.95%"], 2), frac_RH_above = round(mean(d[, "RH"] > z$thr["RH.95%"], na.rm = TRUE), 2),
    fullphaseN_circR = round(circR(d[, "phN"]), 2), fullphaseS_circR = round(circR(d[, "phS"]), 2),
    amp_ratio_N = round(median(d[, "ampN"]) / o$fullN[["amp"]], 2), bp_sd_ratio = round(median(d[, "bpsd"]) / o$bpsd, 2), stringsAsFactors = FALSE) })
dtab <- do.call(rbind, dr); write.csv(dtab, paste0(wd, "out/table_dating.csv"), row.names = FALSE); cat("\nDating-error table\n"); print(dtab, row.names = FALSE)
# cross-site Rayleigh under dating error: per realisation, phases of independent sites
for (cl in c("phN", "phS")) { pp <- sapply(1:200, function(i) rayleigh_p(sapply(indn, function(n) res[[n]]$dstat[i, cl])))
  Rr <- sapply(1:200, function(i) circR(sapply(indn, function(n) res[[n]]$dstat[i, cl])))
  cat(cl, "cross-site Rayleigh under dating errors: median Rbar", round(median(Rr), 2), " median p", round(median(pp), 3), " share p<0.05:", mean(pp < .05), "\n") }

# ---- figures
cols <- c("#1b6ca8", "#c0392b", "#2e8b57", "#8e44ad", "#d68910", "#555555")
nplot <- nm
png(paste0(wd, "fig_spectra.png"), 1800, 1500, res = 130); par(mfrow = c(5, 4), mar = c(3.5, 3.5, 2, 0.8), mgp = c(2.2, 0.7, 0))
for (n in nplot) { s <- res[[n]]$spec; per <- 1 / s$f; ok <- per >= 5 & per <= 60
  plot(per[ok], s$p[ok], type = "l", log = "x", xlab = "period (yr)", ylab = "power", main = sprintf("%s (N=%d)", n, res[[n]]$N), cex.main = 0.85, lwd = 1.5, col = cols[1], ylim = c(0, max(s$p[ok], s$q95[ok]) * 1.05))
  lines(per[ok], s$q95[ok], col = cols[2], lty = 2); lines(per[ok], s$qmean[ok], col = "grey50", lty = 3)
  abline(v = c(18.613, 22), col = c("#2e8b57", "#d68910"), lty = 3); rect(15, -1e9, 25, 1e9, col = rgb(0, 0, 0, 0.05), border = NA)
  mtext(sprintf("peak %.1f yr, pointwise pct %.0f, band p=%.2f", s$peakP, 100 * s$pct_point, s$p_band), 3, -1.2, cex = 0.55) }
plot.new(); legend("center", c("periodogram (10% taper)", "AR(1) 95% (2000 surrogates)", "AR(1) mean", "18.6 yr (green)", "22 yr (orange)", "shaded: 15-25 yr band"), lty = c(1, 2, 3, 3, 3, NA), col = c(cols[1], cols[2], "grey50", "#2e8b57", "#d68910", NA), bty = "n", cex = 0.9)
dev.off()
sel <- c("HonTre_SST", "Fiji_AB_d18O", "Amedee_SrCa", "Rarotonga_SrCa", "DryTortugas_SST", "Ryukyu_SST", "CET", "PAGES2k_GMST", "ERSSTv6_global")
png(paste0(wd, "fig_sliding_period_amp.png"), 1800, 1300, res = 130); par(mfrow = c(3, 3), mar = c(3.5, 3.5, 2, 3.5), mgp = c(2.2, 0.7, 0))
for (n in sel) { o <- res[[n]]$obs; plot(o$cen, o$bestP, type = "b", pch = 19, ylim = c(14.5, 25.5), xlab = "window centre (100-yr windows, step 10)", ylab = "best period (yr)", main = n, col = cols[1], cex.main = 0.9)
  abline(h = c(18.613, 20, 22), lty = 3, col = c("#2e8b57", "grey50", "#d68910")); par(new = TRUE); plot(o$cen, o$bestA, type = "l", col = cols[2], axes = FALSE, xlab = "", ylab = "", ylim = c(0, max(o$bestA) * 1.3)); axis(4, col.axis = cols[2]); mtext("amplitude (SD units)", 4, 2.2, cex = 0.6, col = cols[2]) }
dev.off()
png(paste0(wd, "fig_sliding_phase.png"), 1800, 1300, res = 130); par(mfrow = c(3, 3), mar = c(3.5, 3.5, 2, 0.8), mgp = c(2.2, 0.7, 0))
for (n in sel) { o <- res[[n]]$obs; plot(o$cen, deg(o$phN), type = "b", pch = 19, ylim = c(-180, 180), xlab = "window centre", ylab = "phase (deg)", col = "#2e8b57", main = sprintf("%s  R: nodal %.2f / 22yr %.2f / Hale %.2f", n, o$RN, o$RS, ifelse(is.na(o$RH), NA, o$RH)), cex.main = 0.75)
  lines(o$cen, deg(o$phS), type = "b", pch = 17, col = "#d68910"); lines(o$cen, deg(o$phH), type = "b", pch = 15, col = "#8e44ad"); abline(h = 0, col = "grey80") }
plot.new(); legend("center", c("nodal 18.613 yr (ref. 2006.0)", "22.0-yr sinusoid (ref. 2000.0)", "Hale-signed sunspots"), pch = c(19, 17, 15), col = c("#2e8b57", "#d68910", "#8e44ad"), bty = "n")
dev.off()
# wavelet
morlet <- function(x, per, w0 = 6) { n <- length(x); out <- matrix(0, n, length(per)); x <- x - mean(x)
  for (j in seq_along(per)) { s <- per[j] * w0 / (2 * pi) ; k <- ceiling(4 * s); tt <- -k:k; psi <- exp(1i * w0 * tt / s) * exp(-tt^2 / (2 * s^2)) / sqrt(s)
    cv <- convolve(c(rep(0, k), x, rep(0, k)), Conj(rev(psi)), type = "open"); out[, j] <- abs(cv[(2 * k + 1):(2 * k + n)])^2 }
  out }
per <- seq(8, 40, by = 1); png(paste0(wd, "fig_wavelet.png"), 1800, 1300, res = 130); par(mfrow = c(3, 3), mar = c(3.5, 3.5, 2, 0.8), mgp = c(2.2, 0.7, 0))
for (n in sel) { z <- res[[n]]; yy <- z$yrs[1]:z$yrs[2]; W <- morlet(z$xs, per); image(yy, per, W, col = hcl.colors(30, "YlOrRd", rev = TRUE), xlab = "year", ylab = "period (yr)", main = n, cex.main = 0.9); abline(h = c(18.613, 22), lty = 3, col = "grey30") }
dev.off()
cat("done\n")
