# ---------------------------------------------------------------------------
# 12_regional_analysis.R <REGION_ID> -- Alf's analysis for one region (ERSSTv6 regional mean, 1850-2025)
#  A  maximum-likelihood sinusoid scan + AIC (his method), with three trend choices and two noise nulls
#  B  the 85-yr (Gleissberg) component: amplitude and AIC gain for different trend choices
#  C  sum of the 4 strongest FFT components; SSA periods of the annual series and an eigenvalue test vs AR(1)
#  D  his knot chain (spline at solar minima minus LOESS): real knots vs 150 random knot sets, periods from FFT and SSA
#  E  solar tests: regression on lagged sunspot number (AR(1) errors), 22-yr Hale-signed regressor
#  F  independent proxy records inside the region (corals etc.): same scan, overlap with ERSSTv6
# Output: out/<ID>_*.csv, out/<ID>_*.png, out/<ID>_log.txt
# ---------------------------------------------------------------------------
source("R/regional_lib.R"); source("R/alf_tools.R")
id <- commandArgs(TRUE)[1]; if (is.na(id)) id <- "GOM"
NSIM <- as.integer(Sys.getenv("NSIM", "400")); NKNOT <- as.integer(Sys.getenv("NKNOT", "150"))
sink(sprintf("out/%s_log.txt", id), split = TRUE)
cat("Region", id, "-", REGIONS$name[REGIONS$id == id], "\n")
d <- read.csv("out/regional_monthly_ersst.csv"); y_m <- d[[id]]
a <- annual(d$year, d$month, y_m); n <- nrow(a)
P <- exp(seq(log(150), log(4), length.out = 400))
cat(sprintf("Annual series %d-%d, n = %d. Periods scanned 4-150 yr.\n", min(a$year), max(a$year), n))

## A -- ML sinusoid scan ------------------------------------------------------------------------------------------
cat("\n== A. Maximum-likelihood sinusoid scan (AR(1) errors), AIC gain = n log(RSS0/RSS1) - 6 ==\n")
peaks <- function(s, k = 3) {
  o <- order(-s$dAIC); sel <- integer(0)
  for (j in o) {
    if (length(sel) == 0 || all(abs(log(s$period[j] / s$period[sel])) > 0.25)) sel <- c(sel, j)
    if (length(sel) >= k) break
  }
  s[sel, ]
}
sc <- list(); A_rows <- list()
for (to in 1:3) {
  s <- sin_scan(a$v, a$year, P, to); sc[[to]] <- s
  nl1 <- null_scan(a$v, a$year, P, to, nsim = if (to == 2) NSIM else NSIM / 2)
  nl2 <- if (to == 2) null_scan_pl(a$v, a$year, P, to, nsim = NSIM) else NULL
  pk <- peaks(s, 3)
  for (r in seq_len(nrow(pk))) {
    A_rows[[length(A_rows) + 1]] <- data.frame(region = id, trend_order = to, rank = r, period_yr = round(pk$period[r], 1),
      amp_K = round(pk$amp[r], 3), dAIC = round(pk$dAIC[r], 1),
      p_AR1 = if (r == 1) round((1 + sum(nl1$best >= pk$dAIC[1])) / (length(nl1$best) + 1), 3) else NA,
      p_powerlaw = if (r == 1 && !is.null(nl2)) round((1 + sum(nl2$best >= pk$dAIC[1])) / (length(nl2$best) + 1), 3) else NA)
  }
  if (to == 2) { null_q95_ar1 <- quantile(nl1$best, .95); null_q95_pl <- quantile(nl2$best, .95); beta <- nl2$beta; phi <- nl1$phi
    cat(sprintf("AR(1) phi = %.2f; fitted power-law exponent (2-30 yr) = %.2f; 95%% of the best AIC gain from noise: AR(1) %.1f, power law %.1f\n",
                phi, beta, null_q95_ar1, null_q95_pl)) }
}
A <- do.call(rbind, A_rows); print(A, row.names = FALSE)
write.csv(sc[[2]], sprintf("out/%s_scan.csv", id), row.names = FALSE)

## B -- the 85-yr (Gleissberg) component --------------------------------------------------------------------------
cat("\n== B. Fixed periods: 85 yr (Gleissberg), 40 and 24 yr (Alf's SSA components), 22 yr (Hale), 11 yr ==\n")
fixed_fit <- function(per, to) {
  tt <- (a$year - mean(a$year)) / (diff(range(a$year)) / 2)
  X <- cbind(poly(tt, to, raw = TRUE), s = sin(2 * pi * a$year / per), c = cos(2 * pi * a$year / per))
  f <- tryCatch(arima(a$v, order = c(1, 0, 0), xreg = X), error = function(e) NULL)
  f0 <- tryCatch(arima(a$v, order = c(1, 0, 0), xreg = X[, 1:to, drop = FALSE]), error = function(e) NULL)
  if (is.null(f) || is.null(f0)) return(NULL)
  cf <- coef(f); bs <- cf["s"]; bc <- cf["c"]; v <- f$var.coef
  amp <- sqrt(bs^2 + bc^2); se <- sqrt((bs^2 * v["s", "s"] + bc^2 * v["c", "c"] + 2 * bs * bc * v["s", "c"]) / amp^2)
  sf <- sin_scan(a$v, a$year, per, to)                       # same definition as the scan in A (phi fixed from the trend-only model)
  data.frame(period = per, trend_order = to, amp_K = round(amp, 3), se_K = round(se, 3),
             dAIC_scan_definition = round(sf$dAIC + 2, 1),      # fixed period: 2 free parameters instead of 3
             dAIC_ML_phi_refitted = round(f0$aic - f$aic, 1))
}
B <- do.call(rbind, lapply(c(85, 40, 24, 22, 11), function(per) do.call(rbind, lapply(1:3, function(to) fixed_fit(per, to)))))
print(B, row.names = FALSE)
write.csv(B, sprintf("out/%s_fixed_periods.csv", id), row.names = FALSE)

## C -- FFT top 4 and SSA --------------------------------------------------------------------------------------------
cat("\n== C. FFT (4 strongest components of the quadratic-detrended annual series) and SSA (annual, L = N/3) ==\n")
rq <- resid(lm(a$v ~ poly(a$year, 2)))
ft <- fft_top(rq, 1, 4); cat("FFT top-4 periods (yr):", paste(round(ft$period, 1), collapse = ", "), "; amplitudes (K):", paste(round(ft$amp, 3), collapse = ", "),
  "; share of variance of the sum:", round(1 - var(rq - (ft$recon - mean(rq))) / var(rq), 2), "\n")
L <- floor(n / 3); suppressMessages(library(Rssa))
ss <- ssa(rq, L = L, svd.method = "eigen"); ev <- ss$sigma^2 / sum(ss$sigma^2)
sp <- tryCatch(ssa_periods(rq, L, 12, 1), error = function(e) NULL)
set.seed(2028); evs <- replicate(150, { x <- as.numeric(stats::filter(rnorm(n + 100, 0, sd(rq) * sqrt(1 - phi^2)), phi, method = "recursive"))[-(1:100)]
  ssx <- ssa(x, L = L, svd.method = "eigen"); ssx$sigma[1:12]^2 / sum(ssx$sigma^2) })
q95 <- apply(evs, 1, quantile, .95)
tabss <- data.frame(component = 1:12, share_pct = round(100 * ev[1:12], 1), AR1_q95_pct = round(100 * q95, 1), above = ev[1:12] > q95)
print(tabss, row.names = FALSE); if (!is.null(sp)) { cat("ESPRIT periods of the component pairs (yr):\n"); print(sp, row.names = FALSE) }
write.csv(tabss, sprintf("out/%s_ssa_eigen.csv", id), row.names = FALSE)

## D -- his knot chain on the monthly series --------------------------------------------------------------------------
cat("\n== D. Knot chain on the monthly series: spline at solar minima minus LOESS(0.25) vs", NKNOT, "random knot sets ==\n")
ok <- !is.na(y_m); t_m <- d$t[ok]; y <- y_m[ok]
kc <- knot_chain(y, t_m, SOLAR_MINIMA); h0 <- kc$harm
domp <- function(h) dominant_period(h, dt = 1 / 12, pad = 16384)
set.seed(2026)
rk <- lapply(seq_len(NKNOT), function(i) knot_chain(y, t_m, random_knots(min(t_m), max(t_m)))$harm)
rr <- data.frame(rms = sapply(rk, sd), dom = sapply(rk, domp), cor = sapply(rk, function(h) cor(h, h0)))
cat(sprintf("his knots (SILSO minima): rms %.3f K, dominant period %.1f yr\n", sd(h0), domp(h0)))
cat(sprintf("random knots: rms median %.3f (5-95%%: %.3f-%.3f), share with larger rms %.2f; dominant period median %.1f (5-95%%: %.1f-%.1f); corr with his curve median %.2f\n",
  median(rr$rms), quantile(rr$rms, .05), quantile(rr$rms, .95), mean(rr$rms >= sd(h0)), median(rr$dom), quantile(rr$dom, .05), quantile(rr$dom, .95), median(rr$cor)))
f4 <- fft_top(h0, 1 / 12, 4); cat("FFT top-4 periods of his harmonic part (yr):", paste(round(f4$period, 1), collapse = ", "), "\n")
sp0 <- tryCatch(ssa_periods(h0, 700, 14, 1 / 12), error = function(e) NULL)
sps <- lapply(1:15, function(i) tryCatch(ssa_periods(rk[[i]], 700, 14, 1 / 12), error = function(e) NULL))
top_per <- function(s) { if (is.null(s)) return(NA); s <- s[is.finite(s$period) & s$period > 2 & s$period < 300 & s$modulus > 0.97 & s$modulus < 1.03, ]
  if (nrow(s) == 0) return(NA); round(head(s$period[order(-s$share_pct)], 3), 1) }
cat("SSA (L=700), ESPRIT periods (yr) of the pairs of components 1-14, his knots:\n"); print(sp0, row.names = FALSE)
cat("pairs with modulus 0.97-1.03 (nearly undamped), ordered by variance share: his knots:", paste(top_per(sp0), collapse = ", "), "\n")
cat("same for 15 random knot sets:\n"); for (i in 1:15) cat("  ", paste(top_per(sps[[i]]), collapse = ", "), "\n")
write.csv(rr, sprintf("out/%s_random_knots.csv", id), row.names = FALSE)

## D2 -- what does 'spline at the solar minima minus LOESS(0.25)' do to a pure sinusoid? (gain of his construction) ----
cat("\n== D2. Gain of the construction (spline at solar minima minus LOESS 0.25) on pure sinusoids of unit amplitude ==\n")
gain <- sapply(c(120, 85, 60, 40, 30, 24, 22, 18, 11, 5), function(per) {
  x <- sin(2 * pi * t_m / per); h <- knot_chain(x, t_m, SOLAR_MINIMA)$harm
  fm <- lm(h ~ sin(2 * pi * t_m / per) + cos(2 * pi * t_m / per)); c(period = per, gain = sqrt(sum(coef(fm)[2:3]^2)))
})
print(round(t(gain), 2), row.names = FALSE); write.csv(t(gain), sprintf("out/%s_construction_gain.csv", id), row.names = FALSE)

## C2 -- SSA window length: can SSA see an ~80-yr line? (Alf's 'factor K' question) ------------------------------------------
cat("\n== C2. SSA of the quadratic-detrended annual series with different window lengths L (the window must exceed the period to resolve it) ==\n")
c2 <- lapply(c(30, 58, 70, 88), function(Lw) {
  sx <- ssa(rq, L = Lw, svd.method = "eigen", kind = "1d-ssa"); ev <- sx$sigma^2 / sum(sx$sigma^2)
  pe <- tryCatch(parestimate(sx, groups = list(1:2), method = "esprit"), error = function(e) NULL)
  rc <- reconstruct(sx, groups = list(1:2))[[1]]
  data.frame(L = Lw, share_pairs_1_2_pct = round(100 * sum(ev[1:2]), 1), period_ESPRIT_yr = if (is.null(pe)) NA else round(1 / pe$frequencies[[1]][1], 1),
             rms_K = round(sd(rc), 3), ML_sinusoid_rms_K = round(sc[[2]]$amp[which.max(sc[[2]]$dAIC)] / sqrt(2), 3))
})
c2 <- do.call(rbind, c2); print(c2, row.names = FALSE); write.csv(c2, sprintf("out/%s_ssa_window.csv", id), row.names = FALSE)

## E -- solar tests -----------------------------------------------------------------------------------------------------
cat("\n== E. Regression of the annual anomaly on lagged annual sunspot number (AR(1) errors, quadratic trend) ==\n")
ssn <- annual_ssn(); E <- solar_lag_fit(a$year, a$v, ssn); print(round(E, 4), row.names = FALSE)
cat(sprintf("smallest p = %.3f at lag %d; Bonferroni over 7 lags: %.3f\n", min(E$p, na.rm = TRUE), E$lag[which.min(E$p)], min(1, 7 * min(E$p, na.rm = TRUE))))
ssn$signed <- ssn$ssn * ifelse(solar_cycle(ssn$year + 0.5) %% 2 == 1, 1, -1)
Es <- solar_lag_fit(a$year, a$v, data.frame(year = ssn$year, ssn = ssn$signed)); cat("Hale-signed sunspot number (22-yr):\n"); print(round(Es, 4), row.names = FALSE)
write.csv(E, sprintf("out/%s_solar_lag.csv", id), row.names = FALSE)

## F -- independent proxies ----------------------------------------------------------------------------------------------
proxy <- read.csv("data/proxy_sites.csv", stringsAsFactors = FALSE); pa <- read.csv("data/proxy_annual.csv", stringsAsFactors = FALSE)
recs <- switch(id, GOM = "DryTortugas_SST", WPWP = "Freund_NWP", SCS = "HonTre_SST", NWPAC = "Ryukyu_SST",
               SWPAC = c("Fiji_AB_d18O", "Amedee_SrCa", "Rarotonga_SrCa"), NINO34 = "Wilson_N34_coral", character(0))
Fp <- list()
if (length(recs)) cat("\n== F. Independent proxy records in the region ==\n")
for (rc in recs) {
  p <- pa[pa$record == rc, ]; p <- p[order(p$year), ]; Pp <- exp(seq(log(min(150, diff(range(p$year)) / 2)), log(4), length.out = 400))
  sp_ <- sin_scan(p$value, p$year, Pp, 2); j <- which.max(sp_$dAIC)
  nlp <- null_scan(p$value, p$year, Pp, 2, nsim = NSIM / 2)
  m <- merge(a, p, by = "year"); lo <- function(x) as.numeric(stats::filter(x, rep(1 / 21, 21), sides = 2))
  hi_r <- if (nrow(m) > 30) cor(resid(lm(m$v ~ poly(m$year, 2))), resid(lm(m$value ~ poly(m$year, 2)))) else NA
  per_e <- sc[[2]]$period[which.max(sc[[2]]$dAIC)]                         # ERSST best period (whole record)
  fx <- function(yr, v, per) { tt <- (yr - mean(yr)) / (diff(range(yr)) / 2)
    X <- cbind(poly(tt, 2, raw = TRUE), s = sin(2 * pi * yr / per), c = cos(2 * pi * yr / per))
    f <- tryCatch(arima(v, order = c(1, 0, 0), xreg = X), error = function(e) NULL); if (is.null(f)) return(c(NA, NA, NA))
    f0 <- arima(v, order = c(1, 0, 0), xreg = X[, 1:2]); cf <- coef(f)
    c(amp = sqrt(cf["s"]^2 + cf["c"]^2), phase_deg = (atan2(cf["c"], cf["s"]) * 180 / pi) %% 360, dAIC = f0$aic - f$aic) }
  pz <- p; pz$value <- (pz$value - mean(pz$value)) / sd(pz$value)
  ov <- pz[pz$year >= 1850, ]; eo <- a[a$year %in% ov$year, ]
  fp_full <- fx(pz$year, pz$value, per_e); fp_ov <- fx(ov$year, ov$value, per_e); fe_ov <- fx(eo$year, eo$v, per_e)
  cat(sprintf("  at the ERSSTv6 best period %.1f yr: proxy full record amp %.2f sd, AIC gain %.1f; proxy 1850+ amp %.2f sd phase %.0f deg; ERSSTv6 same years phase %.0f deg (amplitude %.3f K)\n",
      per_e, fp_full[1], fp_full[3], fp_ov[1], fp_ov[2], fe_ov[2], fe_ov[1]))
  Fp[[rc]] <- data.frame(region = id, record = rc, years = paste(range(p$year), collapse = "-"), n = nrow(p),
    best_period = round(sp_$period[j], 1), amp = round(sp_$amp[j], 3), dAIC = round(sp_$dAIC[j], 1),
    p_AR1 = round((1 + sum(nlp$best >= sp_$dAIC[j])) / (length(nlp$best) + 1), 3),
    n_overlap_ERSST = nrow(m), r_highpass_with_ERSST = round(hi_r, 2), ERSST_best_period = round(per_e, 1),
    proxy_amp_sd_at_ERSST_period = round(fp_full[1], 2), proxy_dAIC_at_ERSST_period = round(fp_full[3], 1),
    phase_proxy_1850plus_deg = round(fp_ov[2]), phase_ERSST_same_years_deg = round(fe_ov[2]))
  print(Fp[[rc]], row.names = FALSE)
}
if (length(Fp)) write.csv(do.call(rbind, Fp), sprintf("out/%s_proxies.csv", id), row.names = FALSE)

## figures ------------------------------------------------------------------------------------------------------------------
png(sprintf("out/%s_fig1_series_fits.png", id), 1500, 900, res = 150); par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 1))
plot(a$year, a$v, type = "l", col = "grey40", xlab = "", ylab = "anomaly (K)", main = sprintf("%s: annual ERSSTv6 regional mean, ML sinusoid (best period) and 4 strongest FFT components", id))
to <- 2; tt <- (a$year - mean(a$year)) / (diff(range(a$year)) / 2); j <- which.max(sc[[2]]$dAIC); per <- sc[[2]]$period[j]
fm <- lm(a$v ~ poly(tt, 2, raw = TRUE) + sin(2 * pi * a$year / per) + cos(2 * pi * a$year / per)); lines(a$year, fitted(fm), col = 2, lwd = 2)
lines(a$year, ft$recon - mean(ft$recon) + mean(fitted(lm(a$v ~ poly(tt, 2, raw = TRUE)))), col = 4, lty = 2)
legend("topleft", c("data", sprintf("trend + ML sinusoid %.0f yr", per), "trend + FFT top-4"), col = c("grey40", 2, 4), lty = c(1, 1, 2), bty = "n")
plot(sc[[2]]$period, sc[[2]]$dAIC, log = "x", type = "l", xlab = "period (yr)", ylab = "AIC gain", main = "AIC gain of one fitted sinusoid vs period; dashed = 95% of best gain from noise")
lines(sc[[1]]$period, sc[[1]]$dAIC, col = "grey60"); lines(sc[[3]]$period, sc[[3]]$dAIC, col = "grey80")
abline(h = c(null_q95_ar1, null_q95_pl), lty = c(2, 3)); abline(v = c(85, 22, 11), col = "orange", lty = 3)
legend("topright", c("trend order 2", "order 1", "order 3", "AR(1) null q95", "power-law null q95"), col = c(1, "grey60", "grey80", 1, 1), lty = c(1, 1, 1, 2, 3), bty = "n", cex = .8); dev.off()
sink()
