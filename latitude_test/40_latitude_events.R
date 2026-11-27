# Volcano event study per LATITUDE BAND of the ocean (same method as volcano_test/20_volcano_events.R). Design: PREREG_latitude.md.
# Run from this folder after 39_make_latbands.R: Rscript 40_latitude_events.R
set.seed(20261127); dir.create("out", showWarnings = FALSE)
fe <- function(...) { p <- c(...); p[file.exists(p)][1] }
LB <- fe("data/latband_monthly_ersst.csv")
NINO <- fe("data/nino34.long.anom.data"); AOD <- fe("data/ar6_volcanic_sAOD_monthly_1850-2019.csv")
stopifnot(!anyNA(c(LB, NINO, AOD)))
t0 <- 1870; t1 <- 2019; idx <- function(y, m) (y - 1870) * 12 + m            # month index, 1 = Jan 1870
nm <- idx(t1, 12); yr <- 1870 + (0:(nm - 1)) %/% 12; tt <- 1870 + (0:(nm - 1)) / 12
bd <- read.csv(LB, check.names = FALSE)
get <- function(d, col) { v <- d[[col]]; k <- match(paste(d$year, d$month), paste(yr, ((0:(nm - 1)) %% 12) + 1)); out <- rep(NA, nm); ok <- !is.na(k); out[k[ok]] <- v[ok]; out }
# the match above is year-month -> position; rebuild properly:
pos <- function(d) (d$year - 1870) * 12 + d$month
mk <- function(d, col) { o <- rep(NA_real_, nm); p <- pos(d); ok <- p >= 1 & p <= nm; o[p[ok]] <- d[[col]][ok]; o }
S <- lapply(setNames(setdiff(names(bd), c('year', 'month')), setdiff(names(bd), c('year', 'month'))), function(b) mk(bd, b))
# Nino3.4 (PSL format: year + 12 values, -99.99 missing)
nl <- readLines(NINO); nl <- nl[grepl("^\\s*[0-9]{4}\\s+", nl)]; nl <- nl[sapply(strsplit(trimws(nl), "\\s+"), length) == 13]
nn <- do.call(rbind, lapply(strsplit(trimws(nl), "\\s+"), as.numeric)); nino <- rep(NA_real_, nm)
for (i in seq_len(nrow(nn))) { y <- nn[i, 1]; if (y >= 1870 && y <= t1) nino[idx(y, 1:12)] <- ifelse(nn[i, -1] < -90, NA, nn[i, -1]) }
# sAOD (AR6 monthly, decimal years at month centres)
a <- read.csv(AOD); a <- a[a$year >= 1870 & a$year < 2020, ]; sa <- rep(NA_real_, nm); k <- floor((a$year - 1870) * 12 + 1e-6) + 1; sa[k] <- a[[2]]
sa <- sa - median(sa, na.rm = TRUE)
EV <- data.frame(name = c("Krakatau", "Santa Maria", "Novarupta", "Agung", "El Chichon", "Pinatubo"), y = c(1883, 1902, 1912, 1963, 1982, 1991), m = c(8, 10, 6, 3, 4, 6))
EV$i <- idx(EV$y, EV$m); cat("sAOD peak after each event (AR6 monthly, anomaly from median):\n")
for (j in 1:nrow(EV)) cat(sprintf("  %-12s %.3f\n", EV$name[j], max(sa[EV$i[j] + 0:36], na.rm = TRUE)))
mask <- rep(FALSE, nm); for (i in EV$i) mask[max(1, i - 12):min(nm, i + 72)] <- TRUE
lag <- 3; nino_l <- c(rep(NA, lag), nino[1:(nm - lag)])
prep <- function(x) {
  d <- data.frame(x = x, tt = tt, mask = mask); fit <- loess(x ~ tt, d[!d$mask & !is.na(d$x), ], span = 0.3, degree = 1, control = loess.control(surface = "direct"))
  x1 <- x - predict(fit, data.frame(tt = tt)); ok <- !mask & !is.na(x1) & !is.na(nino_l)
  b <- coef(lm(x1[ok] ~ nino_l[ok]))[2]; x2 <- x1 - b * nino_l
  moy <- ((0:(nm - 1)) %% 12) + 1; mm <- tapply(x2[!mask & !is.na(x2)], moy[!mask & !is.na(x2)], mean); x2 <- x2 - mm[moy]      # remove residual seasonal cycle (amendment i)
  list(x = x2, enso_coef = unname(b), trend = x - x1) }
P <- lapply(S, prep)
# ---- A. composites --------------------------------------------------------------------------------------------------
ep <- -12:60; bl <- -36:-7
curve <- function(x, i) sapply(ep, function(e) x[i + e]) - mean(x[i + bl], na.rm = TRUE)
run3 <- function(v) as.numeric(stats::filter(v, rep(1/3, 3), sides = 2))
metrics <- function(cv) { s <- run3(cv); w <- which(ep >= 0 & ep <= 48); j <- w[which.min(s[w])]; pk <- s[j]
  rec <- NA; if (!is.na(pk) && pk < 0) { after <- which(ep > ep[j] & !is.na(s)); hit <- after[which(s[after] >= pk / exp(1))[1]]; if (!is.na(hit)) rec <- ep[hit] - ep[j] }
  c(peak = pk, month_of_min = ep[j], recovery_months = rec) }
res <- list(); curves <- list(); placebo <- list()
okdates <- which(sapply(1:nm, function(i) i > 40 && i < nm - 65 && all(abs(i - EV$i) > 72)))
for (nmS in names(P)) {
  x <- P[[nmS]]$x; C <- sapply(1:nrow(EV), function(j) curve(x, EV$i[j])); colnames(C) <- EV$name; cm <- rowMeans(C, na.rm = TRUE); curves[[nmS]] <- cbind(ep = ep, C, mean = cm)
  mm <- metrics(cm); loo <- t(sapply(1:nrow(EV), function(j) metrics(rowMeans(C[, -j, drop = FALSE], na.rm = TRUE))))
  moyd <- ((okdates - 1) %% 12) + 1; evm <- EV$m   # placebo dates share the calendar month of the matching eruption (amendment ii)
  pl <- t(replicate(2000, { d <- sapply(evm, function(mo) sample(okdates[moyd == mo], 1)); metrics(rowMeans(sapply(d, function(i) curve(x, i)), na.rm = TRUE)) })); placebo[[nmS]] <- pl
  pval <- mean(pl[, "peak"] <= mm["peak"], na.rm = TRUE)
  per_event <- apply(C, 2, function(c) metrics(c)["peak"])
  res[[nmS]] <- data.frame(series = nmS, peak_K = mm["peak"], month_of_min = mm["month_of_min"], recovery_months = mm["recovery_months"],
    peak_loo_min = min(loo[, "peak"], na.rm = TRUE), peak_loo_max = max(loo[, "peak"], na.rm = TRUE), min_month_loo_min = min(loo[, "month_of_min"]), min_month_loo_max = max(loo[, "month_of_min"]),
    placebo_peak_median = median(pl[, "peak"], na.rm = TRUE), placebo_p = pval, enso_coef = P[[nmS]]$enso_coef, row.names = NULL)
  cat(sprintf("\n== %s: peak %.3f K at month %d, recovery %s months; leave-one-out peak %.3f to %.3f; placebo p = %.3f; placebo median peak %.3f; per event peak: %s\n", nmS, mm["peak"], mm["month_of_min"],
      ifelse(is.na(mm["recovery_months"]), "n/a", mm["recovery_months"]), min(loo[, "peak"], na.rm = TRUE), max(loo[, "peak"], na.rm = TRUE), pval, median(pl[, "peak"], na.rm = TRUE), paste(sprintf("%s %.2f", EV$name, per_event), collapse = ", ")))
}
A <- do.call(rbind, res); write.csv(A, "out/lat_volcano_composite_metrics.csv", row.names = FALSE)
for (nmS in names(curves)) write.csv(curves[[nmS]], sprintf("out/lat_volcano_curve_%s.csv", nmS), row.names = FALSE)
# ---- B. one-box response model --------------------------------------------------------------------------------------
sa0 <- ifelse(is.na(sa), 0, sa)
conv <- function(tau) { a <- exp(-1 / tau); as.numeric(stats::filter(sa0, a, method = "recursive")) * (1 - a) }   # causal exponential response, unit gain per unit sAOD
prew <- function(M, phi) { M <- as.matrix(M); rbind(M[1, , drop = FALSE] * sqrt(1 - phi^2), M[-1, , drop = FALSE] - phi * M[-nrow(M), , drop = FALSE]) }
fit_one <- function(x, tau, phi) { ok <- !is.na(x); xx <- x[ok]; U <- conv(tau)[ok]; f <- lm(prew(xx, phi) ~ prew(U, phi) - 1); list(rss = sum(resid(f)^2), g = unname(coef(f)), r = xx - U * coef(f), n = length(xx)) }
taus <- c(1:36, seq(38, 180, 2)); resB <- list(); PHI <- list()
for (nmS in names(P)) {
  x <- P[[nmS]]$x; ok <- !is.na(x); phi <- 0
  for (it in 1:3) { rss <- sapply(taus, function(tau) fit_one(x, tau, phi)$rss); j <- which.min(rss); fo <- fit_one(x, taus[j], phi); phi <- sum(fo$r[-1] * fo$r[-length(fo$r)]) / sum(fo$r^2) }   # iterate: phi from residuals
  PHI[[nmS]] <- phi; rss <- sapply(taus, function(tau) fit_one(x, tau, phi)$rss); j <- which.min(rss); fo <- fit_one(x, taus[j], phi); n <- fo$n
  dev <- n * log(rss / min(rss)); inside <- taus[dev < 3.84]; xx <- x[ok]; rss0 <- sum(prew(xx, phi)^2)
  resB[[nmS]] <- data.frame(series = nmS, tau_months = taus[j], tau_lo = min(inside), tau_hi = max(inside), gain_K_per_sAOD = fo$g, gain_K_per_Wm2 = fo$g / -20, phi = phi, dAIC_gain_vs_no_response = n * log(rss0 / rss[j]) - 4)
  cat(sprintf("%s: tau = %d months (95%% %d-%d), gain %.2f K per sAOD unit (about %.3f K per W/m2 with ERF = -20 x sAOD), AR(1) %.2f, AIC gain vs no response %.1f\n", nmS, taus[j], min(inside), max(inside), fo$g, fo$g / -20, phi, n * log(rss0 / rss[j]) - 4))
}
B <- do.call(rbind, resB); write.csv(B, "out/lat_volcano_onebox_fit.csv", row.names = FALSE)
cat("\nTwo-timescale variant (fast + slow exponential, same sAOD input, AR(1) pre-whitened with the one-box phi):\n")
t1s <- c(1, 2, 3, 4, 6, 9, 12, 18, 24); t2s <- c(36, 60, 90, 120, 180, 300, 600); resC <- list()
for (nmS in names(P)) { x <- P[[nmS]]$x; ok <- !is.na(x); xx <- x[ok]; phi <- PHI[[nmS]]; n <- length(xx); best <- c(Inf, NA, NA); yw <- prew(xx, phi)
  for (a1 in t1s) for (a2 in t2s) { U <- prew(cbind(conv(a1)[ok], conv(a2)[ok]), phi); f <- lm(yw ~ U - 1); s2 <- sum(resid(f)^2); if (s2 < best[1]) best <- c(s2, a1, a2) }
  fo <- fit_one(x, taus[which.min(sapply(taus, function(tau) fit_one(x, tau, phi)$rss))], phi); aic1 <- n * log(fo$rss / n) + 2 * 1; aic2 <- n * log(best[1] / n) + 2 * 2
  U <- prew(cbind(conv(best[2])[ok], conv(best[3])[ok]), phi); gg <- coef(lm(yw ~ U - 1))
  resC[[nmS]] <- data.frame(series = nmS, tau_fast = best[2], tau_slow = best[3], gain_fast = gg[1], gain_slow = gg[2], dAIC_two_minus_one = aic2 - aic1)
  cat(sprintf("%s: best fast %d / slow %d months (gains %.2f / %.2f K per sAOD unit), AIC(two) - AIC(one) = %.1f (negative = two-timescale preferred)\n", nmS, best[2], best[3], gg[1], gg[2], aic2 - aic1)) }
write.csv(do.call(rbind, resC), "out/lat_volcano_twobox_fit.csv", row.names = FALSE)
# ---- figure ----------------------------------------------------------------------------------------------------------
png("out/lat_volcano_fig_composites.png", 1900, 1100, res = 150); par(mfrow = c(2, 3), mar = c(4.2, 4.4, 3, 1))
for (nmS in names(curves)) { cv <- curves[[nmS]]; plb <- placebo[[nmS]]
  plot(cv[, 1], cv[, "mean"], type = "n", ylim = c(-0.8, 0.6), xlab = "months after eruption", ylab = "SST anomaly (K)", main = nmS)
  abline(h = 0, col = "grey80"); abline(v = 0, lty = 3, col = "grey60")
  for (j in 1:nrow(EV)) lines(cv[, 1], cv[, 1 + j], col = adjustcolor("grey55", 0.8), lwd = 1)
  lines(cv[, 1], cv[, "mean"], lwd = 3, col = "#D55E00"); lines(cv[, 1], run3(cv[, "mean"]), lwd = 2, col = "#0072B2")
  mtext(sprintf("peak %.2f K at month %d, placebo p = %.3f", A$peak_K[A$series == nmS], A$month_of_min[A$series == nmS], A$placebo_p[A$series == nmS]), 3, 0.1, cex = 0.7) }
plot.new(); legend("center", c("6-event mean", "3-month running mean", "single eruptions"), col = c("#D55E00", "#0072B2", "grey55"), lwd = c(3, 2, 1), bty = "n", cex = 1.2)
dev.off(); cat("\nDone.\n")
