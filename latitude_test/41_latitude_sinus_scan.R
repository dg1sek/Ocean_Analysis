# Maximum-likelihood sinusoid scan (4-150 yr, quadratic trend, AR(1) errors) per latitude band, annual means 1850-2025, with AR(1) surrogate null. Design: PREREG_latitude.md.
# Run from this folder after 39_make_latbands.R: Rscript 41_latitude_sinus_scan.R
source("R/alf_tools.R")   # sin_scan(), null_scan() (Alf's steps as functions)
d <- read.csv("data/latband_monthly_ersst.csv", check.names = FALSE); bands <- setdiff(names(d), c("year", "month")); Ps <- seq(4, 150, 0.5); set.seed(20261127)
ann <- function(v) { ok <- tapply(!is.na(v), d$year, sum) == 12; a <- tapply(v, d$year, mean); data.frame(year = as.integer(names(a)), v = as.numeric(a))[ok, ] }
res <- list(); png("out/lat_scan_fig.png", 1800, 1000, res = 150); par(mfrow = c(2, 3), mar = c(4, 4.2, 3, 1))
for (b in bands) { a <- ann(d[[b]]); a <- a[a$year <= 2025, ]; s <- sin_scan(a$v, a$year, Ps, 2); j <- which.max(s$dAIC)
  nl <- null_scan(a$v, a$year, Ps, 2, nsim = 300); p <- mean(nl$best >= s$dAIC[j]); q95 <- quantile(nl$best, .95)
  j80 <- which.min(abs(s$period - 80)); w <- s$period >= 60 & s$period <= 100; j6 <- which(w)[which.max(s$dAIC[w])]
  res[[b]] <- data.frame(band = b, n_years = nrow(a), best_period = s$period[j], amp_K = s$amp[j], AIC_gain = s$dAIC[j], null_q95 = q95, p = p, best_60_100 = s$period[j6], amp_60_100 = s$amp[j6], AIC_gain_60_100 = s$dAIC[j6], phi = s$phi[1])
  cat(sprintf("%-7s n=%d  best %.1f yr, %.2f K, AIC gain %.1f (null 95%% %.1f, p = %.3f); best in 60-100 yr: %.1f yr, %.2f K, gain %.1f\n", b, nrow(a), s$period[j], s$amp[j], s$dAIC[j], q95, p, s$period[j6], s$amp[j6], s$dAIC[j6]))
  plot(s$period, s$dAIC, type = "l", log = "x", xlab = "period (yr)", ylab = "AIC gain", main = b, ylim = c(-8, max(20, s$dAIC[j] + 2))); abline(h = q95, lty = 2); abline(v = 80, col = "grey70", lty = 3) }
dev.off(); R <- do.call(rbind, res); write.csv(R, "out/lat_scan_results.csv", row.names = FALSE); cat("Done.\n")
