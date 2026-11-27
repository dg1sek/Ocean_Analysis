# How strong is a 60-120 yr variation in solar forcing, and what SST amplitude would it imply? (rough scale argument)
# Inputs: NOAA CDR yearly TSI 1610-2025 (nc) and AR6 solar ERF (csv) from ../research/D1/data/
suppressMessages(library(ncdf4)); dir.create("out", showWarnings = FALSE)
fe <- function(...) { p <- c(...); p[file.exists(p)][1] }; D <- if (file.exists("data/ar6_solar_erf.csv")) "data/" else "../research/D1/data/"
nc <- nc_open(paste0(D, "tsi_v03r00_yearly_s1610_e2025_c20260305.nc")); tsi <- as.numeric(ncvar_get(nc, "TSI")); d <- as.numeric(ncvar_get(nc, "time")); nc_close(nc)
yr <- floor(1610 + d / 365.2425 + 0.5)                       # days since 1610-01-01 -> calendar year (mid-year values)
cat("TSI series:", range(yr), "n =", length(tsi), " mean", round(mean(tsi), 2), " range", round(min(tsi), 2), "-", round(max(tsi), 2), "\n")
S <- read.csv(paste0(D, "ar6_solar_erf.csv")); S <- S[S$year >= 1750 & S$year <= 2019, ]
scan <- function(y, t, periods, trend = 1) { out <- t(sapply(periods, function(P) { X <- cbind(1, poly(t, trend), sin(2 * pi * t / P), cos(2 * pi * t / P)); f <- lm.fit(X, y); cf <- f$coefficients[ncol(X) - 1:0]; c(P, sqrt(sum(cf^2)), 1 - sum(f$residuals^2) / sum((y - mean(y))^2)) })); colnames(out) <- c("period", "amp", "R2"); out }
Ps <- seq(50, 150, 1); res <- list()
for (nm in c("NOAA_TSI_1610_2025", "NOAA_TSI_1749_2025", "NOAA_TSI_1850_2025", "AR6_solar_ERF_1750_2019")) {
  if (grepl("NOAA", nm)) { a <- as.integer(sub(".*_(\\d+)_(\\d+)$", "\\1", nm)); b <- as.integer(sub(".*_(\\d+)_(\\d+)$", "\\2", nm)); k <- yr >= a & yr <= b; y <- tsi[k]; t <- yr[k]; unit <- "W/m2 TSI" } else { y <- S$solar_erf; t <- S$year; unit <- "W/m2 ERF" }
  r <- scan(y, t, Ps); j <- which.max(r[, "R2"]); w <- r[, "period"] >= 60 & r[, "period"] <= 120; j2 <- which(w)[which.max(r[w, "R2"])]
  res[[nm]] <- data.frame(series = nm, unit = unit, best_period_50_150 = r[j, 1], amp_best = r[j, 2], R2_best = r[j, 3], best_period_60_120 = r[j2, 1], amp_60_120 = r[j2, 2], R2_60_120 = r[j2, 3], series_sd = sd(y))
  cat(sprintf("%-26s best 60-120 yr: period %d, amplitude %.3f %s (R2 %.2f); best of 50-150: %d yr, %.3f; series sd %.3f\n", nm, r[j2, 1], r[j2, 2], unit, r[j2, 3], r[j, 1], r[j, 2], sd(y)))
}
R <- do.call(rbind, res); write.csv(R, "out/solar_longperiod_amplitude.csv", row.names = FALSE)
# scale argument: ERF amplitude -> SST amplitude for three response strengths (K per W/m2)
# TSI -> ERF: x 0.7/4 = 0.175 (albedo 0.3, geometry /4); the AR6 solar ERF series is used directly as a cross-check
ampT <- R$amp_60_120[R$series == "NOAA_TSI_1749_2025"]; ampE <- R$amp_60_120[R$series == "AR6_solar_ERF_1750_2019"]
erf_from_tsi <- ampT * 0.175; resp <- c(fast_14months = 0.11, transient_TCR = 1.8 / 3.93, equilibrium_ECS = 3 / 3.93)
cat(sprintf("\nERF amplitude of the 60-120 yr line: %.3f W/m2 (from NOAA TSI x 0.175), %.3f W/m2 (AR6 solar ERF series)\n", erf_from_tsi, ampE))
tab <- data.frame(response = names(resp), K_per_Wm2 = resp, SST_amp_from_TSI = erf_from_tsi * resp, SST_amp_from_AR6 = ampE * resp)
print(tab, digits = 3, row.names = FALSE); write.csv(tab, "out/solar_implied_sst_amplitude.csv", row.names = FALSE)
cat("\nObserved 74-80 yr line: Mediterranean 0.23-0.30 K, Gulf of Mexico 0.23 K (ERSSTv6 annual, 1850-2025).\n")
need <- c(0.23, 0.30) / resp["equilibrium_ECS"]; cat(sprintf("ERF amplitude needed to produce 0.23-0.30 K at equilibrium response (0.76 K per W/m2): %.2f-%.2f W/m2, i.e. %.1f-%.1f times the solar ERF amplitude above (AR6 series)\n", need[1], need[2], need[1] / ampE, need[2] / ampE))

# ---- energy view (Alf's heat-capacity argument, 2026-11-23 11:46): how much solar forcing does it take to move a global upper-ocean layer by the observed amplitude? ----
rho_cp <- 4.1e6            # J m-3 K-1 (sea water)
A_earth <- 5.10e14; A_ocean <- 0.71 * A_earth
cat("\nEnergy view: sinusoidal forcing of amplitude F (W/m2 of ocean area), period P, absorbed by a mixed layer of depth h with NO radiative damping (upper bound for the response):\n")
cat("  T amplitude = F / (rho cp h w),  w = 2 pi / P\n")
e <- expand.grid(h = c(50, 75, 100, 200), P = c(11, 80)); e$w <- 2 * pi / (e$P * 3.156e7); e$F_needed_for_0.3K <- 0.3 * rho_cp * e$h * e$w
e$T_from_solar_ERF_0.017 <- 0.017 / (rho_cp * e$h * e$w); e$ratio_needed_over_solar <- e$F_needed_for_0.3K / 0.017
print(e[, c("h", "P", "F_needed_for_0.3K", "T_from_solar_ERF_0.017", "ratio_needed_over_solar")], digits = 3, row.names = FALSE)
Eheat <- 0.3 * rho_cp * 75 * A_ocean; cat(sprintf("Heat content of 0.3 K in a 75 m layer over the global ocean: %.1e J; solar ERF amplitude 0.017 W/m2 over a half cycle of 40 yr over the same area delivers at most %.1e J (everything absorbed, nothing lost)\n",
    Eheat, 0.017 * A_ocean * 40 * 3.156e7))
write.csv(e, "out/solar_energy_view.csv", row.names = FALSE)
