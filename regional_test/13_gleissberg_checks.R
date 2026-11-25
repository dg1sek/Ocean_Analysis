# ---------------------------------------------------------------------------
# 13_gleissberg_checks.R -- is the 80-90 yr line in a region the Atlantic multidecadal oscillation, or the solar Gleissberg cycle?
#  1. AMV proxy from the grid (0-60N, 80W-0E mean minus global-ocean mean), annual 1850-2025: ML sinusoid at 70-100 yr
#  2. Sunspot number 1749-2025: ML sinusoid at 70-100 yr (the Gleissberg modulation of the sunspot record)
#  3. Peak years of the fitted sinusoids at a common period, SST regions vs AMV proxy vs sunspot number
# Output: out/gleissberg_phase_table.csv and the printed table.
# ---------------------------------------------------------------------------
source("R/regional_lib.R"); source("R/alf_tools.R")
d <- read.csv("out/regional_monthly_ersst.csv")
G <- load_ersst_cells()
idx <- which(G$lat >= 0 & G$lat <= 60 & G$lon >= 280 & G$lon <= 360); w <- cos(G$lat[idx] * pi / 180)
na_ <- colSums(G$anom[idx, ] * w, na.rm = TRUE) / colSums((!is.na(G$anom[idx, ])) * w)
wall <- cos(G$lat * pi / 180); gl <- colSums(G$anom * wall, na.rm = TRUE) / colSums((!is.na(G$anom)) * wall)
amv <- annual(G$year, G$month, na_ - gl)
ssn <- annual_ssn()
pk <- function(year, v, per, to = 2) {
  tt <- (year - mean(year)) / (diff(range(year)) / 2)
  X <- cbind(poly(tt, to, raw = TRUE), s = sin(2 * pi * year / per), c = cos(2 * pi * year / per))
  f <- arima(v, order = c(1, 0, 0), xreg = X); f0 <- arima(v, order = c(1, 0, 0), xreg = X[, 1:to])
  cf <- coef(f); amp <- sqrt(cf["s"]^2 + cf["c"]^2); w_ <- 2 * pi / per
  tp <- (atan2(cf["s"], cf["c"]) / w_) %% per                      # time of first maximum within a period (mod P)
  last <- max(year); peak_last <- last - ((last - tp) %% per)       # most recent peak at or before the last year
  data.frame(period = per, amp = round(amp, 3), dAIC = round(f0$aic - f$aic, 1), most_recent_peak_year = round(peak_last))
}
rows <- list()
for (per in c(70, 80, 85, 90, 100)) {
  for (nm in c("AMV_proxy", "SSN_annual", "MED", "GOM", "WPWP")) {
    dat <- if (nm == "AMV_proxy") amv else if (nm == "SSN_annual") data.frame(year = ssn$year, v = ssn$ssn / 100) else annual(d$year, d$month, d[[nm]])
    r <- pk(dat$year, dat$v, per); rows[[length(rows) + 1]] <- cbind(series = nm, years = paste(range(dat$year), collapse = "-"), r)
  }
}
tab <- do.call(rbind, rows); write.csv(tab, "out/gleissberg_phase_table.csv", row.names = FALSE); print(tab, row.names = FALSE)
cat("\nBest single sinusoid at 60-120 yr (AR(1) errors, quadratic trend):\n")
P2 <- exp(seq(log(120), log(60), length.out = 80))
for (nm in c("AMV_proxy", "SSN_annual", "MED", "GOM", "WPWP")) {
  dat <- if (nm == "AMV_proxy") amv else if (nm == "SSN_annual") data.frame(year = ssn$year, v = ssn$ssn / 100) else annual(d$year, d$month, d[[nm]])
  s <- sin_scan(dat$v, dat$year, P2, 2); j <- which.max(s$dAIC)
  cat(sprintf("  %-11s best period %.0f yr, amplitude %.3f, AIC gain %.1f\n", nm, s$period[j], s$amp[j], s$dAIC[j]))
}
cat("\nCorrelation of 31-yr running means (1865-2010):\n")
lo <- function(x) as.numeric(stats::filter(x, rep(1 / 31, 31), sides = 2))
M <- data.frame(year = amv$year, AMV = lo(amv$v))
for (nm in c("MED", "GOM", "WPWP")) { a <- annual(d$year, d$month, d[[nm]]); M[[nm]] <- lo(a$v[match(M$year, a$year)]) }
M <- na.omit(M); print(round(cor(M[, -1]), 2))
