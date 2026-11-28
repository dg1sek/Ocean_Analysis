# 90_sealevel.R -- tide gauges as an independent test of the 70-80 yr line (PREREG_sealevel.md). Run: Rscript 90_sealevel.R
source("R/alf_tools.R"); set.seed(2026)
ps <- read.csv("data/psmsl_rlr_annual.csv"); ps <- ps[ps$year <= 2025, ]
cont <- function(yr, v) {                  # fill gaps <= 3 yr, longest continuous stretch
  all <- min(yr):max(yr); x <- rep(NA_real_, length(all)); x[match(yr, all)] <- v
  i <- which(!is.na(x)); nb <- length(all)
  for (k in seq_along(i)[-1]) { g <- i[k] - i[k-1] - 1; if (g >= 1 && g <= 3) x[(i[k-1]+1):(i[k]-1)] <- approx(c(i[k-1], i[k]), c(x[i[k-1]], x[i[k]]), (i[k-1]+1):(i[k]-1))$y }
  ok <- !is.na(x); r <- rle(ok); e <- cumsum(r$lengths); s <- e - r$lengths + 1; j <- which(r$values); j <- j[which.max(r$lengths[j])]
  idx <- s[j]:e[j]; list(t = all[idx], y = x[idx])
}
Sst <- split(ps, ps$station); S <- lapply(Sst, function(d) cont(d$year, d$rlr_mm))
# Mediterranean composite
med <- do.call(rbind, lapply(c("TRIESTE", "GENOVA", "MARSEILLE"), function(n) { d <- Sst[[n]]; d$anom <- d$rlr_mm - mean(d$rlr_mm[d$year >= 1890 & d$year <= 2020]); d[, c("year", "anom")] }))
cm <- tapply(med$anom, med$year, function(z) if (length(z) >= 2) mean(z) else NA); cm <- cm[!is.na(cm)]
S$MED_COMPOSITE <- cont(as.integer(names(cm)), as.numeric(cm))
for (n in names(S)) cat(sprintf("%-14s %d-%d (n = %d)\n", n, min(S[[n]]$t), max(S[[n]]$t), length(S[[n]]$t)))
# ---- P1
per <- 60:100; P1 <- do.call(rbind, lapply(names(S), function(n) { y <- S[[n]]$y; t <- S[[n]]$t; y <- (y - mean(y))/sd(y)   # standardised so that amplitudes are in sd
  s <- sin_scan(y, t, per, 2); j <- which.max(s$dAIC); nl <- null_scan(y, t, per, 2, nsim = 300)
  data.frame(series = n, n = length(y), from = min(t), to = max(t), best_P = s$period[j], amp_sd = round(s$amp[j], 2), dAIC = round(s$dAIC[j], 1), null95 = round(quantile(nl$best, .95), 1), p = round(mean(nl$best >= s$dAIC[j]), 3), phi = round(s$phi[1], 2)) }))
cat("\n== P1: line at 60-100 yr ==\n"); print(P1, row.names = FALSE); write.csv(P1, "out/sealevel_P1_line.csv", row.names = FALSE)
# ---- P2
design <- function(t, P = NULL) { x <- (t - 1988)/100; cbind(1, x, x^2, if (!is.null(P)) cbind(sin(2*pi*t/P), cos(2*pi*t/P))) }
pwfit <- function(X, y, phi) { n <- length(y); pw <- function(M) rbind(M[1, , drop = FALSE]*sqrt(1-phi^2), M[-1, , drop = FALSE] - phi*M[-n, , drop = FALSE]); lm.fit(pw(X), pw(matrix(y)))$coefficients }
skill1 <- function(y, t, Y) { tr <- t <= Y; ytr <- y[tr]; ttr <- t[tr]; tf <- t[!tr]; yf <- y[!tr]
  s <- sin_scan(ytr, ttr, per, 2); P <- s$period[which.max(s$dAIC)]; phi <- s$phi[1]
  c0 <- pwfit(design(ttr), ytr, phi); c1 <- pwfit(design(ttr, P), ytr, phi)
  f0 <- as.numeric(design(tf) %*% c0); f1 <- as.numeric(design(tf, P) %*% c1); c(skill = 1 - sum((yf - f1)^2)/sum((yf - f0)^2), P = P) }
P2 <- list()
for (n in names(S)) { y <- S[[n]]$y; t <- S[[n]]$t; y <- (y - mean(y))/sd(y); N <- length(y)
  for (Y in c(1900, 1950, 1975)) { if (sum(t <= Y) < 70 || sum(t > Y) < 25) next
    ob <- skill1(y, t, Y); tr <- t <= Y; s0 <- sin_scan(y[tr], t[tr], per, 2); phi <- s0$phi[1]
    c0 <- pwfit(design(t[tr]), y[tr], phi); mu <- as.numeric(design(t) %*% c0); r <- y[tr] - as.numeric(design(t[tr]) %*% c0); sdi <- sd(r)*sqrt(1 - phi^2)
    nul <- replicate(150, { e <- as.numeric(stats::filter(rnorm(N + 100, 0, sdi), phi, method = "recursive"))[-(1:100)]; skill1(mu + e, t, Y)[["skill"]] })
    P2[[length(P2)+1]] <- data.frame(series = n, cut = Y, P_train = ob[["P"]], skill = round(ob[["skill"]], 3), null95 = round(quantile(nul, .95), 3), p = round(mean(nul >= ob[["skill"]]), 3)); print(P2[[length(P2)]], row.names = FALSE) } }
P2 <- do.call(rbind, P2); write.csv(P2, "out/sealevel_P2_forecast.csv", row.names = FALSE)
cat("\nforecast value (skill > null95 at >= 2 usable cut-offs):\n"); print(tapply(P2$skill > P2$null95, P2$series, function(z) length(z) >= 2 && sum(z) >= 2))
# ---- P3 phase at fixed P = 75 over 1890-2025
rg <- read.csv("data/regional_monthly_ersst.csv"); rg <- rg[rg$year <= 2025, ]; mt <- tapply(rg$MED, rg$year, mean)
S$MED_SST_ERSST <- list(t = as.integer(names(mt)), y = as.numeric(mt))
Pfix <- 75; w <- 2*pi/Pfix
tmax <- function(cf) { phi <- atan2(cf[4], cf[5]); k <- ceiling((1951*w - phi)/(2*pi)); (phi + 2*pi*k)/w }   # maximum of a sin + b cos in (1950, 2026]
phase_fit <- function(y, t) { X <- cbind(1, (t - 1988)/100, ((t - 1988)/100)^2, sin(w*t), cos(w*t)); cf <- lm.fit(X, y)$coefficients; list(cf = cf, res = y - X %*% cf, X = X) }
wrap <- function(d) ((d + Pfix/2) %% Pfix) - Pfix/2
P3 <- do.call(rbind, lapply(names(S), function(n) { k <- S[[n]]$t >= 1890; t <- S[[n]]$t[k]; y <- S[[n]]$y[k]; if (length(t) < 70) return(NULL)
  f <- phase_fit(y, t); tm <- tmax(f$cf); amp <- sqrt(f$cf[4]^2 + f$cf[5]^2); N <- length(y); BL <- 5
  bs <- replicate(200, { st <- sample.int(N - BL + 1, ceiling(N/BL), TRUE); rb <- unlist(lapply(st, function(s) f$res[s:(s + BL - 1)]))[1:N]; yb <- as.numeric(f$X %*% f$cf) + rb; b <- phase_fit(yb, t); c(tmax(b$cf), sqrt(b$cf[4]^2 + b$cf[5]^2)) })
  d <- wrap(bs[1, ] - tm); data.frame(series = n, from = min(t), n = N, last_max_year = round(tm, 1), lo = round(tm + quantile(d, .025), 1), hi = round(tm + quantile(d, .975), 1), amp_sd_or_unit = round(amp/sd(y), 2)) }))
cat("\n== P3: year of the last maximum of a 75-yr sinusoid (quadratic trend removed), 1890-2025 ==\n"); print(P3, row.names = FALSE); write.csv(P3, "out/sealevel_P3_phase.csv", row.names = FALSE)
# figure: standardised series with quadratic + best sinusoid
png("out/sealevel_fig.png", 2400, 1500, res = 200); par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
for (n in c("BREST", "TRIESTE", "MED_COMPOSITE", "STOCKHOLM")) { y <- S[[n]]$y; t <- S[[n]]$t; y <- (y - mean(y))/sd(y); s <- sin_scan(y, t, per, 2); P <- s$period[which.max(s$dAIC)]
  X <- design(t, P); cf <- pwfit(X, y, s$phi[1]); plot(t, y, pch = 16, cex = .5, col = "grey50", xlab = "Jahr", ylab = "Meeresspiegel (standardisiert)", main = sprintf("%s, Periode %d J.", n, P)); lines(t, X[, 1:3] %*% cf[1:3], col = "#1f77b4", lwd = 2); lines(t, X %*% cf, col = "#d62728", lwd = 2) }
dev.off()
