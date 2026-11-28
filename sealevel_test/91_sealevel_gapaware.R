# 91_sealevel_gapaware.R -- gap-aware version of the sea-level tests (amendment in PREREG_sealevel.md). Run: Rscript 91_sealevel_gapaware.R
set.seed(2026)
ps <- read.csv("data/psmsl_rlr_annual.csv"); ps <- ps[ps$year <= 2025, ]; Sst <- split(ps, ps$station)
S <- lapply(Sst, function(d) list(t = d$year, y = d$rlr_mm))
med <- do.call(rbind, lapply(c("TRIESTE", "GENOVA", "MARSEILLE"), function(n) { d <- Sst[[n]]; data.frame(year = d$year, anom = d$rlr_mm - mean(d$rlr_mm[d$year >= 1890 & d$year <= 2020])) }))
cm <- tapply(med$anom, med$year, function(z) if (length(z) >= 2) mean(z) else NA); cm <- cm[!is.na(cm)]; S$MED_COMPOSITE <- list(t = as.integer(names(cm)), y = as.numeric(cm))
rg <- read.csv("data/regional_monthly_ersst.csv"); rg <- rg[rg$year <= 2025, ]; mt <- tapply(rg$MED, rg$year, mean); S$MED_SST_ERSST <- list(t = as.integer(names(mt)), y = as.numeric(mt))
per <- 60:100
dm <- function(t, P = NULL) { x <- (t - 1988)/100; cbind(1, x, x^2, if (!is.null(P)) cbind(sin(2*pi*t/P), cos(2*pi*t/P))) }
rss <- function(X, y) sum(lm.fit(X, y)$residuals^2)
scan <- function(y, t) { r0 <- rss(dm(t), y); g <- sapply(per, function(P) length(y)*log(r0/rss(dm(t, P), y)) - 6); j <- which.max(g); c(P = per[j], dAIC = g[j]) }
phi_of <- function(y, t) { f <- lm.fit(dm(t), y); r <- f$residuals; k <- which(diff(t) == 1); sum(r[k]*r[k+1])/sum(r[c(k, k+1)]^2)*2 }   # lag-1 over adjacent observed years
surr <- function(y, t, nsim, FUN) {
  f <- lm.fit(dm(t), y); mu <- as.numeric(dm(t) %*% f$coefficients); r <- f$residuals; ph <- max(min(phi_of(y, t), .95), 0); sdi <- sd(r)*sqrt(1 - ph^2)
  all <- min(t):max(t); keep <- match(t, all)
  replicate(nsim, { e <- as.numeric(stats::filter(rnorm(length(all) + 100, 0, sdi), ph, method = "recursive"))[-(1:100)][keep]; FUN(mu + e) }) }
P1 <- do.call(rbind, lapply(names(S), function(n) { y <- S[[n]]$y; t <- S[[n]]$t; y <- (y - mean(y))/sd(y); o <- scan(y, t)
  nl <- surr(y, t, 300, function(z) scan(z, t)[["dAIC"]])
  data.frame(series = n, n = length(y), from = min(t), to = max(t), best_P = o[["P"]], dAIC = round(o[["dAIC"]], 1), null95 = round(quantile(nl, .95), 1), p = round(mean(nl >= o[["dAIC"]]), 3), phi = round(phi_of(y, t), 2)) }))
cat("== P1b gap-aware: line at 60-100 yr ==\n"); print(P1, row.names = FALSE); write.csv(P1, "out/sealevel_P1b_line_gapaware.csv", row.names = FALSE)
sk <- function(y, t, Y) { tr <- t <= Y; s <- scan(y[tr], t[tr]); P <- s[["P"]]; c0 <- lm.fit(dm(t[tr]), y[tr])$coefficients; c1 <- lm.fit(dm(t[tr], P), y[tr])$coefficients
  f0 <- as.numeric(dm(t[!tr]) %*% c0); f1 <- as.numeric(dm(t[!tr], P) %*% c1); c(1 - sum((y[!tr] - f1)^2)/sum((y[!tr] - f0)^2), P) }
P2 <- list()
for (n in setdiff(names(S), "MED_SST_ERSST")) { y <- S[[n]]$y; t <- S[[n]]$t; y <- (y - mean(y))/sd(y)
  for (Y in c(1900, 1950, 1975)) { if (sum(t <= Y) < 60 || sum(t > Y) < 25) next
    ob <- sk(y, t, Y); tr <- t <= Y; f <- lm.fit(dm(t[tr]), y[tr]); mu <- as.numeric(dm(t) %*% f$coefficients); r <- f$residuals; ph <- max(min(phi_of(y[tr], t[tr]), .95), 0); sdi <- sd(r)*sqrt(1 - ph^2); all <- min(t):max(t); keep <- match(t, all)
    nul <- replicate(150, { e <- as.numeric(stats::filter(rnorm(length(all) + 100, 0, sdi), ph, method = "recursive"))[-(1:100)][keep]; sk(mu + e, t, Y)[1] })
    P2[[length(P2)+1]] <- data.frame(series = n, cut = Y, n_train = sum(t <= Y), n_fore = sum(t > Y), P_train = ob[2], skill = round(ob[1], 3), null95 = round(quantile(nul, .95), 3), p = round(mean(nul >= ob[1]), 3)); print(P2[[length(P2)]], row.names = FALSE) } }
P2 <- do.call(rbind, P2); write.csv(P2, "out/sealevel_P2b_forecast_gapaware.csv", row.names = FALSE)
cat("\nforecast value (skill > null95 at >= 2 usable cut-offs):\n"); print(tapply(P2$skill > P2$null95, P2$series, function(z) length(z) >= 2 && sum(z) >= 2))
Pfix <- 75; w <- 2*pi/Pfix; wrap <- function(d) ((d + Pfix/2) %% Pfix) - Pfix/2
tmax <- function(cf) { ph <- atan2(cf[4], cf[5]); k <- ceiling((1951*w - ph)/(2*pi)); (ph + 2*pi*k)/w }
pf <- function(y, t) { X <- cbind(dm(t), sin(w*t), cos(w*t)); X <- X[, c(1:3, 4, 5)]; X <- cbind(1, (t-1988)/100, ((t-1988)/100)^2, sin(w*t), cos(w*t)); f <- lm.fit(X, y); list(cf = f$coefficients, res = f$residuals, X = X) }
P3 <- do.call(rbind, lapply(names(S), function(n) { k <- S[[n]]$t >= 1850; t <- S[[n]]$t[k]; y <- S[[n]]$y[k]; if (length(t) < 70) return(NULL); f <- pf(y, t); tm <- tmax(f$cf); N <- length(y); BL <- 5
  bs <- replicate(200, { st <- sample.int(N - BL + 1, ceiling(N/BL), TRUE); rb <- unlist(lapply(st, function(s) f$res[s:(s + BL - 1)]))[1:N]; b <- pf(as.numeric(f$X %*% f$cf) + rb, t); tmax(b$cf) })
  d <- wrap(bs - tm); data.frame(series = n, from = min(t), n = N, last_max_year = round(tm, 1), lo = round(tm + quantile(d, .025), 1), hi = round(tm + quantile(d, .975), 1), amp_over_sd = round(sqrt(f$cf[4]^2 + f$cf[5]^2)/sd(y), 2)) }))
cat("\n== P3b gap-aware: year of the last maximum of a 75-yr sinusoid (quadratic removed), all valid years ==\n"); print(P3, row.names = FALSE); write.csv(P3, "out/sealevel_P3b_phase_gapaware.csv", row.names = FALSE)
