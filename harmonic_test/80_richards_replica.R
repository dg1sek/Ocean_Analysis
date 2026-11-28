# 80_richards_replica.R -- see PREREG_richards.md. Run: Rscript 80_richards_replica.R
suppressMessages(library(minpack.lm)); source("R/alf_tools.R"); set.seed(2026)
rg <- read.csv("data/regional_monthly_ersst.csv"); rg <- rg[rg$year <= 2025, ]   # complete years only
y <- as.numeric(tapply(rg$MED, rg$year, mean)); t <- as.numeric(names(tapply(rg$MED, rg$year, mean))); n <- length(y)
ctl <- nls.lm.control(maxiter = 500, ftol = 1e-10, ptol = 1e-10)
tryfit <- function(...) tryCatch(nlsLM(...), error = function(e) NULL)
fit_chain <- function(Kfix = NA, plo = 60, phi = 100) {
  free <- is.na(Kfix); d <- data.frame(t = t, y = y); if (!free) d$Kc <- Kfix
  # step 1: logistic (nu = 1), B free; K free or fixed (column Kc)
  m1 <- if (free) tryfit(y ~ A + (K - A)/(1 + exp(-B*(t - t0))), d, start = list(A = min(y), K = max(y) + .5, B = .1, t0 = 2000), lower = c(A = -3, K = .2, B = 0, t0 = 1850), upper = c(A = 1, K = 8, B = 2, t0 = 2076), control = ctl)
        else tryfit(y ~ A + (Kc - A)/(1 + exp(-B*(t - t0))), d, start = list(A = min(y), B = .1, t0 = 2000), lower = c(A = -3, B = 0, t0 = 1850), upper = c(A = 1, B = 2, t0 = 2076), control = ctl)
  if (is.null(m1)) return(NULL); d$Bc <- coef(m1)[["B"]]
  # step 2: Richards with B fixed (column Bc), several nu starts
  st0 <- as.list(coef(m1)); st0$B <- NULL; best2 <- NULL
  for (nu0 in c(.3, .6, 1.5, 3)) { st <- if (free) list(A = st0$A, K = st0$K, nu = nu0, t0 = st0$t0) else list(A = st0$A, nu = nu0, t0 = st0$t0)   # order must match lower/upper
    mm <- if (free) tryfit(y ~ A + (K - A)/(1 + nu*exp(-Bc*(t - t0)))^(1/nu), d, start = st, lower = c(A = -3, K = .2, nu = .05, t0 = 1850), upper = c(A = 1, K = 8, nu = 10, t0 = 2076), control = ctl)
          else tryfit(y ~ A + (Kc - A)/(1 + nu*exp(-Bc*(t - t0)))^(1/nu), d, start = st, lower = c(A = -3, nu = .05, t0 = 1850), upper = c(A = 1, nu = 10, t0 = 2076), control = ctl)
    if (!is.null(mm) && (is.null(best2) || sum(resid(mm)^2) < sum(resid(best2)^2))) best2 <- mm }
  m2 <- best2; if (is.null(m2)) return(NULL)
  # step 3: add harmonic, scan period (columns cw, sw)
  do.call(rbind, lapply(seq(plo, phi, by = .1), function(P) { w <- 2*pi/P; d$cw <- cos(w*t); d$sw <- sin(w*t); s <- as.list(coef(m2)); s$alpha <- .02; s$beta <- .02
    m3 <- if (free) tryfit(y ~ A + (K - A)/(1 + nu*exp(-Bc*(t - t0)))^(1/nu) + alpha*cw + beta*sw, d, start = s, lower = c(A = -3, K = .2, nu = .05, t0 = 1850, alpha = -2, beta = -2), upper = c(A = 1, K = 8, nu = 10, t0 = 2076, alpha = 2, beta = 2), control = ctl)
          else tryfit(y ~ A + (Kc - A)/(1 + nu*exp(-Bc*(t - t0)))^(1/nu) + alpha*cw + beta*sw, d, start = s, lower = c(A = -3, nu = .05, t0 = 1850, alpha = -2, beta = -2), upper = c(A = 1, nu = 10, t0 = 2076, alpha = 2, beta = 2), control = ctl)
    if (is.null(m3)) return(NULL); cf <- coef(m3); data.frame(P = P, rss = sum(resid(m3)^2), amp = sqrt(cf[["alpha"]]^2 + cf[["beta"]]^2), alpha = cf[["alpha"]], beta = cf[["beta"]], A = cf[["A"]], K = if (free) cf[["K"]] else Kfix, nu = cf[["nu"]], t0 = cf[["t0"]], B = d$Bc[1]) }))
}
aic <- function(rss, k) n*log(rss/n) + 2*k
# free-K chain (k = 8: A, K, nu, t0, B, alpha, beta, P)
fr <- fit_chain(); b <- fr[which.min(fr$rss), ]; fr$aic <- aic(fr$rss, 8)
cat("== Free K (bounds 0.2-8) ==\n"); print(round(b, 3)); cat("K on its bound (>= 7.99):", b$K >= 7.99, "\n")
ok <- fr$P[fr$aic - min(fr$aic) < 2]; cat(sprintf("R1 period: best %.1f yr; AIC within 2: %.1f - %.1f yr\n", b$P, min(ok), max(ok)))
write.csv(fr, "out/richards_free_K_period_scan.csv", row.names = FALSE)
# R2 K profile (k = 7)
Ks <- c(2, 3, 4, 6, 8, 12, 20, 50)
prof <- do.call(rbind, lapply(Ks, function(K) { r <- fit_chain(K); if (is.null(r)) return(data.frame(Kfix = K, P = NA, rss = NA, amp = NA, nu = NA, t0 = NA)); bb <- r[which.min(r$rss), ]; data.frame(Kfix = K, P = bb$P, rss = bb$rss, amp = round(bb$amp, 3), nu = round(bb$nu, 3), t0 = round(bb$t0, 1)) }))
prof$aic <- aic(prof$rss, 7); prof$dAIC <- round(prof$aic - min(prof$aic, na.rm = TRUE), 2); print(prof, row.names = FALSE); write.csv(prof, "out/richards_K_profile.csv", row.names = FALSE)
# R3 quadratic + sinus (ML, AR(1)) on the same data: AIC on the same RSS scale (plain OLS RSS, k = 7: 3 trend + a,b,P + intercept counted in 3)
sc <- sin_scan(y, t, seq(60, 100, by = .1), 2); P2 <- sc$period[which.max(sc$dAIC)]
X <- cbind(1, (t-1988)/100, ((t-1988)/100)^2, sin(2*pi*t/P2), cos(2*pi*t/P2)); rss2 <- sum(lm.fit(X, y)$residuals^2)
cat(sprintf("\nR3 quadratic+sinus (OLS): P = %.1f, RSS = %.3f (k = 6), AIC = %.1f;  best Richards+harmonic (K free): RSS = %.3f (k = 8), AIC = %.1f\n", P2, rss2, aic(rss2, 6), b$rss, min(fr$aic)))
# R4 residuals of the best free-K model
w <- 2*pi/b$P; fit_best <- b$A + (b$K - b$A)/(1 + b$nu*exp(-b$B*(t - b$t0)))^(1/b$nu) + b$alpha*cos(w*t) + b$beta*sin(w*t)
r <- y - fit_best; m <- NULL
phi <- sum(r[-1]*r[-n])/sum(r^2)
pk <- function(z) { z <- z - mean(z); f <- (0:(n-1))/n; I <- Mod(fft(z))^2/n; pe <- 1/f[2:(n %/% 2)]; ii <- pe >= 10 & pe <= 40; c(max(I[2:(n %/% 2)][ii]), pe[ii][which.max(I[2:(n %/% 2)][ii])]) }
ob <- pk(r); sdi <- sd(r)*sqrt(1 - phi^2)
nul <- replicate(300, pk(as.numeric(stats::filter(rnorm(n + 100, 0, sdi), phi, method = "recursive"))[-(1:100)])[1])
cat(sprintf("R4 residuals: AR(1) = %.2f; largest periodogram peak 10-40 yr at %.1f yr, p = %.3f against AR(1) surrogates\n", phi, ob[2], mean(nul >= ob[1])))
write.csv(data.frame(t = t, y = y, fitted = fit_best, resid = r), "out/richards_best_fit.csv", row.names = FALSE)
png("out/richards_fig.png", 2200, 900, res = 200); par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
plot(t, y, pch = 16, cex = .5, col = "grey50", xlab = "Jahr", ylab = "Anomalie (K)", main = "ERSST Mittelmeer, Richards + Harmonische"); lines(t, fit_best, col = "#d62728", lwd = 2)
plot(prof$Kfix, prof$dAIC, type = "b", log = "x", xlab = "K fest (K)", ylab = "AIC - Minimum", main = "Profil über K"); abline(h = 2, lty = 3); dev.off()
