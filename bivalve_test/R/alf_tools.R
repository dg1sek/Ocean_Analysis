# ---------------------------------------------------------------------------
# alf_tools.R -- Alf's analysis steps as functions (base R + Rssa).
#   sin_scan()      maximum-likelihood fit of a*sin + b*cos at every trial period (Alf's "10 lines" idea) + AIC
#   null_scan()     the same scan on AR(1) red-noise surrogates: how large is the best AIC gain from noise alone?
#   fft_top()       sum of the k strongest Fourier components (his "FFT index 6 7 8 9" curve)
#   knot_chain()    his "solar-locked trend": cubic spline with knots at the solar minima, minus LOESS span 0.25
#   ssa_periods()   SSA eigen-triples, periods and moduli of oscillatory pairs (Rssa::parestimate, ESPRIT)
#   solar_lag_fit() regression of the annual anomaly on lagged sunspot number with AR(1) errors
# ---------------------------------------------------------------------------
suppressMessages(library(splines))

# Maximum likelihood = least squares after AR(1) pre-whitening (phi from the trend-only model).
# AIC gain = n*log(RSS0/RSS1) - 2*3 : amplitude a, b and the trial period count as three fitted parameters.
sin_scan <- function(y, t, periods, trend_order = 2, phi = NULL) {
  n <- length(y); tt <- (t - mean(t)) / (diff(range(t)) / 2)
  Xt <- cbind(1, poly(tt, trend_order, raw = TRUE)[, , drop = FALSE])
  if (is.null(phi)) { r0 <- resid(lm.fit(Xt, y)); phi <- sum(r0[-1] * r0[-n]) / sum(r0^2) }
  pw <- function(M) rbind(M[1, , drop = FALSE] * sqrt(1 - phi^2), M[-1, , drop = FALSE] - phi * M[-n, , drop = FALSE])
  yw <- pw(matrix(y)); Xw <- pw(Xt)
  rss0 <- sum(lm.fit(Xw, yw)$residuals^2)
  res <- t(vapply(periods, function(P) {
    S <- cbind(sin(2 * pi * t / P), cos(2 * pi * t / P)); f <- lm.fit(cbind(Xw, pw(S)), yw)
    cf <- f$coefficients[ncol(Xw) + 1:2]
    c(amp = sqrt(sum(cf^2)), phase = atan2(cf[2], cf[1]), dAIC = n * log(rss0 / sum(f$residuals^2)) - 6)
  }, numeric(3)))
  data.frame(period = periods, res, phi = phi)
}

null_scan <- function(y, t, periods, trend_order = 2, nsim = 300, seed = 2026) {
  set.seed(seed); n <- length(y); tt <- (t - mean(t)) / (diff(range(t)) / 2)
  Xt <- cbind(1, poly(tt, trend_order, raw = TRUE)[, , drop = FALSE])
  fit <- lm.fit(Xt, y); r <- resid(fit); phi <- sum(r[-1] * r[-n]) / sum(r^2); sdi <- sd(r) * sqrt(1 - phi^2)
  best <- numeric(nsim); bestP <- numeric(nsim); amp_at <- NULL
  for (i in 1:nsim) {
    ys <- fitted(fit) + as.numeric(stats::filter(rnorm(n + 100, 0, sdi), phi, method = "recursive"))[-(1:100)]
    s <- sin_scan(ys, t, periods, trend_order); j <- which.max(s$dAIC); best[i] <- s$dAIC[j]; bestP[i] <- s$period[j]
  }
  list(best = best, bestP = bestP, phi = phi)
}

fft_top <- function(x, dt = 1, k = 4, pad = NULL) {
  n <- length(x); X <- fft(x - mean(x)); amp <- Mod(X)[2:(n %/% 2 + 1)] * 2 / n
  idx <- order(amp, decreasing = TRUE)[1:k] + 1
  Y <- rep(0 + 0i, n); Y[idx] <- X[idx]; Y[n - idx + 2] <- X[n - idx + 2]
  list(recon = Re(fft(Y, inverse = TRUE)) / n + mean(x), period = n * dt / (idx - 1), amp = Mod(X[idx]) * 2 / n)
}

# Alf's construction: cubic B-spline regression with knots at the solar minima ("solar-locked trend"),
# harmonic part = that spline minus LOESS (span 0.25) of the same series.
knot_chain <- function(y, t, knots, span = 0.25) {
  k <- knots[knots > min(t) & knots < max(t)]
  sp <- as.numeric(fitted(lm(y ~ bs(t, knots = k, degree = 3))))
  lo <- as.numeric(fitted(loess(y ~ t, span = span, degree = 2)))
  list(spline = sp, loess = lo, harm = sp - lo)
}

# ESPRIT period and modulus of each pair of adjacent eigen-triples (1:2, 3:4, ...); modulus ~1 = undamped
ssa_periods <- function(x, L, nmax = 20, dt = 1) {
  suppressMessages(library(Rssa))
  s <- ssa(x, L = L, svd.method = "eigen"); ev <- s$sigma^2 / sum(s$sigma^2)
  pairs <- lapply(seq(1, nmax - 1, by = 2), function(i) c(i, i + 1))
  out <- lapply(pairs, function(g) {
    pe <- tryCatch(parestimate(s, groups = list(g), method = "esprit"), error = function(e) NULL)
    if (is.null(pe)) return(NULL)
    data.frame(pair = paste(g, collapse = "+"), share_pct = round(100 * sum(ev[g]), 1),
               period = (1 / pe$frequencies[[1]] * dt)[1], modulus = pe$moduli[[1]][1])
  })
  do.call(rbind, out)
}

# annual sunspot number from the SILSO monthly file (calendar years with 12 months)
annual_ssn <- function(f = SILSO_FILE) {
  s <- read.table(f, sep = ";", header = FALSE); names(s)[1:4] <- c("y", "m", "t", "sn")
  s <- s[s$sn >= 0, ]; a <- tapply(s$sn, s$y, mean); n <- tapply(s$sn, s$y, length); a <- a[n == 12]
  data.frame(year = as.integer(names(a)), ssn = as.numeric(a))
}

solar_lag_fit <- function(year, y, ssn, lags = 0:6, trend_order = 2) {
  out <- lapply(lags, function(L) {
    d <- merge(data.frame(year = year, y = y), data.frame(year = ssn$year + L, ssn = ssn$ssn), by = "year")
    tt <- (d$year - mean(d$year)) / 100
    X <- cbind(poly(tt, trend_order, raw = TRUE), ssn = d$ssn / 100)
    f <- tryCatch(arima(d$y, order = c(1, 0, 0), xreg = X), error = function(e) NULL)
    if (is.null(f)) return(data.frame(lag = L, K_per_100ssn = NA, se = NA, p = NA))
    b <- tail(coef(f), 1); se <- sqrt(tail(diag(f$var.coef), 1))
    data.frame(lag = L, K_per_100ssn = b, se = se, p = 2 * pnorm(-abs(b / se)))
  })
  do.call(rbind, out)
}

# Second null: spectral-synthesis surrogates with a power-law spectrum fitted ONLY to periods of 2-30 yr and
# extrapolated to long periods (so a genuine multidecadal oscillation cannot inflate its own null).
null_scan_pl <- function(y, t, periods, trend_order = 2, nsim = 300, seed = 2027) {
  set.seed(seed); n <- length(y); tt <- (t - mean(t)) / (diff(range(t)) / 2)
  Xt <- cbind(1, poly(tt, trend_order, raw = TRUE)[, , drop = FALSE])
  fit <- lm.fit(Xt, y); r <- resid(fit)
  Fr <- fft(r - mean(r)); k <- 1:((n - 1) %/% 2); f <- k / n; P <- Mod(Fr[k + 1])^2 / n   # E|F_k|^2 = n * spectral density
  sel <- f >= 1 / 30 & f <= 1 / 2.2
  b <- coef(lm(log(P[sel]) ~ log(f[sel]))); beta <- -b[2]; lvl <- exp(b[1])
  best <- numeric(nsim); bestP <- numeric(nsim)
  for (i in 1:nsim) {
    Z <- rep(0 + 0i, n)
    Z[k + 1] <- sqrt(n * lvl * f^(-beta) * rexp(length(k))) * complex(modulus = 1, argument = runif(length(k), 0, 2 * pi))
    Z[n - k + 1] <- Conj(Z[k + 1])
    xs <- Re(fft(Z, inverse = TRUE)) / n
    s <- sin_scan(fitted(fit) + xs, t, periods, trend_order); j <- which.max(s$dAIC); best[i] <- s$dAIC[j]; bestP[i] <- s$period[j]
  }
  list(best = best, bestP = bestP, beta = beta)
}
