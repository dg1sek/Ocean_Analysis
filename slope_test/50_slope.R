# 50_slope.R -- slope analysis (see PREREG_slope.md). Run from this folder: Rscript 50_slope.R
source("R/alf_tools.R"); set.seed(2026)
rg <- read.csv("data/regional_monthly_ersst.csv"); gl <- read.csv("data/global_series.csv")
ann <- function(df, col) { d <- df[df$year <= 2025, ]; tapply(d[[col]], d$year, mean) }   # complete years only (2026 has 10 months)
S <- list(GLOBAL = ann(gl, "ersst"), MED = ann(rg, "MED"), GOM = ann(rg, "GOM"), WPWP = ann(rg, "WPWP"))
periods <- 60:100; NB <- 500; BL <- 5

fit1 <- function(y, t) {
  sc <- sin_scan(y, t, periods, trend_order = 2); P <- sc$period[which.max(sc$dAIC)]
  x <- (t - 1988) / 100; n <- length(y); phi <- sc$phi[1]
  X <- cbind(1, x, x^2, sin(2*pi*t/P), cos(2*pi*t/P))
  pw <- function(M) rbind(M[1, , drop = FALSE]*sqrt(1-phi^2), M[-1, , drop = FALSE] - phi*M[-n, , drop = FALSE])
  f <- lm.fit(pw(X), pw(matrix(y))); cf <- f$coefficients
  sinus <- as.numeric(X[, 4:5] %*% cf[4:5]); fitted <- as.numeric(X %*% cf)
  win <- function(z, a, b) { i <- t >= a & t < b + 1; 10 * coef(lm(z[i] ~ t[i]))[2] }
  W <- function(z) c(w1 = win(z, 1936, 1965), w2 = win(z, 1966, 1995), w3 = win(z, 1996, 2025))
  list(P = P, c2 = unname(cf[3]), b = unname(cf[2]), A = sqrt(sum(cf[4:5]^2)), dAIC = max(sc$dAIC),
       raw = W(y), nosin = W(y - sinus), fitted = fitted, resid = y - fitted, phi = phi)
}
boot <- function(y, t) {
  f0 <- fit1(y, t); n <- length(y); r <- f0$resid
  out <- t(replicate(NB, { st <- sample.int(n - BL + 1, ceiling(n/BL), TRUE)
    rb <- unlist(lapply(st, function(s) r[s:(s + BL - 1)]))[1:n]
    f <- fit1(f0$fitted + rb, t); c(c2 = f$c2, P = f$P, rawd = unname(f$raw[3] - f$raw[2]), nosd = unname(f$nosin[3] - f$nosin[2]),
      raw3 = unname(f$raw[3]), nos3 = unname(f$nosin[3])) }))
  q <- function(v) quantile(v, c(.025, .975))
  list(f0 = f0, ci = apply(out, 2, q), out = out)
}
rows <- list(); wins <- list()
for (per in list(c(1850, 2025), c(1900, 2025))) for (nm in names(S)) {
  y0 <- S[[nm]]; tt <- as.numeric(names(y0)); k <- tt >= per[1] & tt <= per[2]; y <- as.numeric(y0[k]); t <- tt[k]
  b <- boot(y, t); f <- b$f0; ci <- b$ci
  verdict <- function(lo, hi) if (hi < 0) "slows" else if (lo > 0) "accelerates" else "undetermined"
  rows[[length(rows)+1]] <- data.frame(series = nm, from = per[1], P = f$P, dAIC = round(f$dAIC, 1), sin_amp_K = round(f$A, 2),
    c2 = round(f$c2, 2), c2_lo = round(ci[1, "c2"], 2), c2_hi = round(ci[2, "c2"], 2), Q1 = verdict(ci[1, "c2"], ci[2, "c2"]),
    slope_trend_1996_2025_K_dec = round(10*(f$b/100 + 2*f$c2*(2010-1988)/100^2), 3),
    sinus_slope_amp_K_dec = round(10*2*pi*f$A/f$P, 3),
    raw_w1 = round(f$raw[1], 3), raw_w2 = round(f$raw[2], 3), raw_w3 = round(f$raw[3], 3),
    raw_d = round(f$raw[3]-f$raw[2], 3), raw_d_lo = round(ci[1, "rawd"], 3), raw_d_hi = round(ci[2, "rawd"], 3), Q2raw = verdict(ci[1, "rawd"], ci[2, "rawd"]),
    nos_w1 = round(f$nosin[1], 3), nos_w2 = round(f$nosin[2], 3), nos_w3 = round(f$nosin[3], 3),
    nos_d = round(f$nosin[3]-f$nosin[2], 3), nos_d_lo = round(ci[1, "nosd"], 3), nos_d_hi = round(ci[2, "nosd"], 3), Q2nosin = verdict(ci[1, "nosd"], ci[2, "nosd"]))
  cat(nm, per[1], "done\n")
  if (per[1] == 1850) { x <- (t-1988)/100; wins[[nm]] <- data.frame(series = nm, t = t, y = y, slope_trend = 10*(f$b + 2*f$c2*x)/100,
      slope_sin = NA) }
}
R <- do.call(rbind, rows); write.csv(R, "out/slope_results.csv", row.names = FALSE); print(t(R), quote = FALSE)
