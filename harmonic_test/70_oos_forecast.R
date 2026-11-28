# 70_oos_forecast.R -- out-of-sample test of the harmonic (PREREG_oos_ssa.md, part A). Run: Rscript 70_oos_forecast.R
source("R/alf_tools.R"); set.seed(2026)
rg <- read.csv("data/regional_monthly_ersst.csv"); gl <- read.csv("data/global_series.csv")
ann <- function(df, col) { d <- df[df$year <= 2025, ]; tapply(d[[col]], d$year, mean) }   # complete years only (2026 has 10 months)
S <- list(GLOBAL = ann(gl, "ersst"), MED = ann(rg, "MED"), GOM = ann(rg, "GOM")); cuts <- c(1950, 1975, 2000); NS <- 300
design <- function(t, P = NULL) { x <- (t - 1988) / 100; cbind(1, x, x^2, if (!is.null(P)) cbind(sin(2*pi*t/P), cos(2*pi*t/P))) }
pwfit <- function(X, y, phi) { n <- length(y); pw <- function(M) rbind(M[1, , drop = FALSE]*sqrt(1-phi^2), M[-1, , drop = FALSE] - phi*M[-n, , drop = FALSE]); lm.fit(pw(X), pw(matrix(y)))$coefficients }
skill1 <- function(y, t, Y) {
  tr <- t <= Y; ytr <- y[tr]; ttr <- t[tr]; tf <- t[!tr]; yf <- y[!tr]
  s <- sin_scan(ytr, ttr, 60:100, 2); P <- s$period[which.max(s$dAIC)]; phi <- s$phi[1]
  c0 <- pwfit(design(ttr), ytr, phi); c1 <- pwfit(design(ttr, P), ytr, phi)
  f0 <- as.numeric(design(tf) %*% c0); f1 <- as.numeric(design(tf, P) %*% c1)
  c(skill = 1 - sum((yf - f1)^2) / sum((yf - f0)^2), P = P)
}
res <- list()
for (nm in names(S)) { y <- as.numeric(S[[nm]]); t <- as.numeric(names(S[[nm]])); n <- length(y)
  for (Y in cuts) { obs <- skill1(y, t, Y); skill_obs <- as.numeric(obs["skill"]); P <- as.numeric(obs["P"])
    tr <- t <= Y; ytr <- y[tr]; ttr <- t[tr]; s0 <- sin_scan(ytr, ttr, 60:100, 2); phi <- s0$phi[1]
    c0 <- pwfit(design(ttr), ytr, phi); mu <- as.numeric(design(t) %*% c0)               # M0 trend extended over all years
    r <- ytr - as.numeric(design(ttr) %*% c0); sdi <- sd(r) * sqrt(1 - phi^2)
    null <- replicate(NS, { e <- as.numeric(stats::filter(rnorm(n + 100, 0, sdi), phi, method = "recursive"))[-(1:100)]; as.numeric(skill1(mu + e, t, Y)["skill"]) })
    res[[length(res)+1]] <- data.frame(series = nm, cut = Y, P_train = P, skill = round(skill_obs, 3), null95 = round(quantile(null, .95), 3), p = round(mean(null >= skill_obs), 3))
    print(res[[length(res)]], row.names = FALSE) } }
R <- do.call(rbind, res); write.csv(R, "out/oos_skill.csv", row.names = FALSE)
cat("\nharmonic has forecast value (skill > null95 at >= 2 of 3 cut-offs):\n"); print(tapply(R$skill > R$null95, R$series, sum) >= 2)
