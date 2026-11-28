# 72_oos_richards.R -- forecast test with Richards trend (PREREG_oos2_richards.md). Run: Rscript 72_oos_richards.R
suppressMessages(library(minpack.lm)); suppressMessages(library(parallel)); set.seed(2026, kind = "L'Ecuyer-CMRG")
rg <- read.csv("data/regional_monthly_ersst.csv"); gl <- read.csv("data/global_series.csv")
ann <- function(df, col) { d <- df[df$year <= 2025, ]; tapply(d[[col]], d$year, mean) }
S <- list(GLOBAL = ann(gl, "ersst"), MED = ann(rg, "MED"), GOM = ann(rg, "GOM")); cuts <- c(1950, 1975, 2000); NS <- 100; KF <- 8
ctl <- nls.lm.control(maxiter = 300, ftol = 1e-9, ptol = 1e-9)
tf <- function(...) tryCatch(suppressWarnings(nlsLM(...)), error = function(e) NULL)
rich <- function(t, A, nu, t0, B) A + (KF - A)/(1 + nu*exp(-B*(t - t0)))^(1/nu)
fit0 <- function(ttr, ytr) {     # M0 on training: logistic step for B, then Richards with B fixed
  d <- data.frame(t = ttr, y = ytr)
  m1 <- tf(y ~ A + (KF - A)/(1 + exp(-B*(t - t0))), d, start = list(A = min(ytr), B = .05, t0 = max(ttr) + 20), lower = c(A = -3, B = 0, t0 = 1850), upper = c(A = 1, B = 2, t0 = 2076), control = ctl)
  if (is.null(m1)) return(NULL); B <- coef(m1)[["B"]]; d$Bc <- B; best <- NULL
  for (nu0 in c(.5, 1.5)) { m <- tf(y ~ A + (KF - A)/(1 + nu*exp(-Bc*(t - t0)))^(1/nu), d, start = list(A = coef(m1)[["A"]], nu = nu0, t0 = coef(m1)[["t0"]]), lower = c(A = -3, nu = .05, t0 = 1850), upper = c(A = 1, nu = 10, t0 = 2076), control = ctl)
    if (!is.null(m) && (is.null(best) || sum(resid(m)^2) < sum(resid(best)^2))) best <- m }
  if (is.null(best)) NULL else c(as.list(coef(best)), B = B)
}
fit1 <- function(ttr, ytr, p0, P) {  # M1 for one period, start from M0 parameters
  w <- 2*pi/P; d <- data.frame(t = ttr, y = ytr, cw = cos(w*ttr), sw = sin(w*ttr), Bc = p0$B)
  m <- tf(y ~ A + (KF - A)/(1 + nu*exp(-Bc*(t - t0)))^(1/nu) + alpha*cw + beta*sw, d, start = list(A = p0$A, nu = p0$nu, t0 = p0$t0, alpha = .02, beta = .02), lower = c(A = -3, nu = .05, t0 = 1850, alpha = -2, beta = -2), upper = c(A = 1, nu = 10, t0 = 2076, alpha = 2, beta = 2), control = ctl)
  if (is.null(m)) NULL else list(cf = as.list(coef(m)), rss = sum(resid(m)^2), P = P)
}
skill <- function(y, t, Y) {
  tr <- t <= Y; ttr <- t[tr]; ytr <- y[tr]; tt <- t[!tr]; yt <- y[!tr]
  p0 <- fit0(ttr, ytr); if (is.null(p0)) return(c(NA, NA, NA))
  best <- NULL; for (P in 60:100) { f <- fit1(ttr, ytr, p0, P); if (!is.null(f) && (is.null(best) || f$rss < best$rss)) best <- f }
  if (is.null(best)) return(c(NA, NA, NA))
  f0 <- rich(tt, p0$A, p0$nu, p0$t0, p0$B); cf <- best$cf; w <- 2*pi/best$P
  f1 <- rich(tt, cf$A, cf$nu, cf$t0, p0$B) + cf$alpha*cos(w*tt) + cf$beta*sin(w*tt)
  c(1 - sum((yt - f1)^2)/sum((yt - f0)^2), best$P, sqrt(mean((yt - f0)^2)))
}
res <- list()
for (nm in names(S)) { y <- as.numeric(S[[nm]]); t <- as.numeric(names(S[[nm]])); n <- length(y)
  for (Y in cuts) {
    ob <- skill(y, t, Y)
    if (is.na(ob[1])) { res[[length(res)+1]] <- data.frame(series = nm, cut = Y, P_train = NA, skill = NA, null95 = NA, p = NA, n_null = 0, note = "not fittable"); print(res[[length(res)]], row.names = FALSE); next }
    tr <- t <= Y; p0 <- fit0(t[tr], y[tr]); mu <- rich(t, p0$A, p0$nu, p0$t0, p0$B); r <- y[tr] - mu[tr]; phi <- sum(r[-1]*r[-length(r)])/sum(r^2); sdi <- sd(r)*sqrt(1 - phi^2)
    nul <- unlist(mclapply(1:NS, function(i) { e <- as.numeric(stats::filter(rnorm(n + 100, 0, sdi), phi, method = "recursive"))[-(1:100)]; skill(mu + e, t, Y)[1] }, mc.cores = 8))
    nul <- nul[!is.na(nul)]
    res[[length(res)+1]] <- data.frame(series = nm, cut = Y, P_train = ob[2], skill = round(ob[1], 3), null95 = round(quantile(nul, .95), 3), p = round(mean(nul >= ob[1]), 3), n_null = length(nul), note = "")
    print(res[[length(res)]], row.names = FALSE) } }
R <- do.call(rbind, res); write.csv(R, "out/oos_richards_skill.csv", row.names = FALSE)
cat("\nharmonic has forecast value with Richards trend (skill > null95 at >= 2 of 3 cut-offs):\n"); print(tapply(!is.na(R$skill) & R$skill > R$null95, R$series, sum) >= 2)
