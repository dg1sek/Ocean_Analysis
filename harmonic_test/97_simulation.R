# 97_simulation.R -- coupled stochastic oscillators (PREREG_sim.md). Run: Rscript 97_simulation.R
suppressMessages(library(parallel)); source("R/alf_tools.R"); RNGkind("L'Ecuyer-CMRG"); set.seed(2026)
t <- 1850:2025; n <- length(t); per <- 60:100
rg <- read.csv("data/regional_monthly_ersst.csv"); rg <- rg[rg$year <= 2025, ]; ymed <- as.numeric(tapply(rg$MED, rg$year, mean))
X0 <- cbind(1, (t - 1988)/100, ((t - 1988)/100)^2); f0 <- lm.fit(X0, ymed); trend <- as.numeric(X0 %*% f0$coefficients); sd_res <- sd(f0$residuals)
cat(sprintf("real MED: residual sd about the quadratic = %.3f K\n", sd_res))
ar2sim <- function(P, Q, n) { w <- 2*pi/P; r <- exp(-w/(2*Q)); x <- as.numeric(stats::filter(rnorm(n + 300), c(2*r*cos(w), -r^2), method = "recursive"))[-(1:300)]; x/sd(x) }
draw_period <- function(scen) { repeat { p <- exp(runif(1, log(10), log(150))); if (scen != "S1a" || p < 55 || p > 110) return(p) } }
make_series <- function(scen) {
  if (scen == "S0") { e <- as.numeric(stats::filter(rnorm(n + 100), 0.3, method = "recursive"))[-(1:100)]; return(trend + sd_res*e/sd(e)) }
  Ps <- replicate(3, draw_period(scen)); if (scen == "S2") Ps[1] <- runif(1, 65, 85)
  Qs <- runif(3, 2, 10); f <- runif(1, .3, .7); osc <- rowSums(sapply(1:3, function(k) ar2sim(Ps[k], Qs[k], n)))/sqrt(3)
  bg <- as.numeric(stats::filter(rnorm(n + 100), 0.3, method = "recursive"))[-(1:100)]; bg <- bg/sd(bg)
  trend + sd_res*(sqrt(f)*osc/sd(osc) + sqrt(1 - f)*bg) }
scanb <- function(y) { s <- sin_scan(y, t, per, 2); j <- which.max(s$dAIC); c(P = s$period[j], g = s$dAIC[j]) }
null_max <- function(y, kind, nsim = 100) {
  f <- lm.fit(X0, y); r <- f$residuals; mu <- as.numeric(X0 %*% f$coefficients)
  if (kind == "AR1") { ph <- sum(r[-1]*r[-n])/sum(r^2); sdi <- sd(r)*sqrt(1 - ph^2); gen <- function() as.numeric(stats::filter(rnorm(n + 100, 0, sdi), ph, method = "recursive"))[-(1:100)] }
  else { a <- ar(r, order.max = 2, aic = FALSE, method = "yule-walker"); gen <- function() as.numeric(stats::filter(rnorm(n + 100, 0, sqrt(a$var.pred)), a$ar, method = "recursive"))[-(1:100)] }
  vapply(1:nsim, function(i) scanb(mu + gen())[["g"]], 0) }
design <- function(tt, P = NULL) { x <- (tt - 1988)/100; cbind(1, x, x^2, if (!is.null(P)) cbind(sin(2*pi*tt/P), cos(2*pi*tt/P))) }
skills <- function(y) sapply(c(1950, 1975, 2000), function(Y) { tr <- t <= Y; s <- sin_scan(y[tr], t[tr], per, 2); P <- s$period[which.max(s$dAIC)]
  c0 <- lm.fit(design(t[tr]), y[tr])$coefficients; c1 <- lm.fit(design(t[tr], P), y[tr])$coefficients; yf <- y[!tr]
  1 - sum((yf - design(t[!tr], P) %*% c1)^2)/sum((yf - design(t[!tr]) %*% c0)^2) })
thr <- c(0.53, 0.50, 0.39)
one <- function(scen) { y <- make_series(scen); o <- scanb(y); p1 <- mean(null_max(y, "AR1") >= o[["g"]]); p2 <- mean(null_max(y, "AR2") >= o[["g"]]); sk <- skills(y)
  c(P = o[["P"]], gain = o[["g"]], p_ar1 = p1, p_ar2 = p2, pass_fc = as.numeric(sum(sk > thr) >= 2), sk1 = sk[1], sk2 = sk[2], sk3 = sk[3]) }
NS <- 400; res <- list()
for (scen in c("S0", "S1a", "S1b", "S2")) { m <- do.call(rbind, mclapply(1:NS, function(i) one(scen), mc.cores = 8)); m <- as.data.frame(m); m$scen <- scen; res[[scen]] <- m
  fl1 <- m$P >= 60 & m$P <= 90 & m$p_ar1 < .05; fl2 <- m$P >= 60 & m$P <= 90 & m$p_ar2 < .05
  cat(sprintf("%-4s n = %d | flagged (60-90 yr, p<0.05): AR(1) null %.1f %%, AR(2) null %.1f %% | median AIC gain %.1f (90th pct %.1f) | forecast rule passes %.1f %%\n", scen, NS, 100*mean(fl1), 100*mean(fl2), median(m$gain), quantile(m$gain, .9), 100*mean(m$pass_fc))) }
R <- do.call(rbind, res); write.csv(R, "out/sim_results.csv", row.names = FALSE)
# real Mediterranean series
o <- scanb(ymed); pa1 <- mean(null_max(ymed, "AR1", 300) >= o[["g"]]); pa2 <- mean(null_max(ymed, "AR2", 300) >= o[["g"]])
cat(sprintf("\nreal MED (1850-2025): best period %d yr, AIC gain %.1f; p (AR(1) null) = %.3f, p (AR(2) null) = %.3f\n", o[["P"]], o[["g"]], pa1, pa2))
write.csv(data.frame(P = o[["P"]], gain = o[["g"]], p_ar1 = pa1, p_ar2 = pa2), "out/sim_real_med.csv", row.names = FALSE)
