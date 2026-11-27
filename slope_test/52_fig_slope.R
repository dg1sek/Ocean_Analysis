# 52_fig_slope.R -- figure: series with quadratic trend + sinusoid, and the slope split into trend and sinusoid parts (extended to 2060 for MED/GOM).
source("R/alf_tools.R")
rg <- read.csv("data/regional_monthly_ersst.csv"); gl <- read.csv("data/global_series.csv")
ann <- function(df, col) tapply(df[[col]], df$year, mean)
S <- list(Global = ann(gl, "ersst"), `Mittelmeer (MED)` = ann(rg, "MED"), `Golf von Mexiko (GOM)` = ann(rg, "GOM"))
png("out/slope_fig.png", 2400, 1500, res = 200); par(mfcol = c(2, 3), mar = c(4, 4.2, 2.5, 1), cex = 0.9, cex.main = 0.95)
for (nm in names(S)) {
  y <- as.numeric(S[[nm]]); t <- as.numeric(names(S[[nm]])); n <- length(y)
  sc <- sin_scan(y, t, 60:100, 2); P <- sc$period[which.max(sc$dAIC)]; phi <- sc$phi[1]
  X <- function(tt) { x <- (tt - 1988)/100; cbind(1, x, x^2, sin(2*pi*tt/P), cos(2*pi*tt/P)) }
  Xm <- X(t); pw <- function(M) rbind(M[1, , drop = FALSE]*sqrt(1-phi^2), M[-1, , drop = FALSE] - phi*M[-n, , drop = FALSE])
  cf <- lm.fit(pw(Xm), pw(matrix(y)))$coefficients
  tt <- seq(1850, 2060, .25); Xt <- X(tt)
  tr <- as.numeric(Xt[, 1:3] %*% cf[1:3]); si <- as.numeric(Xt[, 4:5] %*% cf[4:5])
  plot(t, y, pch = 16, cex = .5, col = "grey50", xlim = c(1850, 2060), xlab = "Jahr", ylab = "Anomalie (K)", main = sprintf("%s, Periode %d J.", nm, P))
  lines(tt, tr, col = "#1f77b4", lwd = 2); lines(tt, tr + si, col = "#d62728", lwd = 2); abline(v = 2025, lty = 3)
  if (nm == names(S)[1]) legend("topleft", c("Jahresmittel ERSSTv6", "Trend (quadratisch)", "Trend + Sinus"), col = c("grey50", "#1f77b4", "#d62728"), pch = c(16, NA, NA), lwd = c(NA, 2, 2), bty = "n", cex = .8)
  x <- (tt - 1988)/100; st <- 10*(cf[2] + 2*cf[3]*x)/100
  ss <- 10*diff(si)/diff(tt); ss <- c(ss, tail(ss, 1))
  plot(tt, st + ss, type = "n", xlab = "Jahr", ylab = "Steigung (K/Jahrzehnt)", ylim = range(c(st, ss, st + ss, 0)), main = "Steigung")
  abline(h = 0, col = "grey70"); abline(v = 2025, lty = 3)
  lines(tt, st, col = "#1f77b4", lwd = 2); lines(tt, ss, col = "#2ca02c", lwd = 2); lines(tt, st + ss, col = "#d62728", lwd = 2)
  if (nm == names(S)[1]) legend("topleft", c("Trend", "Sinus", "Summe"), col = c("#1f77b4", "#2ca02c", "#d62728"), lwd = 2, bty = "n", cex = .8)
}
dev.off()
