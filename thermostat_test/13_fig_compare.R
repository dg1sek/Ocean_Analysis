# Comparison figure for step 2: three SST products on the common period 1982-2022 (reads out/step2_bins.csv, out/step2_products.csv)
b <- read.csv("out/step2_bins.csv"); p <- read.csv("out/step2_products.csv")
b <- b[b$period == "P5" & b$ncell >= 5, ]; p <- p[p$period == "P5", ]
b$mid <- as.numeric(sub("^[\\[(]([-0-9.]+),.*", "\\1", b$bin)) + 0.5
prods <- c("ERSSTv6", "OISSTv2", "HadSST4"); col <- c(ERSSTv6 = "#0072B2", OISSTv2 = "#D55E00", HadSST4 = "#009E73")
lab <- c(ERSSTv6 = "ERSSTv6 (2°, filled)", OISSTv2 = "OISST v2 (1°, satellite)", HadSST4 = "HadSST4.2 (5°, unfilled)")
png("out/step2_fig_compare.png", 1900, 800, res = 170); par(mar = c(4.6, 4.6, 3.2, 1))
layout(matrix(1:2, 1), widths = c(1.15, 1))
plot(NA, xlim = c(15, 31), ylim = c(-0.1, 0.32), xlab = "mean SST of the cell, 1982-2001 (°C)", ylab = "SST trend 1982-2022 (K per decade)",
     main = "A  Warming rate vs. warmth of the cell", cex.main = 1)
abline(v = 27, lty = 3, col = "grey60"); abline(h = 0, col = "grey85")
for (i in seq_along(prods)) { d <- b[b$product == prods[i], ]; off <- (i - 2) * 0.18
  arrows(d$mid + off, d$lo, d$mid + off, d$hi, 0.02, 90, 3, col = adjustcolor(col[i], 0.55)); lines(d$mid + off, d$trend, col = col[i], lwd = 1.5); points(d$mid + off, d$trend, pch = 19, col = col[i], cex = 0.8) }
legend("topleft", lab[prods], col = col[prods], pch = 19, lty = 1, bty = "n", cex = 0.85)
text(27.3, -0.075, "a thermostat would\nbend the curve down\nright of this line", adj = 0, cex = 0.7, col = "grey35")
# B: forest plot of the two scored tests
par(mar = c(4.6, 9.5, 3.2, 2)); plot(NA, xlim = c(-0.045, 0.085), ylim = c(0.5, 6.5), yaxt = "n", ylab = "", xlab = "K per decade (95 % block-bootstrap interval)", main = "B  Scored tests (< 0 = thermostat)", cex.main = 1)
abline(v = 0, lwd = 1.5, col = "grey40"); y <- 6
for (tst in c("A2", "B1")) { for (pr in prods) { r <- p[p$product == pr, ]; est <- r[[tst]]; lo <- r[[paste0(tst, "_lo")]]; hi <- r[[paste0(tst, "_hi")]]
    segments(lo, y, hi, y, col = col[pr], lwd = 3); points(est, y, pch = 19, col = col[pr], cex = 1.3); y <- y - 1 }; y <- y - 0.0 }
axis(2, at = c(5, 2), labels = c("A2: slope of trend\nvs starting SST\n(cells >= 27 °C)", "B1: trend of the 99 %\nquantile minus trend\nof the median"), las = 1, tick = FALSE, cex.axis = 0.75)

dev.off()
