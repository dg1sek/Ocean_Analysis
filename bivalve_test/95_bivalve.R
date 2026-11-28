# 95_bivalve.R -- shell chronology test (PREREG_bivalve.md). Run: Rscript 95_bivalve.R
suppressMessages(library(terra)); source("R/alf_tools.R"); source("R/read_nc.R"); set.seed(2026)
L <- readLines("data/butler2013-raw.txt", encoding = "latin1"); i0 <- grep("^age_AD", L)[1]; D <- read.delim(text = L[(i0 + 1):length(L)], header = FALSE, na.strings = "-999", stringsAsFactors = FALSE)
ns <- ncol(D) %/% 2; shells <- lapply(1:ns, function(k) { a <- D[[2*k - 1]]; v <- D[[2*k]]; ok <- !is.na(a) & !is.na(v) & a > 0; data.frame(shell = k, year = a[ok], inc = v[ok]) })
SH <- do.call(rbind, shells); SH <- SH[!is.na(SH$inc) & SH$inc > 0, ]; SH$age <- ave(SH$year, SH$shell, FUN = function(z) z - min(z) + 1)
cat(sprintf("%d shells, %d increments, years %d-%d\n", length(unique(SH$shell)), nrow(SH), min(SH$year), max(SH$year)))
# RCS: mean increment by biological age, smoothed (88-yr spline)
ma <- tapply(SH$inc, SH$age, mean); ages <- as.numeric(names(ma)); sp <- smooth.spline(ages, as.numeric(ma), spar = NULL, df = max(4, round(length(ages)/88)*2)); curve <- predict(sp, SH$age)$y
SH$idx <- SH$inc / curve
chr <- aggregate(idx ~ year, SH, function(z) c(m = mean(z), n = length(z))); chr <- data.frame(year = chr$year, idx = chr$idx[, "m"], n = chr$idx[, "n"]); chr <- chr[order(chr$year), ]
cat(sprintf("RCS chronology %d-%d; years with >= 3 shells from %d\n", min(chr$year), max(chr$year), min(chr$year[chr$n >= 3 & chr$year > 1000])))
write.csv(chr, "out/bivalve_rcs_chronology.csv", row.names = FALSE)
per <- 40:120
scan <- function(y, t, nm) { y <- (y - mean(y))/sd(y); s <- sin_scan(y, t, per, 1); j <- which.max(s$dAIC); nl <- null_scan(y, t, per, 1, nsim = 300)
  data.frame(series = nm, from = min(t), to = max(t), n = length(y), best_P = s$period[j], amp_sd = round(s$amp[j], 2), dAIC = round(s$dAIC[j], 1), null95 = round(quantile(nl$best, .95), 1), p = round(mean(nl$best >= s$dAIC[j]), 3), phi = round(s$phi[1], 2)) }
seg <- function(a, b) chr[chr$year >= a & chr$year <= b, ]
B1 <- list(); for (w in list(c(1175, 2005, "main"), c(1175, 1590, "half 1"), c(1590, 2005, "half 2"), c(649, 2005, "full (1-2 shells before 1175)"))) { d <- seg(as.numeric(w[1]), as.numeric(w[2])); B1[[length(B1)+1]] <- scan(d$idx, d$year, paste("RCS own,", w[3])) }
bc <- read.delim(text = { l <- readLines("data/butler2013-chron.txt", encoding = "latin1"); l[(grep("^age_AD", l)[1]):length(l)] }, stringsAsFactors = FALSE)
for (nm in c("NE.ars", "RCS.ars")) B1[[length(B1)+1]] <- scan(bc[[nm]], bc$age_AD, paste("published", nm, "(pre-whitened), 649-2005"))
B1 <- do.call(rbind, B1); cat("\n== B1/B2: best period 40-120 yr ==\n"); print(B1, row.names = FALSE); write.csv(B1, "out/bivalve_B1_line.csv", row.names = FALSE)
# top peaks of the main scan
d <- seg(1175, 2005); s <- sin_scan((d$idx - mean(d$idx))/sd(d$idx), d$year, per, 1); top <- s[order(-s$dAIC), ]; pk <- NULL; for (r in seq_len(nrow(top))) { if (all(abs(top$period[r] - pk) > 8)) pk <- c(pk, top$period[r]); if (length(pk) == 4) break }
cat("\nmain series: four strongest separate peaks (period, dAIC):\n"); print(round(s[match(pk, s$period), c("period", "dAIC")], 1), row.names = FALSE)
# B3: SST in a box around Grimsey
E <- c("data/ersst.v6.sst.mnmean.nc", "../followup/data/ersst.v6.sst.mnmean.nc")[file.exists(c("data/ersst.v6.sst.mnmean.nc", "../followup/data/ersst.v6.sst.mnmean.nc"))][1]; r <- tryCatch(rast(E, subds = "sst"), error = function(e) tryCatch(rast(E), error = function(e2) rast_nc(E, "sst")))
n <- nlyr(r); tm <- seq(as.Date("1850-01-01"), by = "month", length.out = n); yr <- as.integer(format(tm, "%Y")); mo <- as.integer(format(tm, "%m")); S <- values(r); S[S < -5 | S > 45] <- NA
xy <- xyFromCell(r, 1:ncell(r)); k <- xy[, 2] >= 62 & xy[, 2] <= 70 & (xy[, 1] >= 332 | xy[, 1] <= 350) & xy[, 1] >= 332 & xy[, 1] <= 350
w <- cos(xy[k, 2]*pi/180); X <- S[k, , drop = FALSE]; ts <- colSums(X * w, na.rm = TRUE)/colSums((!is.na(X)) * w); ts[colSums(!is.na(X)) < 0.5*sum(k)] <- NA
clim <- sapply(1:12, function(m) mean(ts[mo == m & yr >= 1951 & yr <= 1980], na.rm = TRUE)); an <- ts - clim[mo]
ann <- tapply(an, yr, function(z) if (sum(!is.na(z)) >= 10) mean(z, na.rm = TRUE) else NA); ann <- ann[!is.na(ann) & as.integer(names(ann)) <= 2005 & as.integer(names(ann)) >= 1854]
j <- merge(data.frame(year = as.integer(names(ann)), sst = as.numeric(ann)), chr, by = "year"); rr <- cor(j$sst, j$idx)
a1 <- function(z) sum(z[-1]*z[-length(z)])/sum(z^2); ne <- nrow(j)*(1 - a1(j$sst - mean(j$sst))*a1(j$idx - mean(j$idx)))/(1 + a1(j$sst - mean(j$sst))*a1(j$idx - mean(j$idx)))
tt <- rr*sqrt((ne - 2)/(1 - rr^2)); cat(sprintf("\nB3: shell index vs ERSST box 62-70N,28-10W annual (n = %d, 1854-2005): r = %.2f, effective n = %.0f, p = %.3f\n", nrow(j), rr, ne, 2*pt(-abs(tt), ne - 2)))
write.csv(data.frame(r = rr, n = nrow(j), n_eff = ne), "out/bivalve_B3_sst.csv", row.names = FALSE)
png("out/bivalve_fig.png", 2400, 1500, res = 200); par(mfrow = c(2, 1), mar = c(4, 4, 2.5, 1))
plot(chr$year, chr$idx, type = "l", col = "grey60", xlab = "Jahr", ylab = "Wachstumsindex (RCS)", main = "Arctica islandica, Nordisländischer Schelf (eigene RCS-Chronologie)"); lines(chr$year, stats::filter(chr$idx, rep(1/31, 31), sides = 2), col = "#d62728", lwd = 2); abline(v = 1175, lty = 3)
pg <- spec.pgram(ts(d$idx - mean(d$idx)), taper = 0.1, plot = FALSE); plot(1/pg$freq, pg$spec, type = "l", log = "y", xlim = c(10, 200), xlab = "Periode (Jahre)", ylab = "Spektrum", main = "Periodogramm 1175-2005"); abline(v = c(65, 80), lty = 3); dev.off()
