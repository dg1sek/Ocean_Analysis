# 60_eof.R -- EOF/PCA of annual ERSSTv6 SST (see PREREG_eof.md). Needs data/ersst.v6.sst.mnmean.nc (or a path found below).
suppressMessages(library(terra)); source("R/alf_tools.R"); source("R/read_nc.R"); set.seed(2026)
fe <- function(...) { p <- c(...); p[file.exists(p)][1] }
E <- fe("data/ersst.v6.sst.mnmean.nc", "../followup/data/ersst.v6.sst.mnmean.nc", "../thermostat/data/ersst.v6.sst.mnmean.nc", "../regional_test/data/ersst.v6.sst.mnmean.nc"); stopifnot(!is.na(E))
r <- tryCatch(rast(E, subds = "sst"), error = function(e) tryCatch(rast(E), error = function(e2) rast_nc(E, "sst")))
n <- nlyr(r); tm <- seq(as.Date("1850-01-01"), by = "month", length.out = n); yr <- as.integer(format(tm, "%Y")); mo <- as.integer(format(tm, "%m"))
S <- values(r); S[S < -5 | S > 45] <- NA; xy <- xyFromCell(r, 1:ncell(r)); inb <- yr >= 1991 & yr <= 2020
A <- S; for (m in 1:12) A[, mo == m] <- S[, mo == m] - rowMeans(S[, inb & mo == m], na.rm = TRUE)
years <- sort(unique(yr)); years <- years[years <= 2025]
Y <- sapply(years, function(y) rowMeans(A[, yr == y, drop = FALSE]))      # cell x year (NA if any month missing -> NA)
nino <- tapply(read.csv("data/regional_monthly_ersst.csv")$NINO34, read.csv("data/regional_monthly_ersst.csv")$year, mean)

run <- function(y0, remove_mean, tag) {
  k <- years >= y0; yy <- years[k]; X <- Y[, k]; ok <- which(rowSums(is.na(X)) == 0); X <- X[ok, ]; lat <- xy[ok, 2]; lon <- xy[ok, 1]
  w <- sqrt(cos(lat*pi/180)); if (remove_mean) X <- sweep(X, 2, colSums(X * w^2) / sum(w^2))
  X <- X - rowMeans(X); Z <- t(X * w)                                      # year x cell
  sv <- svd(Z, nu = 8, nv = 8); lam <- sv$d^2; frac <- lam / sum(lam); N <- length(yy)
  pcs <- sv$u[, 1:8] * rep(sv$d[1:8], each = N); sdp <- apply(pcs, 2, sd); pcs <- sweep(pcs, 2, sdp, "/")
  eofs <- X %*% pcs / (N - 1)                                              # pattern: K per 1 sd of the PC (regression)
  phi <- apply(pcs, 2, function(p) sum(p[-1]*p[-N]) / sum(p^2)); Ns <- N * (1 - phi^2) / (1 + phi^2)
  sep <- sapply(1:7, function(i) lam[i] - lam[i+1] > lam[i] * sqrt(2 / Ns[i]))
  cat(sprintf("\n== %s: %d cells, %d yr (%d-2025) ==\n", tag, length(ok), N, y0))
  cat("variance %:", round(100*frac[1:8], 1), "\nseparable from next (North):", sep, "\n")
  sc <- lapply(1:5, function(i) { s <- sin_scan(pcs[, i], yy, 60:100, 2); j <- which.max(s$dAIC); nl <- null_scan(pcs[, i], yy, 60:100, 2, nsim = 300)
    data.frame(tag = tag, pc = i, var_pct = round(100*frac[i], 1), separable = sep[i], best_P = s$period[j], dAIC = round(s$dAIC[j], 1), null95 = round(quantile(nl$best, .95), 1), p = round(mean(nl$best >= s$dAIC[j]), 3),
      r_nino = round(cor(pcs[, i], nino[as.character(yy)]), 2), phi = round(phi[i], 2)) })
  sc <- do.call(rbind, sc); print(sc, row.names = FALSE)
  # E5 lagged cross-correlations
  pairs <- t(combn(1:4, 2)); lagcc <- do.call(rbind, lapply(1:nrow(pairs), function(q) { a <- pcs[, pairs[q, 1]]; b <- pcs[, pairs[q, 2]]
    cc <- sapply(-15:15, function(L) if (L >= 0) cor(a[1:(N-L)], b[(1+L):N]) else cor(a[(1-L):N], b[1:(N+L)])); j <- which.max(abs(cc))
    data.frame(tag = tag, pair = sprintf("PC%d-PC%d", pairs[q,1], pairs[q,2]), lag_of_max = (-15:15)[j], r_max = round(cc[j], 2), r_lag0 = round(cc[16], 2)) })); print(lagcc, row.names = FALSE)
  list(sc = sc, lagcc = lagcc, pcs = pcs, yy = yy, eofs = eofs, lon = lon, lat = lat, frac = frac, ok = ok, X = X, w = w)
}
res <- list(); for (v in list(list(1900, FALSE, "V1_1900"), list(1900, TRUE, "V2_1900_noglobalmean"), list(1850, FALSE, "V1_1850"), list(1850, TRUE, "V2_1850_noglobalmean"))) res[[v[[3]]]] <- run(v[[1]], v[[2]], v[[3]])
write.csv(do.call(rbind, lapply(res, `[[`, "sc")), "out/eof_modes.csv", row.names = FALSE)
write.csv(do.call(rbind, lapply(res, `[[`, "lagcc")), "out/eof_lagcc.csv", row.names = FALSE)
for (nm in names(res)) write.csv(data.frame(year = res[[nm]]$yy, round(res[[nm]]$pcs[, 1:5], 4)), sprintf("out/eof_pcs_%s.csv", nm), row.names = FALSE)

# E3: where does the line live? per-cell sinusoid amplitude, 1900-2025
rr <- res[["V1_1900"]]; best <- rr$sc$best_P[rr$sc$p < .05]; Pst <- if (length(best)) best[1] else 76
cat("\nE3: period used for the map:", Pst, "yr", if (!length(best)) "(no significant PC; default 76)" else "", "\n")
k <- years >= 1900; yy <- years[k]; X <- Y[, k]; ok <- rr$ok; X <- X[ok, ]; x <- (yy - 1988)/100
M <- cbind(1, x, x^2, sin(2*pi*yy/Pst), cos(2*pi*yy/Pst)); cf <- lm.fit(M, t(X))$coefficients; amp <- sqrt(cf[4, ]^2 + cf[5, ]^2)
write.csv(data.frame(lon = rr$lon, lat = rr$lat, amp = amp), "out/eof_line_amp_map.csv", row.names = FALSE)

# figures
mp <- function(lon, lat, z, main, zl) { zz <- pmin(pmax(z, -zl), zl); plot(lon, lat, pch = 15, cex = 0.45, col = hcl.colors(101, "Blue-Red 3")[round((zz + zl)/(2*zl)*100) + 1], main = main, xlab = "", ylab = "", asp = 1, axes = FALSE); box() }
png("out/eof_fig.png", 2600, 2000, res = 200); par(mfrow = c(4, 2), mar = c(1, 1, 2.3, 1))
for (nm in c("V1_1900", "V2_1900_noglobalmean")) { rs <- res[[nm]]; for (i in 1:2) {
  mp(rs$lon, rs$lat, rs$eofs[, i] * 1, sprintf("%s  EOF%d  (%.0f%% var)", nm, i, 100*rs$frac[i]), max(abs(rs$eofs[, i]))) } }
for (nm in c("V1_1900", "V2_1900_noglobalmean")) { rs <- res[[nm]]; for (i in 3:4) {
  mp(rs$lon, rs$lat, rs$eofs[, i] * 1, sprintf("%s  EOF%d  (%.0f%% var)", nm, i, 100*rs$frac[i]), max(abs(rs$eofs[, i]))) } }
dev.off()
png("out/eof_pcs_fig.png", 2400, 1400, res = 200); par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
for (nm in c("V1_1900", "V2_1900_noglobalmean")) { rs <- res[[nm]]; matplot(rs$yy, rs$pcs[, 1:4], type = "l", lty = 1, lwd = 1.5, col = c("black", "#d62728", "#1f77b4", "#2ca02c"), main = paste(nm, "PC1-4 (unit sd)"), xlab = "", ylab = "") ; legend("topleft", paste0("PC", 1:4), col = c("black", "#d62728", "#1f77b4", "#2ca02c"), lwd = 2, horiz = TRUE, bty = "n") }
dev.off()
png("out/eof_line_map.png", 2200, 1100, res = 200); par(mar = c(1, 1, 2.3, 1)); zl <- quantile(amp, .98)
plot(rr$lon, rr$lat, pch = 15, cex = 0.45, col = hcl.colors(101, "YlOrRd", rev = TRUE)[round(pmin(amp, zl)/zl*100) + 1], asp = 1, axes = FALSE, xlab = "", ylab = "", main = sprintf("Amplitude (K) of a %d-yr sinusoid per 2-deg cell, 1900-2025 (scale max %.2f)", Pst, zl)); box(); dev.off()
