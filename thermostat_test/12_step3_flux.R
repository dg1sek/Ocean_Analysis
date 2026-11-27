# Tropical thermostat, step 3b: sign of the surface-flux feedback in warm vs cool tropical ocean cells.
# Design fixed in PREREG_step2_3.md (section 3b). Fluxes: NCEP/NCAR Reanalysis 1 (a model product, caveat!); SST: ERSSTv6 (nearest 2-degree cell).
# Heat flux into the ocean Q = (SW down - SW up) + (LW down - LW up) - latent - sensible (W/m2, positive = ocean heats).
# Run from this folder: Rscript 12_step3_flux.R   (needs data/ncep_{dswrf,uswrf,dlwrf,ulwrf,lhtfl,shtfl}.nc and ERSSTv6)
suppressMessages({library(terra); library(ncdf4)})
fe <- function(...) { p <- c(...); p[file.exists(p)][1] }
source("R/read_nc.R"); dir.create("out", showWarnings = FALSE); set.seed(20261127)
Y0 <- 1982; Y1 <- 2025; NB <- 1000; BL <- 3
rd <- function(v) { nc <- nc_open(sprintf("data/ncep_%s.nc", v)); a <- ncvar_get(nc, v); lon <- nc$dim$lon$vals; lat <- nc$dim$lat$vals
  tt <- as.Date("1800-01-01") + nc$dim$time$vals / 24; nc_close(nc); list(a = a, lon = lon, lat = lat, tt = tt) }
F <- lapply(c(dswrf = "dswrf", uswrf = "uswrf", dlwrf = "dlwrf", ulwrf = "ulwrf", lhtfl = "lhtfl", shtfl = "shtfl"), rd)
lon <- F$dswrf$lon; lat <- F$dswrf$lat; tt <- F$dswrf$tt; yr <- as.integer(format(tt, "%Y")); mo <- as.integer(format(tt, "%m"))
cat(sprintf("NCEP R1 grid %d x %d, months %s to %s\n", length(lon), length(lat), format(min(tt)), format(max(tt))))
sel <- which(yr >= Y0 & yr <= Y1); nm <- length(sel); cat("months used:", nm, "\n")
nx <- length(lon); ny <- length(lat); asm <- function(v) matrix(F[[v]]$a[, , sel], nx * ny, nm)
SW <- asm("dswrf") - asm("uswrf"); LW <- asm("dlwrf") - asm("ulwrf"); LH <- asm("lhtfl"); SH <- asm("shtfl"); Q <- SW + LW - LH - SH
latv <- rep(lat, each = nx); lonv <- rep(lon, times = ny)
# SST: ERSSTv6, nearest 2-degree cell
E <- fe("data/ersst.v6.sst.mnmean.nc", "../regional/data/ersst.v6.sst.mnmean.nc", "../followup/data/ersst.v6.sst.mnmean.nc", "../regional_test/data/ersst.v6.sst.mnmean.nc")
r <- tryCatch(rast(E, subds = "sst"), error = function(e) tryCatch(rast(E), error = function(e2) rast_nc(E, "sst")))
xy <- xyFromCell(r, 1:ncell(r)); elon <- (xy[, 1] + 360) %% 360; elat <- xy[, 2]
lonu <- sort(unique(round(elon, 3))); latu <- sort(unique(round(elat, 3)))
id <- function(a, b) which(abs(elon - a) < 1e-3 & abs(elat - b) < 1e-3)[1]
ilon <- sapply(lonv, function(x) lonu[which.min(pmin(abs(lonu - x), 360 - abs(lonu - x)))]); ilat <- sapply(latv, function(x) latu[which.min(abs(latu - x))])
key <- match(paste(round(ilon, 3), round(ilat, 3)), paste(round(elon, 3), round(elat, 3)))
ntE <- nlyr(r); tmE <- seq(as.Date("1850-01-01"), by = "month", length.out = ntE); yE <- as.integer(format(tmE, "%Y"))
SSTall <- values(r); SSTall[SSTall < -5 | SSTall > 45] <- NA
selE <- which(yE >= Y0 & yE <= Y1); SST <- SSTall[key, selE]; stopifnot(ncol(SST) == nm)
rm(SSTall)
# cell classes from the 1991-2020 climatology of SST (annual mean), tropics |lat|<=30, ocean with complete SST in the period
clim_idx <- which(yr[sel] >= 1991 & yr[sel] <= 2020); sstm <- rowMeans(SST[, clim_idx], na.rm = TRUE)
valid <- rowSums(is.na(SST)) == 0 & abs(latv) <= 30
grp <- ifelse(!valid, NA, ifelse(sstm >= 28, "warm(>=28C)", ifelse(sstm < 26, "cool(<26C)", "mid(26-28C)")))
w <- cos(latv * pi / 180)
cat("cells:", paste(names(table(grp)), table(grp), collapse = "; "), "\n")
cat(sprintf("Sanity: mean Q over valid tropical ocean cells %.1f W/m2 (SW %.1f, LW %.1f, LH %.1f, SH %.1f)\n", mean(Q[valid, ]), mean(SW[valid, ]), mean(LW[valid, ]), mean(LH[valid, ]), mean(SH[valid, ])))
# monthly anomalies (calendar-month climatology over the period) and linear detrending per cell
anom <- function(M) { M <- M[valid, , drop = FALSE]; for (m in 1:12) { k <- mo[sel] == m; M[, k] <- M[, k] - rowMeans(M[, k, drop = FALSE]) }
  tc <- seq_len(nm) - mean(seq_len(nm)); M - outer(as.numeric(M %*% tc) / sum(tc^2), tc) }
A <- list(SST = anom(SST), Q = anom(Q), SW = anom(SW), LH = anom(-LH), LW = anom(LW), SH = anom(-SH))   # LH, SH entered as heating of the ocean (minus sign)
g <- grp[valid]; wv <- w[valid]; yy <- yr[sel]; yrs <- Y0:Y1
groups <- c("warm(>=28C)", "mid(26-28C)", "cool(<26C)")
# per-year sufficient statistics for the pooled slope: Sxy[year, group], Sxx
stat <- function(var) { out <- list()
  for (gg in groups) { k <- which(g == gg); x <- A$SST[k, , drop = FALSE]; y <- A[[var]][k, , drop = FALSE]; ww <- wv[k]
    sxy <- sapply(yrs, function(Y) { m <- yy == Y; sum(ww * rowSums(x[, m, drop = FALSE] * y[, m, drop = FALSE])) }); sxx <- sapply(yrs, function(Y) { m <- yy == Y; sum(ww * rowSums(x[, m, drop = FALSE]^2)) })
    out[[gg]] <- cbind(sxy, sxx) }
  out }
block_idx <- function(n, bl) { nb <- ceiling(n / bl); st <- sample.int(n - bl + 1, nb, replace = TRUE); (rep(st, each = bl) + rep(0:(bl - 1), nb))[1:n] }
ci <- function(x) quantile(x, c(.025, .975)); res <- list()
for (var in c("Q", "SW", "LH", "LW", "SH")) {
  S <- stat(var); sl <- function(ix) sapply(groups, function(gg) sum(S[[gg]][ix, 1]) / sum(S[[gg]][ix, 2]))
  obs <- sl(seq_along(yrs)); bs <- t(replicate(NB, sl(block_idx(length(yrs), BL))))
  d_wc <- obs["warm(>=28C)"] - obs["cool(<26C)"]; bd <- bs[, 1] - bs[, 3]
  res[[var]] <- data.frame(flux = var, group = groups, slope_W_m2_per_K = as.numeric(obs), lo = apply(bs, 2, function(z) ci(z)[1]), hi = apply(bs, 2, function(z) ci(z)[2]))
  res[[paste0(var, "_diff")]] <- data.frame(flux = var, group = "warm - cool", slope_W_m2_per_K = as.numeric(d_wc), lo = ci(bd)[1], hi = ci(bd)[2])
}
R <- do.call(rbind, res); rownames(R) <- NULL; write.csv(R, "out/step3_flux_slopes.csv", row.names = FALSE)
cat("\nPooled slope of the flux anomaly (heating of the ocean, W/m2) on the SST anomaly (K), detrended monthly anomalies, 1982-2025;\nblock bootstrap (3-yr blocks, 1000) 95 % intervals. Negative = damping of SST anomalies.\n"); print(R, digits = 3, row.names = FALSE)
cat("\nDone.\n")
