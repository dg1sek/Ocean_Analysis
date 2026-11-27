# Tropical thermostat, step 2: the step-1 tests on three SST products. Design and decision rule: PREREG_step1.md and PREREG_step2_3.md.
# Run from this folder: Rscript 11_step2_products.R   (needs data/oisst.v2.sst.mnmean.nc, data/HadSST.4.2.0.0_actuals_median.nc, ERSSTv6)
suppressMessages({library(terra); library(ncdf4)})
fe <- function(...) { p <- c(...); p[file.exists(p)][1] }
source("R/read_nc.R"); dir.create("out", showWarnings = FALSE); set.seed(20261127)
NB <- 1000; BL <- 5; thrs <- c(27, 28, 29, 30)
# ---- loaders: each returns list(S = cells x months, lat, yr, mo) for |lat| <= 30 ----
load_ersst <- function() {
  f <- fe("data/ersst.v6.sst.mnmean.nc", "../regional/data/ersst.v6.sst.mnmean.nc", "../followup/data/ersst.v6.sst.mnmean.nc", "../regional_test/data/ersst.v6.sst.mnmean.nc")
  r <- tryCatch(rast(f, subds = "sst"), error = function(e) tryCatch(rast(f), error = function(e2) rast_nc(f, "sst")))
  n <- nlyr(r); tm <- seq(as.Date("1850-01-01"), by = "month", length.out = n); S <- values(r); S[S < -5 | S > 45] <- NA
  lat <- xyFromCell(r, 1:ncell(r))[, 2]; k <- abs(lat) <= 30
  list(S = S[k, ], lat = lat[k], yr = as.integer(format(tm, "%Y")), mo = as.integer(format(tm, "%m")))
}
load_nc <- function(f, var, tconv, mask = NULL) {
  nc <- nc_open(f); a <- ncvar_get(nc, var); if (!is.null(mask)) { m <- ncvar_get(nc_open(mask), "mask"); for (i in seq_len(dim(a)[3])) { z <- a[, , i]; z[m < 0.5] <- NA; a[, , i] <- z } }  # OISST v2 fills land: apply NOAA land-sea mask
  lat <- nc$dim$lat$vals %||% nc$dim$latitude$vals; tt <- tconv(nc); nc_close(nc)
  nx <- dim(a)[1]; ny <- dim(a)[2]; nt <- dim(a)[3]; S <- matrix(a, nx * ny, nt); latv <- rep(lat, each = nx)
  S[S < -5 | S > 45] <- NA; k <- abs(latv) <= 30
  list(S = S[k, ], lat = latv[k], yr = as.integer(format(tt, "%Y")), mo = as.integer(format(tt, "%m")))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
t_oisst <- function(nc) as.Date("1800-01-01") + nc$dim$time$vals          # days since 1800-01-01
t_had <- function(nc) { u <- nc$dim$time$units; as.Date(sub(".*since ([0-9-]+).*", "\\1", u)) + nc$dim$time$vals }
DATA <- list(
  ERSSTv6 = list(load = load_ersst, periods = c("P1", "P2", "P4", "P5")),
  OISSTv2 = list(load = function() load_nc("data/oisst.v2.sst.mnmean.nc", "sst", t_oisst, mask = "data/lsmask.nc"), periods = "P5"),
  HadSST4 = list(load = function() load_nc("data/HadSST.4.2.0.0_actuals_median.nc", "tos", t_had), periods = c("P1", "P2", "P4", "P5")))
PER <- list(P1 = c(1950, 2025), P2 = c(1980, 2025), P4 = c(1982, 2025), P5 = c(1982, 2022)); BASEN <- c(P1 = 30, P2 = 20, P4 = 20, P5 = 20)
wq <- function(x, w, p) { o <- order(x); cw <- cumsum(w[o]) / sum(w[o]); x[o][vapply(p, function(pp) which(cw >= pp)[1], 1L)] }
slope_vec <- function(Y, tt) { tc <- tt - mean(tt); as.numeric((Y - rowMeans(Y)) %*% tc) / sum(tc^2) }
slope1 <- function(y, tt) { tc <- tt - mean(tt); sum((y - mean(y)) * tc) / sum(tc^2) }
block_idx <- function(ny, bl) { nb <- ceiling(ny / bl); st <- sample.int(ny - bl + 1, nb, replace = TRUE); (rep(st, each = bl) + rep(0:(bl - 1), nb))[1:ny] }
ci <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
rows <- list(); bins <- list()
for (dn in names(DATA)) {
  D <- DATA[[dn]]$load(); cat(sprintf("\n######## %s: %d tropical cells, %d-%d\n", dn, nrow(D$S), min(D$yr), max(D$yr)))
  for (pn in DATA[[dn]]$periods) {
    y0 <- PER[[pn]][1]; y1 <- PER[[pn]][2]; yrs <- y0:y1; ny <- length(yrs)
    if (min(D$yr) > y0) { cat(pn, ": product starts later, skipped\n"); next }
    # annual means need 12 valid months; cell complete = every year valid
    A <- sapply(yrs, function(y) { m <- D$yr == y; v <- D$S[, m, drop = FALSE]; ifelse(rowSums(!is.na(v)) == 12, rowMeans(v), NA) })
    ok <- rowSums(is.na(A)) == 0; if (sum(ok) < 20) { cat(pn, ": only", sum(ok), "complete cells, skipped\n"); next }
    A <- A[ok, , drop = FALSE]; lat <- D$lat[ok]; w <- cos(lat * pi / 180); Sok <- D$S[ok, , drop = FALSE]
    cat(sprintf("\n== %s %s %d-%d: %d complete cells ==\n", dn, pn, y0, y1, nrow(A)))
    pred <- rowMeans(A[, 1:BASEN[pn]]); tr <- slope_vec(A, yrs) * 10; hi <- pred >= 27
    a2 <- function(trv) { ww <- w[hi]; x <- pred[hi]; yv <- trv[hi]; xm <- sum(ww * x) / sum(ww); sum(ww * (x - xm) * (yv - sum(ww * yv) / sum(ww))) / sum(ww * (x - xm)^2) }
    Q <- t(vapply(yrs, function(y) { v <- Sok[, D$yr == y, drop = FALSE]; ww <- rep(w, ncol(v)); v <- as.numeric(v); k <- !is.na(v)
           c(wq(v[k], ww[k], c(.5, .9, .99)), vapply(thrs, function(th) sum(ww[k][v[k] > th]) / sum(ww[k]), 1)) }, numeric(7)))
    colnames(Q) <- c("q50", "q90", "q99", paste0("frac", thrs)); bm <- cut(pred, seq(floor(min(pred)), ceiling(max(pred)), 1), include.lowest = TRUE)
    b1 <- (slope1(Q[, "q99"], yrs) - slope1(Q[, "q50"], yrs)) * 10; a2o <- a2(tr)
    bs_a2 <- numeric(NB); bs_b1 <- numeric(NB); bs_f30 <- numeric(NB); bs_bin <- matrix(NA, NB, nlevels(bm), dimnames = list(NULL, levels(bm)))
    for (b in 1:NB) { ix <- block_idx(ny, BL); tt <- yrs[ix]; trb <- slope_vec(A[, ix, drop = FALSE], tt) * 10; bs_a2[b] <- a2(trb)
      bs_b1[b] <- (slope1(Q[ix, "q99"], tt) - slope1(Q[ix, "q50"], tt)) * 10; bs_f30[b] <- slope1(Q[ix, "frac30"], tt) * 10
      bs_bin[b, ] <- tapply(trb * w, bm, sum) / tapply(w, bm, sum) }
    bo <- tapply(tr * w, bm, sum) / tapply(w, bm, sum); nb <- table(bm)
    bins[[paste(dn, pn)]] <- data.frame(product = dn, period = pn, bin = levels(bm), ncell = as.integer(nb), trend = as.numeric(bo), lo = apply(bs_bin, 2, function(z) ci(z)[1]), hi = apply(bs_bin, 2, function(z) ci(z)[2]))
    rows[[paste(dn, pn)]] <- data.frame(product = dn, period = pn, ncell = nrow(A), A2 = a2o, A2_lo = ci(bs_a2)[1], A2_hi = ci(bs_a2)[2],
      B1 = b1, B1_lo = ci(bs_b1)[1], B1_hi = ci(bs_b1)[2], frac30_start = Q[1, "frac30"], frac30_end = Q[ny, "frac30"], frac30_trend_dec = slope1(Q[, "frac30"], yrs) * 10,
      q99_start = Q[1, "q99"], q99_end = Q[ny, "q99"], q50_start = Q[1, "q50"], q50_end = Q[ny, "q50"], mean_trend = sum(w * tr) / sum(w))
    print(t(round(unlist(rows[[paste(dn, pn)]][-(1:2)]), 4)))
  }
}
R <- do.call(rbind, rows); write.csv(R, "out/step2_products.csv", row.names = FALSE); write.csv(do.call(rbind, bins), "out/step2_bins.csv", row.names = FALSE)
cat("\n\n===== SUMMARY (A2, B1 with 95 % block-bootstrap intervals) =====\n"); print(R[, c("product", "period", "ncell", "A2", "A2_lo", "A2_hi", "B1", "B1_lo", "B1_hi", "frac30_start", "frac30_end")], digits = 3, row.names = FALSE)
cat("\nDone.\n")
