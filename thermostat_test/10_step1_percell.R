# Tropical thermostat, step 1: does ERSSTv6 show a warm-water ceiling signature?  See PREREG_step1.md (design fixed beforehand).
# Run from this folder: Rscript 10_step1_percell.R
suppressMessages(library(terra))
fe <- function(...) { p <- c(...); p[file.exists(p)][1] }
source("R/read_nc.R")
ERSST <- fe("data/ersst.v6.sst.mnmean.nc", "../regional/data/ersst.v6.sst.mnmean.nc", "../followup/data/ersst.v6.sst.mnmean.nc",
            "../regional_test/data/ersst.v6.sst.mnmean.nc")
if (is.na(ERSST)) stop("ERSSTv6 file not found; run regional_test/00_download_ersst.R")
dir.create("out", showWarnings = FALSE); set.seed(20261126)
r <- tryCatch(rast(ERSST, subds = "sst"), error = function(e) tryCatch(rast(ERSST), error = function(e2) rast_nc(ERSST, "sst")))
n <- nlyr(r); tm <- seq(as.Date("1850-01-01"), by = "month", length.out = n)
yr <- as.integer(format(tm, "%Y")); mo <- as.integer(format(tm, "%m"))
S <- values(r); S[S < -5 | S > 45] <- NA
xy <- xyFromCell(r, 1:ncell(r)); lat <- xy[, 2]
keep <- abs(lat) <= 30; S <- S[keep, ]; lat <- lat[keep]; w <- cos(lat * pi / 180)
cat(sprintf("Tropical cells (|lat|<=30): %d, months %d (%d-%d)\n", nrow(S), n, min(yr), max(yr)))
PER <- list(P1 = c(1950, 2025), P2 = c(1980, 2025), P3 = c(1950, 2019)); BASEN <- c(P1 = 30, P2 = 20, P3 = 30)
NB <- 1000; BL <- 5

wq <- function(x, w, p) { o <- order(x); cw <- cumsum(w[o]) / sum(w[o]); x[o][vapply(p, function(pp) which(cw >= pp)[1], 1L)] }
slope_vec <- function(Y, tt) { tc <- tt - mean(tt); as.numeric((Y - rowMeans(Y)) %*% tc) / sum(tc^2) }  # Y: cells x years
slope1 <- function(y, tt) { tc <- tt - mean(tt); sum((y - mean(y)) * tc) / sum(tc^2) }
block_idx <- function(ny, bl) { nb <- ceiling(ny / bl); st <- sample.int(ny - bl + 1, nb, replace = TRUE); (rep(st, each = bl) + rep(0:(bl - 1), nb))[1:ny] }
wmean_by <- function(x, w, g) tapply(x * w, g, sum) / tapply(w, g, sum)

out_bins <- list(); out_sum <- list(); out_q <- list(); out_thr <- list(); thrs <- c(27, 28, 29, 30)
for (pn in names(PER)) {
  y0 <- PER[[pn]][1]; y1 <- PER[[pn]][2]; yrs <- y0:y1; ny <- length(yrs)
  sel <- yr >= y0 & yr <= y1; Ssel <- S[, sel]
  ok <- rowSums(is.na(Ssel)) == 0                        # cell has all months in the period
  A <- t(vapply(yrs, function(y) rowMeans(S[, yr == y, drop = FALSE]), numeric(nrow(S))))  # years x cells (annual means)
  A <- t(A)[ok, , drop = FALSE]; wk <- w[ok]; cat(sprintf("\n== %s %d-%d: %d valid cells ==\n", pn, y0, y1, nrow(A)))
  pred <- rowMeans(A[, 1:BASEN[pn]])                      # warmth predictor: mean SST of the first years of the period
  tr <- slope_vec(A, yrs) * 10                            # K per decade
  # --- A1 bins, A2 slope ----
  bin <- cut(pred, breaks = seq(floor(min(pred)), ceiling(max(pred)), 1), include.lowest = TRUE)
  b1 <- wmean_by(tr, wk, bin); nb_cell <- table(bin)
  hi <- pred >= 27
  a2 <- function(trv, idx = seq_along(trv)) { h <- hi[idx]; ww <- wk[idx][h]; x <- pred[idx][h]; yv <- trv[h]
                                              xm <- sum(ww * x) / sum(ww); sum(ww * (x - xm) * (yv - sum(ww * yv) / sum(ww))) / sum(ww * (x - xm)^2) }
  a2_obs <- a2(tr)
  bs_a2 <- numeric(NB); bs_bin <- matrix(NA, NB, length(b1), dimnames = list(NULL, names(b1)))
  # --- B quantiles and thresholds per year (pooled monthly SST, area-weighted) ----
  Q <- t(vapply(yrs, function(y) { v <- S[, yr == y, drop = FALSE]; ww <- rep(w, ncol(v)); v <- as.numeric(v); k <- !is.na(v)
                                   c(wq(v[k], ww[k], c(.5, .9, .99)), vapply(thrs, function(th) sum(ww[k][v[k] > th]) / sum(ww[k]), 1)) }, numeric(7)))
  colnames(Q) <- c("q50", "q90", "q99", paste0("frac", thrs))
  b1_obs <- slope1(Q[, "q99"], yrs) * 10 - slope1(Q[, "q50"], yrs) * 10
  b1b_obs <- slope1(Q[, "q90"], yrs) * 10 - slope1(Q[, "q50"], yrs) * 10
  bs_b1 <- numeric(NB); bs_b1b <- numeric(NB); bs_fr <- matrix(NA, NB, length(thrs))
  for (b in 1:NB) {
    ix <- block_idx(ny, BL); tt <- yrs[ix]
    trb <- slope_vec(A[, ix], tt) * 10
    bs_a2[b] <- a2(trb); bs_bin[b, ] <- wmean_by(trb, wk, bin)
    bs_b1[b] <- (slope1(Q[ix, "q99"], tt) - slope1(Q[ix, "q50"], tt)) * 10
    bs_b1b[b] <- (slope1(Q[ix, "q90"], tt) - slope1(Q[ix, "q50"], tt)) * 10
    for (j in seq_along(thrs)) bs_fr[b, j] <- slope1(Q[ix, paste0("frac", thrs[j])], tt) * 10
  }
  ci <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
  out_bins[[pn]] <- data.frame(period = pn, bin = names(b1), ncell = as.integer(nb_cell), trend = as.numeric(b1),
                               lo = apply(bs_bin, 2, function(z) ci(z)[1]), hi = apply(bs_bin, 2, function(z) ci(z)[2]))
  out_sum[[pn]] <- data.frame(period = pn, test = c("A2 slope (K/dec per K, cells>=27C)", "B1 trend q99 - q50 (K/dec)", "B1b trend q90 - q50 (K/dec)"),
                              est = c(a2_obs, b1_obs, b1b_obs), lo = c(ci(bs_a2)[1], ci(bs_b1)[1], ci(bs_b1b)[1]), hi = c(ci(bs_a2)[2], ci(bs_b1)[2], ci(bs_b1b)[2]))
  out_q[[pn]] <- data.frame(period = pn, year = yrs, Q)
  out_thr[[pn]] <- data.frame(period = pn, threshold = thrs, mean_frac = colMeans(Q[, paste0("frac", thrs)]),
                              trend_per_dec = vapply(thrs, function(th) slope1(Q[, paste0("frac", th)], yrs) * 10, 1),
                              lo = apply(bs_fr, 2, function(z) ci(z)[1]), hi = apply(bs_fr, 2, function(z) ci(z)[2]))
  cat("Trend per predictor bin (K/decade):\n"); print(out_bins[[pn]][, -1], digits = 3, row.names = FALSE)
  cat("\nScored tests:\n"); print(out_sum[[pn]][, -1], digits = 3, row.names = FALSE)
  cat(sprintf("Tropical-ocean mean trend: %.3f K/decade; q50 %.3f, q90 %.3f, q99 %.3f K/decade\n", sum(wk * tr) / sum(wk),
              slope1(Q[, "q50"], yrs) * 10, slope1(Q[, "q90"], yrs) * 10, slope1(Q[, "q99"], yrs) * 10))
  cat("\nArea fraction above threshold (mean, trend per decade, 95% CI):\n"); print(out_thr[[pn]][, -1], digits = 3, row.names = FALSE)
}
for (nm in c("bins", "sum", "q", "thr")) write.csv(do.call(rbind, get(paste0("out_", nm))), sprintf("out/step1_%s.csv", nm), row.names = FALSE)

# --- figures ----
bb <- do.call(rbind, out_bins); bb$mid <- as.numeric(sub("^[\\[(]([-0-9.]+),.*", "\\1", bb$bin)) + 0.5
png("out/step1_fig_trend_vs_warmth.png", 1500, 700, res = 150); par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
for (pn in c("P1", "P2")) { d <- bb[bb$period == pn & bb$ncell >= 5, ]
  plot(d$mid, d$trend, type = "b", pch = 19, ylim = range(c(d$lo, d$hi)), xlab = "mean SST at start of period (C)", ylab = "SST trend (K/decade)",
       main = sprintf("%s %d-%d", pn, PER[[pn]][1], PER[[pn]][2])); arrows(d$mid, d$lo, d$mid, d$hi, 0.03, 90, 3, col = "grey40"); abline(v = 27, lty = 3) }
dev.off()
qq <- do.call(rbind, out_q); q1 <- qq[qq$period == "P1", ]
png("out/step1_fig_quantiles.png", 1500, 700, res = 150); par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
matplot(q1$year, q1[, c("q50", "q90", "q99")], type = "l", lty = 1, col = c("black", "orange", "red"), xlab = "year", ylab = "SST quantile (C)", main = "Tropical SST quantiles (P1)")
legend("topleft", c("median", "90 %", "99 %"), col = c("black", "orange", "red"), lty = 1, bty = "n")
matplot(q1$year, q1[, paste0("frac", thrs)], type = "l", lty = 1, col = 1:4, xlab = "year", ylab = "area fraction", main = "Area above threshold (P1)")
legend("topleft", paste0(">", thrs, " C"), col = 1:4, lty = 1, bty = "n"); dev.off()
cat("\nDone.\n")
