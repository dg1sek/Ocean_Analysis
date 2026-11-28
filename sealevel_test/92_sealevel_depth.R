# 92_sealevel_depth.R -- apparent heat-reservoir depth from sea level vs SST (PREREG_sealevel.md part 2). Run: Rscript 92_sealevel_depth.R
set.seed(2026)
ps <- read.csv("data/psmsl_rlr_annual.csv"); ps <- ps[ps$year <= 2025, ]; Sst <- split(ps, ps$station)
anom <- function(n) { d <- Sst[[n]]; data.frame(year = d$year, h = d$rlr_mm - mean(d$rlr_mm[d$year >= 1890 & d$year <= 2020])) }
med <- do.call(rbind, lapply(c("TRIESTE", "GENOVA", "MARSEILLE"), anom)); cm <- tapply(med$h, med$year, function(z) if (length(z) >= 2) mean(z) else NA); cm <- cm[!is.na(cm)]
H <- list(MED_COMPOSITE = data.frame(year = as.integer(names(cm)), h = as.numeric(cm)), TRIESTE = anom("TRIESTE"), GENOVA = anom("GENOVA"), MARSEILLE = anom("MARSEILLE"))
rg <- read.csv("data/regional_monthly_ersst.csv"); rg <- rg[rg$year <= 2025, ]; mt <- tapply(rg$MED, rg$year, mean); SS <- data.frame(year = as.integer(names(mt)), sst = as.numeric(mt))
alpha <- 2.3e-4; dep <- function(sl, a = alpha) sl/(a*1000)
out <- list()
for (nm in names(H)) { d <- merge(H[[nm]], SS, by = "year"); d <- d[order(d$year), ]
  for (k in c(1, 5, 10)) { i <- which(d$year %in% (d$year + k)); pair <- do.call(rbind, lapply(d$year, function(y) { j <- match(y - k, d$year); if (is.na(j)) NULL else data.frame(year = y, dh = d$h[d$year == y] - d$h[j], ds = d$sst[d$year == y] - d$sst[j]) }))
    f <- lm(dh ~ ds, pair); sl <- coef(f)[["ds"]]; N <- nrow(pair); BL <- max(5, k); pair <- pair[order(pair$year), ]
    bs <- replicate(500, { st <- sample.int(N - BL + 1, ceiling(N/BL), TRUE); idx <- unlist(lapply(st, function(s) s:(s + BL - 1)))[1:N]; coef(lm(dh ~ ds, pair[idx, ]))[["ds"]] })
    ci <- quantile(bs, c(.025, .975)); out[[length(out)+1]] <- data.frame(series = nm, k_yr = k, n = N, slope_mm_per_K = round(sl, 1), lo = round(ci[1], 1), hi = round(ci[2], 1), R2 = round(summary(f)$r.squared, 3), H_m = round(dep(sl)), H_lo_alpha_plus20 = round(dep(sl, alpha*1.2)), H_hi_alpha_minus20 = round(dep(sl, alpha*0.8)), excludes_zero = ci[1] > 0 | ci[2] < 0)
    print(out[[length(out)]], row.names = FALSE) } }
R <- do.call(rbind, out); write.csv(R, "out/sealevel_depth.csv", row.names = FALSE)
d <- merge(H$MED_COMPOSITE, SS, by = "year"); d <- d[order(d$year), ]; dh <- diff(d$h)[diff(d$year) == 1]; ds <- diff(d$sst)[diff(d$year) == 1]
cat("\nlag of the strongest correlation of 1-yr differences (composite; positive = SST leads sea level):\n"); yrs <- d$year[-1][diff(d$year) == 1]
dd <- data.frame(year = yrs, dh = dh, ds = ds); cc <- sapply(-3:3, function(L) { a <- dd; b <- dd; b$year <- b$year + L; m <- merge(a, b, by = "year"); cor(m$dh.x, m$ds.y) }); print(round(setNames(cc, -3:3), 2))
