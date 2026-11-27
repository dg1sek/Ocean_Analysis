# Build monthly ocean-only latitude-band means from ERSSTv6 (needs the ERSST file; see ../regional_test/00_download_ersst.R). Writes data/latband_monthly_ersst.csv
suppressMessages(library(terra)); fe <- function(...) { p <- c(...); p[file.exists(p)][1] }
source("R/read_nc.R")
E <- fe("data/ersst.v6.sst.mnmean.nc", "../regional/data/ersst.v6.sst.mnmean.nc", "../followup/data/ersst.v6.sst.mnmean.nc", "../thermostat/data/ersst.v6.sst.mnmean.nc", "../regional_test/data/ersst.v6.sst.mnmean.nc"); stopifnot(!is.na(E))
r <- tryCatch(rast(E, subds = "sst"), error = function(e) tryCatch(rast(E), error = function(e2) rast_nc(E, "sst")))
n <- nlyr(r); tm <- seq(as.Date("1850-01-01"), by = "month", length.out = n); yr <- as.integer(format(tm, "%Y")); mo <- as.integer(format(tm, "%m"))
S <- values(r); S[S < -5 | S > 45] <- NA; lat <- xyFromCell(r, 1:ncell(r))[, 2]
inb <- yr >= 1991 & yr <= 2020; A <- S; for (m in 1:12) A[, mo == m] <- S[, mo == m] - rowMeans(S[, inb & mo == m], na.rm = TRUE)
bands <- list("60-40S" = c(-60, -40), "40-20S" = c(-40, -20), "20S-0" = c(-20, 0), "0-20N" = c(0, 20), "20-40N" = c(20, 40), "40-60N" = c(40, 60))
out <- data.frame(year = yr, month = mo)
for (b in names(bands)) { k <- which(lat >= bands[[b]][1] & lat < bands[[b]][2]); w <- cos(lat[k] * pi / 180); X <- A[k, , drop = FALSE]
  out[[b]] <- colSums(X * w, na.rm = TRUE) / colSums((!is.na(X)) * w); cat(b, ": cells", length(k), ", ocean cells with data in Jan 2000:", sum(!is.na(X[, which(yr == 2000 & mo == 1)])), "\n") }
dir.create("data", showWarnings = FALSE); write.csv(out, "data/latband_monthly_ersst.csv", row.names = FALSE); cat("written\n")
