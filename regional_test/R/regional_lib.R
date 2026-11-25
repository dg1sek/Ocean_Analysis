# ---------------------------------------------------------------------------
# regional_lib.R -- shared helpers for the regional analyses (base R + terra/ncdf4)
# Paths are relative to the folder regional/ (run scripts from there).
# ---------------------------------------------------------------------------
suppressMessages(library(terra))
first_existing <- function(...) { p <- c(...); p[file.exists(p)][1] }
source(first_existing("R/read_nc.R", "../followup/R/read_nc.R"))   # rast_nc(): ncdf4 fallback when GDAL has no netCDF driver
source(first_existing("R/helpers.R", "../followup/R/helpers.R"))   # BASE_PERIOD, dec_year, SOLAR_MIN_TABLE, random_knots, ...
ERSST_FILE <- first_existing("data/ersst.v6.sst.mnmean.nc", "../followup/data/ersst.v6.sst.mnmean.nc")
SILSO_FILE <- first_existing("data/silso_monthly.csv", "../followup/data/silso_monthly.csv")
if (is.na(ERSST_FILE)) stop("ERSSTv6 file missing: run Rscript 00_download_ersst.R")
dir.create("out", showWarnings = FALSE); dir.create("data", showWarnings = FALSE)

# Candidate regions: lat1, lat2, lon1, lon2 (degrees east, 0-360; lon1 > lon2 means the box wraps over 0 E)
REGIONS <- data.frame(
  id   = c("NINO34","WPWP","SCS","SWPAC","TIND","MED","GOM","NWPAC","NPAC","NATL_SUBPOLAR","NATL_SUBTROP","REDSEA"),
  name = c("Nino3.4 (central equatorial Pacific)","Western Pacific warm pool","South China Sea",
           "South-west Pacific / SPCZ (Fiji, New Caledonia, Rarotonga)","Tropical Indian Ocean","Mediterranean Sea",
           "Gulf of Mexico","North-west Pacific (Kuroshio, Ryukyu)","North Pacific mid-latitudes",
           "Subpolar North Atlantic","Subtropical North Atlantic (Sargasso)","Red Sea"),
  lat1 = c(-5, -10,   5, -25, -20, 30, 18, 20, 30, 50, 25, 13),
  lat2 = c( 5,  10,  22, -10,  10, 46, 30, 35, 45, 65, 40, 28),
  lon1 = c(190, 130, 105, 160,  50, 354, 262, 120, 150, 300, 280, 32),
  lon2 = c(240, 170, 120, 205, 100,  36, 279, 145, 210, 350, 320, 44),
  stringsAsFactors = FALSE)

in_box <- function(lat, lon, b) {
  ok_lat <- lat >= b$lat1 & lat <= b$lat2
  ok_lon <- if (b$lon1 <= b$lon2) lon >= b$lon1 & lon <= b$lon2 else lon >= b$lon1 | lon <= b$lon2
  ok_lat & ok_lon
}

# ERSSTv6 monthly anomaly cube (cells x months), anomalies against the 1991-2020 calendar-month climatology, per cell
load_ersst_cells <- function(f = ERSST_FILE) {
  r <- tryCatch(rast(f, subds = "sst"), error = function(e) tryCatch(rast(f), error = function(e2) rast_nc(f, "sst")))
  n <- nlyr(r); tm <- seq(as.Date("1850-01-01"), by = "month", length.out = n)
  year <- as.integer(format(tm, "%Y")); month <- as.integer(format(tm, "%m"))
  sst <- values(r); sst[sst < -5 | sst > 45] <- NA
  xy <- xyFromCell(r, 1:ncell(r)); lon <- (xy[, 1] + 360) %% 360; lat <- xy[, 2]
  inb <- year >= BASE_PERIOD[1] & year <= BASE_PERIOD[2]
  anom <- sst
  for (m in 1:12) {
    clim <- rowMeans(sst[, inb & month == m, drop = FALSE], na.rm = TRUE)
    anom[, month == m] <- sst[, month == m] - clim
  }
  list(anom = anom, sst = sst, lat = lat, lon = lon, year = year, month = month, t = dec_year(year, month))
}

regional_series <- function(G, b) {
  idx <- which(in_box(G$lat, G$lon, b)); w <- cos(G$lat[idx] * pi / 180)
  A <- G$anom[idx, , drop = FALSE]
  v <- colSums(A * w, na.rm = TRUE) / colSums((!is.na(A)) * w)
  list(v = v, ncell = length(idx), nvalid = colSums(!is.na(A)))
}

annual <- function(year, month, v) {
  ok <- tapply(!is.na(v), year, sum) == 12
  a <- tapply(v, year, mean); yrs <- as.integer(names(a))
  data.frame(year = yrs[ok], v = as.numeric(a[ok]))
}

# Daily-mean top-of-atmosphere insolation (W m-2) by day of year (standard formula, solar constant 1361)
toa_insolation <- function(lat, S0 = 1361) {
  doy <- 1:365; phi <- lat * pi / 180
  decl <- 23.44 * pi / 180 * sin(2 * pi * (doy - 81) / 365)
  e <- 1 + 0.033 * cos(2 * pi * doy / 365)
  h0 <- acos(pmax(-1, pmin(1, -tan(phi) * tan(decl))))
  S0 / pi * e * (h0 * sin(phi) * sin(decl) + cos(phi) * cos(decl) * sin(h0))
}

# AR(1) fit and surrogate
ar1_fit <- function(x) { x <- x - mean(x); phi <- sum(x[-1] * x[-length(x)]) / sum(x^2); c(phi = phi, sd = sd(x) * sqrt(1 - phi^2)) }
