# ---------------------------------------------------------------------------
# helpers.R  -- shared settings and functions for the follow-up tests
# Base R only. Sourced by every numbered script.
# ---------------------------------------------------------------------------

BASE_PERIOD <- c(1991, 2020)      # common reference period for all anomalies
SERIES      <- Sys.getenv("SERIES", "ersst")  # column of output/global_series.csv used by
                                  # scripts 02-05: "ersst", "ersst_noice",
                                  # "ersst_median", "noaa", "hadsst" or "alf"
ALF_RDS     <- "../data/alf_series_rebuilt.rds"   # $data of NOAA.OCEAN.ANOMALIES.rds (Alf's scripts, rebuilt 22 Nov)
ALF_COLUMN  <- "anoma.mean"
SEED        <- 2026

dir.create("data",   showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

dec_year <- function(year, month) year + (month - 0.5) / 12

# Months of solar-cycle minima, cycles 9 to 25, from the SILSO table "Table of
# minima, maxima and cycle durations" (13-month smoothed sunspot number V2.0,
# https://www.sidc.be/SILSO/cyclesmm; read on 15 Nov 2026). Row i is the start
# of cycle 8 + i. The spline test (02) and the composite (03) take their knots
# and epochs from here, and 03 takes the cycle number for the polarity sign.
# Replace the table with the dates used in the paper if they differ.
SOLAR_MIN_TABLE <- data.frame(
  cycle = 9:25,
  year  = c(1843, 1855, 1867, 1878, 1890, 1902, 1913, 1923, 1933,
            1944, 1954, 1964, 1976, 1986, 1996, 2008, 2019),
  month = c(   7,   12,    3,   12,    3,    1,    7,    8,    9,
                2,    4,   10,    3,    9,    8,   12,   12))
SOLAR_MINIMA      <- dec_year(SOLAR_MIN_TABLE$year, SOLAR_MIN_TABLE$month)
SOLAR_CYCLE_FIRST <- SOLAR_MIN_TABLE$cycle[1]

# Solar-cycle number of decimal-year times (NA before the first minimum)
solar_cycle <- function(t, minima = SOLAR_MINIMA, first = SOLAR_CYCLE_FIRST) {
  i <- findInterval(t, minima)
  ifelse(i == 0, NA_integer_, first + i - 1L)
}
# Hale polarity sign: +1 in odd cycles, -1 in even cycles. The convention is
# arbitrary; Hale's polarity law only says that the sign alternates from one
# cycle to the next, so that the magnetic pattern repeats every second cycle.
hale_sign <- function(t) ifelse(solar_cycle(t) %% 2 == 1, 1, -1)

# Subtract the calendar-month means of the base period (year, month, y)
rebaseline <- function(year, month, y, base = BASE_PERIOD) {
  inb  <- year >= base[1] & year <= base[2] & !is.na(y)
  clim <- tapply(y[inb], month[inb], mean)
  y - unname(clim[as.character(month)])
}

# Monthly series chosen by SERIES: data.frame(year, month, t, y), no gaps kept
load_monthly <- function(col = SERIES) {
  f <- "output/global_series.csv"
  if (!file.exists(f)) stop("Run 01_global_mean.R first (", f, " missing).")
  d <- read.csv(f)
  if (!col %in% names(d)) stop("Column '", col, "' not in ", f)
  d <- data.frame(year = d$year, month = d$month, t = d$t, y = d[[col]])
  d[!is.na(d$y), ]
}

# Annual means, complete years only: data.frame(t = year, y)
annual_mean <- function(d) {
  n <- tapply(d$y, d$year, length)
  m <- tapply(d$y, d$year, mean)
  ok <- n == 12
  data.frame(t = as.numeric(names(m))[ok], y = as.numeric(m)[ok])
}

# AR(1) noise with innovation s.d. 'sd_innov'
ar1_sim <- function(n, phi, sd_innov) {
  e <- rnorm(n + 200, 0, sd_innov)
  x <- as.numeric(stats::filter(e, phi, method = "recursive"))
  x[-(1:200)]
}

# Random knot set with the same spacing statistics as the real solar minima:
# intervals are resampled from the observed ones, the start is random.
random_knots <- function(tmin, tmax, minima = SOLAR_MINIMA) {
  iv <- diff(minima)
  k  <- tmin - runif(1, 0, mean(iv))
  while (k[length(k)] < tmax) k <- c(k, k[length(k)] + sample(iv, 1))
  k
}

# Period (in units of dt) of the strongest spectral peak of a vector
dominant_period <- function(v, dt = 1, pad = 8192) {
  v <- v - mean(v)
  z <- c(v, rep(0, max(0, pad - length(v))))
  n <- length(z)
  P <- Mod(fft(z))^2
  k <- 2:(n %/% 2)
  f <- (k - 1) / (n * dt)
  1 / f[which.max(P[k])]
}

# --- minimal SSA (no package needed) ---------------------------------------
# Trajectory matrix X is K x L, row i = x[i], ..., x[i+L-1].
ssa_decomp <- function(x, L) {
  X  <- embed(x, L)[, L:1, drop = FALSE]
  sv <- svd(X)
  list(X = X, u = sv$u, d = sv$d, v = sv$v, L = L, K = nrow(X), N = length(x))
}
# Reconstruction of the components 'idx' by diagonal averaging
ssa_recon <- function(s, idx) {
  Xr <- s$u[, idx, drop = FALSE] %*% (s$d[idx] * t(s$v[, idx, drop = FALSE]))
  as.numeric(tapply(Xr, row(Xr) + col(Xr) - 1, mean))
}
