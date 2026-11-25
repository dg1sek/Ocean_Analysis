# ---------------------------------------------------------------------------
# 10_region_screening.R -- screen candidate regions BEFORE looking at any solar or 20-yr signal.
# Criteria are chosen a priori so that the selection cannot be tuned to the result:
#   proxy   independent long records inside the box (coral SST etc., already downloaded; see data/proxy_sites.csv)
#   early   ERSSTv6 early-period reliability: ratio of annual detrended sd 1850-99 / 1950-99
#           (values far below 1 mean the early decades are mostly statistical infill)
#   enso    absolute correlation of the region with the Nino3.4 box (annual, high-passed, 1870-2025); low = cleaner
#   noise   annual high-passed standard deviation of the regional mean (low = better signal to noise)
#   alf     named by Alf in his mails (Mediterranean, Gulf of Mexico, Indian Ocean, mid-latitude Pacific, equatorial ENSO area)
# Solar geometry (annual-mean and seasonal range of TOA insolation at the box centre) is reported, not scored.
# Output: out/regions_screening.csv, out/regional_monthly_ersst.csv
# ---------------------------------------------------------------------------
source("R/regional_lib.R")
G <- load_ersst_cells()
cat("cells with a value in every month:", sum(rowSums(!is.na(G$anom)) == ncol(G$anom)), "\n")

ser <- list(); tab <- list()
for (i in seq_len(nrow(REGIONS))) {
  b <- REGIONS[i, ]; s <- regional_series(G, b); ser[[b$id]] <- s$v
  a <- annual(G$year, G$month, s$v)
  hp <- function(y, t) resid(lm(y ~ poly(t, 2)))                 # quadratic detrend = "high-passed"
  a$hp <- hp(a$v, a$year)
  sd_all <- sd(a$hp)
  sd_early <- sd(a$hp[a$year <= 1899]); sd_late <- sd(a$hp[a$year >= 1950 & a$year <= 1999])
  ## Nino3.4 from the same grid (box 5S-5N, 170W-120W)
  n34 <- regional_series(G, REGIONS[REGIONS$id == "NINO34", ])$v
  an34 <- annual(G$year, G$month, n34); an34$hp <- hp(an34$v, an34$year)
  m <- merge(a, an34, by = "year", suffixes = c("", ".n34")); m <- m[m$year >= 1870, ]
  r_enso <- if (b$id == "NINO34") NA else cor(m$hp, m$hp.n34)
  lat0 <- mean(c(b$lat1, b$lat2)); S <- toa_insolation(lat0)
  tab[[i]] <- data.frame(id = b$id, name = b$name, lat_centre = lat0, ncell = s$ncell,
    insol_mean_Wm2 = round(mean(S)), insol_seasonal_range_Wm2 = round(diff(range(S))),
    sd_annual_hp_K = round(sd_all, 3), early_late_sd_ratio = round(sd_early / sd_late, 2),
    abs_r_nino34 = round(abs(r_enso), 2))
}
tab <- do.call(rbind, tab)

# Independent proxy records (all downloaded in research/A2/data and W2; coordinates from the NOAA file headers)
proxy <- read.csv("data/proxy_sites.csv", stringsAsFactors = FALSE)
tab$proxy_records <- sapply(tab$id, function(id) sum(proxy$region == id & proxy$years >= 270))
tab$proxy_names   <- sapply(tab$id, function(id) paste(proxy$site[proxy$region == id], collapse = "; "))

# Scores (0-3 proxy, 0-2 early, 0-2 enso, 0-2 noise, 0-1 alf); proxy scores 1 if only catalogued, not downloaded
cat_only <- c(TIND = 1, NATL_SUBTROP = 1, REDSEA = 1)             # Mentawai/Ifaty; Bahamas; Umh Huj/Aqaba exist at NOAA NCEI
tab$s_proxy <- pmin(3, ifelse(tab$proxy_records >= 2, 3, ifelse(tab$proxy_records == 1, 2, cat_only[tab$id])))
tab$s_proxy[is.na(tab$s_proxy)] <- 0
tab$s_early <- ifelse(tab$early_late_sd_ratio >= 0.8, 2, ifelse(tab$early_late_sd_ratio >= 0.6, 1, 0))
tab$s_enso  <- ifelse(is.na(tab$abs_r_nino34), 0, ifelse(tab$abs_r_nino34 < 0.3, 2, ifelse(tab$abs_r_nino34 < 0.6, 1, 0)))
tab$s_noise <- ifelse(tab$sd_annual_hp_K <= 0.15, 2, ifelse(tab$sd_annual_hp_K <= 0.25, 1, 0))
tab$s_alf   <- as.integer(tab$id %in% c("MED", "GOM", "TIND", "NPAC", "NWPAC", "NINO34", "WPWP"))
tab$score   <- tab$s_proxy + tab$s_early + tab$s_enso + tab$s_noise + tab$s_alf
tab$band    <- cut(abs(tab$lat_centre), c(0, 15, 30, 90), labels = c("tropical 0-15", "subtropical 15-30", "mid/high 30-65"), right = FALSE)
tab <- tab[order(-tab$score, tab$id), ]
write.csv(tab, "out/regions_screening.csv", row.names = FALSE)
write.csv(data.frame(t = G$t, year = G$year, month = G$month, as.data.frame(ser)), "out/regional_monthly_ersst.csv", row.names = FALSE)
print(tab[, c("id", "band", "lat_centre", "ncell", "insol_mean_Wm2", "insol_seasonal_range_Wm2", "sd_annual_hp_K",
              "early_late_sd_ratio", "abs_r_nino34", "proxy_records", "score")], row.names = FALSE)
cat("\nBest in each band:\n")
for (bd in levels(tab$band)) { x <- tab[tab$band == bd, ][1, ]; cat(sprintf("  %-18s %s (score %d)\n", bd, x$id, x$score)) }
