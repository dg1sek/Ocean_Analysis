# Where does "-20 W/m2 per unit stratospheric AOD" come from, and how much does it matter?  (uses data/ in this folder or the original working tree)
fe <- function(...) { p <- c(...); p[file.exists(p)][1] }
E <- fe("data/ar6_volcanic_erf.csv", "../research/D1/data/ar6_volcanic_erf.csv"); A <- fe("data/ar6_volcanic_sAOD_monthly_1850-2019.csv", "../research/D1/data/ar6_volcanic_sAOD_monthly_-50001-201912.csv")
e <- read.csv(E); a <- read.csv(A); names(a) <- c("t", "aod"); a$year <- floor(a$t); am <- aggregate(aod ~ year, a, mean); m <- merge(e, am, by = "year"); m <- m[m$year >= 1850 & m$year <= 2019, ]
f <- lm(volcanic_erf ~ aod, m)
cat(sprintf("AR6 volcanic ERF (annual) regressed on annual-mean AR6 sAOD, %d years: slope %.3f W/m2 per unit sAOD, intercept %.3f, R2 = %.6f\n", nrow(m), coef(f)[2], coef(f)[1], summary(f)$r.squared))
cat("-> the AR6 ERF series is exactly -20 x (sAOD - baseline): the factor is an input of the construction, not an independent measurement.\n")
# effect on our fitted gains: gain (K per unit sAOD, from out/volcano_onebox_fit.csv) divided by the assumed factor
B <- read.csv(fe("out/volcano_onebox_fit.csv")); fac <- c(-13, -17, -20, -25)
t <- sapply(fac, function(k) B$gain_K_per_sAOD / k); colnames(t) <- paste0("factor", fac); rownames(t) <- B$series
cat("\nSST response per W/m2 (K per W/m2) for different assumed conversion factors (W/m2 per unit sAOD):\n"); print(round(t, 3))
write.csv(data.frame(series = B$series, gain_K_per_sAOD = B$gain_K_per_sAOD, round(t, 3)), "out/volcano_gain_vs_factor.csv", row.names = FALSE)

# ---- effective heat capacity implied by the fitted memory time and gain (one-box: tau = C x gain, gain in K per W/m2) -------------------
# C = tau / gain (J m-2 K-1); equivalent sea-water layer depth h = C / (rho cp). Compare with the 50-100 m mixed layer quoted by Alf (2026-11-23).
rho_cp <- 4.1e6; sec_month <- 3.156e7 / 12
cat("\nEffective heat capacity implied by the one-box fit, as an equivalent layer of sea water (m):\n")
O <- read.csv(fe("out/volcano_onebox_fit.csv")); res <- list()
for (k in c(-13, -17, -20, -25)) { g <- O$gain_K_per_sAOD / k; h <- (O$tau_months * sec_month / g) / rho_cp; res[[as.character(k)]] <- data.frame(series = O$series, factor = k, tau_months = O$tau_months, gain_K_per_Wm2 = round(g, 3), equiv_depth_m = round(h)) }
H <- do.call(rbind, res); print(H[H$factor == -20, ], row.names = FALSE); cat("global and warm pool for all factors:\n"); print(H[H$series %in% c("GLOBAL", "WPWP"), ], row.names = FALSE)
write.csv(H, "out/volcano_effective_depth.csv", row.names = FALSE)
