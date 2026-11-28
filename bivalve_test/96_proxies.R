# 96_proxies.R -- Barents Sea shells (Mette 2021) and tree-ring AMO reconstruction (Gray 2004), see PREREG_bivalve.md part 2. Run: Rscript 96_proxies.R
source("R/alf_tools.R"); set.seed(2026)
m <- read.delim(text = { l <- readLines("data/mette2021.txt", encoding = "latin1"); l[grep("^age_calCE", l)[1]:length(l)] }, na.strings = c("-999.9", "-999"), stringsAsFactors = FALSE)
g <- read.delim("data/gray2004.txt", comment.char = "#", na.strings = "NA", stringsAsFactors = FALSE); g <- g[grep("^age$", names(g)) , drop = FALSE]  # placeholder
l <- readLines("data/gray2004.txt", encoding = "latin1"); g <- read.delim(text = l[grep("^age\\tssta", l)[1]:length(l)], na.strings = "NA", stringsAsFactors = FALSE)
P <- list(BARENTS_growth = data.frame(t = m$age_calCE, y = m$MShellCrn), BARENTS_d18O = data.frame(t = m$age_calCE, y = m$d18Oa.isla), GRAY_ssta = data.frame(t = g$age, y = g$ssta))
P <- lapply(P, function(d) { d <- d[!is.na(d$y), ]; d[order(d$t), ] }); for (n in names(P)) cat(sprintf("%-15s %d-%d (n = %d)\n", n, min(P[[n]]$t), max(P[[n]]$t), nrow(P[[n]])))
per <- 40:120
scan <- function(y, t, nm) { y <- (y - mean(y))/sd(y); s <- sin_scan(y, t, per, 1); j <- which.max(s$dAIC); nl <- null_scan(y, t, per, 1, nsim = 300)
  data.frame(series = nm, from = min(t), to = max(t), n = length(y), best_P = s$period[j], amp_sd = round(s$amp[j], 2), dAIC = round(s$dAIC[j], 1), null95 = round(quantile(nl$best, .95), 1), p = round(mean(nl$best >= s$dAIC[j]), 3), band_65_80_dAIC = round(max(s$dAIC[s$period >= 65 & s$period <= 80]), 1), phi = round(s$phi[1], 2)) }
sub <- function(d, a, b) d[d$t >= a & d$t <= b, ]
C1 <- list()
for (n in c("BARENTS_growth", "BARENTS_d18O")) { d <- P[[n]]; a <- min(d$t); b <- max(d$t); mid <- round((a + b)/2)
  C1[[length(C1)+1]] <- with(d, scan(y, t, paste(n, "full"))); h1 <- sub(d, a, mid); h2 <- sub(d, mid + 1, b)
  C1[[length(C1)+1]] <- scan(h1$y, h1$t, paste(n, "first half")); C1[[length(C1)+1]] <- scan(h2$y, h2$t, paste(n, "second half")) }
d <- P$GRAY_ssta; C1[[length(C1)+1]] <- scan(d$y, d$t, "GRAY_ssta full 1567-1990"); h <- sub(d, 1567, 1855); C1[[length(C1)+1]] <- scan(h$y, h$t, "GRAY_ssta pre-instrumental 1567-1855"); h <- sub(d, 1856, 1990); C1[[length(C1)+1]] <- scan(h$y, h$t, "GRAY_ssta instrumental 1856-1990")
C1 <- do.call(rbind, C1); cat("\n== C1: best period 40-120 yr ==\n"); print(C1, row.names = FALSE); write.csv(C1, "out/proxies_C1_line.csv", row.names = FALSE)
# C2 phase at 75 yr
Pfix <- 75; w <- 2*pi/Pfix; wrap <- function(d) ((d + Pfix/2) %% Pfix) - Pfix/2
pf <- function(y, t) { x <- (t - mean(t))/100; X <- cbind(1, x, x^2, sin(w*t), cos(w*t)); f <- lm.fit(X, y); list(cf = f$coefficients, res = f$residuals, X = X) }
tmax <- function(cf, tend) { ph <- atan2(cf[4], cf[5]); k <- floor((tend*w - ph)/(2*pi)); (ph + 2*pi*k)/w }
rg <- read.csv("data/regional_monthly_ersst.csv"); rg <- rg[rg$year <= 2025, ]; mt <- tapply(rg$MED, rg$year, mean)
sets <- list(BARENTS_growth = P$BARENTS_growth, BARENTS_d18O = P$BARENTS_d18O, GRAY_ssta_pre = sub(P$GRAY_ssta, 1567, 1855), GRAY_ssta_full = P$GRAY_ssta, MED_SST_ERSST_1850_2025 = data.frame(t = as.integer(names(mt)), y = as.numeric(mt)))
C2 <- do.call(rbind, lapply(names(sets), function(n) { d <- sets[[n]]; y <- (d$y - mean(d$y))/sd(d$y); t <- d$t; f <- pf(y, t); tend <- 2022; tm <- tmax(f$cf, tend); N <- length(y); BL <- 5
  bs <- replicate(200, { st <- sample.int(N - BL + 1, ceiling(N/BL), TRUE); rb <- unlist(lapply(st, function(s) f$res[s:(s + BL - 1)]))[1:N]; b <- pf(as.numeric(f$X %*% f$cf) + rb, t); tmax(b$cf, tend) })
  dd <- wrap(bs - tm); data.frame(series = n, from = min(t), to = max(t), last_max_year_before_2022 = round(tm, 1), lo = round(tm + quantile(dd, .025), 1), hi = round(tm + quantile(dd, .975), 1), offset_to_MED_max_2022_yr = round(wrap(tm - 2022.2), 1), amp_sd = round(sqrt(f$cf[4]^2 + f$cf[5]^2)/sd(y), 2)) }))
cat("\n== C2: year of the last maximum at fixed 75 yr (<= 2022), offset to the ERSST Mediterranean maximum (2022.2) modulo 75 yr ==\n"); print(C2, row.names = FALSE); write.csv(C2, "out/proxies_C2_phase.csv", row.names = FALSE)
png("out/proxies_fig.png", 2400, 1800, res = 200); par(mfrow = c(3, 1), mar = c(4, 4, 2.5, 1))
for (n in c("BARENTS_growth", "BARENTS_d18O", "GRAY_ssta")) { d <- P[[n]]; plot(d$t, d$y, type = "l", col = "grey60", xlab = "Jahr", ylab = n, main = n); k <- 31; lines(d$t, stats::filter(d$y, rep(1/k, k), sides = 2), col = "#d62728", lwd = 2) }
dev.off()
