# 51_sinus_phase.R -- extra check (labelled as such, after seeing 50_slope.R): where is the fitted sinusoid in its cycle now?
source("R/alf_tools.R")
rg <- read.csv("data/regional_monthly_ersst.csv"); gl <- read.csv("data/global_series.csv")
ann <- function(df, col) tapply(df[[col]], df$year, mean)
S <- list(GLOBAL = ann(gl, "ersst"), MED = ann(rg, "MED"), GOM = ann(rg, "GOM"))
g <- expand.grid(nm = names(S), from = c(1850, 1900), stringsAsFactors = FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) { nm <- g$nm[i]; from <- g$from[i]
  y0 <- S[[nm]]; t0 <- as.numeric(names(y0)); k <- t0 >= from; y <- as.numeric(y0[k]); t <- t0[k]
  sc <- sin_scan(y, t, 60:100, 2); j <- which.max(sc$dAIC); P <- sc$period[j]; ph <- sc$phase[j]; A <- sc$amp[j]
  tt <- seq(1850, 2100, .25); s <- A * sin(2*pi*tt/P + ph); ds <- diff(s)/diff(tt)
  mxs <- tt[which(diff(sign(ds)) < 0)]
  data.frame(series = nm, from = from, P = P, amp = round(A, 2), sinus_slope_2025_K_dec = round(10*ds[which.min(abs(tt-2025))], 3),
    last_maxima = paste(round(mxs[mxs < 2060 & mxs > 1900], 0), collapse = ","), next_max_after_2000 = round(mxs[mxs > 2000][1], 0))
}))
print(res); write.csv(res, "out/sinus_phase.csv", row.names = FALSE)
