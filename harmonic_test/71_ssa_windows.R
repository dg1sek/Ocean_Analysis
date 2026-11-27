# 71_ssa_windows.R -- SSA window length diagnostic (PREREG_oos_ssa.md, part B). Run: Rscript 71_ssa_windows.R
suppressMessages(library(Rssa))
rg <- read.csv("data/regional_monthly_ersst.csv"); gl <- read.csv("data/global_series.csv")
ann <- function(df, col) tapply(df[[col]], df$year, mean)
S <- list(GLOBAL = ann(gl, "ersst"), MED = ann(rg, "MED")); N <- 176; Ls <- c(N/2, floor(N/3), N/4, 40, 60, 80)
per_of <- function(v) { f <- Mod(fft(v - mean(v)))[2:(length(v) %/% 2)]; k <- which.max(f); length(v) / k }
out <- list()
for (nm in names(S)) for (L in Ls) { y <- as.numeric(S[[nm]]); s <- ssa(y, L = L); ne <- 8
  W <- wcor(s, groups = as.list(1:ne)); sh <- round(100 * s$sigma[1:ne]^2 / sum(s$sigma^2), 1)
  U <- s$U[, 1:ne, drop = FALSE]; pe <- round(apply(U, 2, per_of), 1)
  W13 <- W[1:3, 1:3]; m13 <- max(abs(W13[upper.tri(W13)])); wmax <- apply(abs(W[1:ne, 1:ne] - diag(ne)), 1, max)
  has <- any(pe >= 60 & pe <= 100 & seq_len(ne) > 1)
  out[[length(out)+1]] <- data.frame(series = nm, L = L, max_wcor_1to3 = round(m13, 2), pair_60_100 = has, share_1to3 = paste(sh[1:3], collapse = "/"), periods_1to8 = paste(pe, collapse = ","))
  print(out[[length(out)]], row.names = FALSE) }
write.csv(do.call(rbind, out), "out/ssa_windows.csv", row.names = FALSE)
