# 11_export_proxies.R -- annual proxy records (from the NOAA NCEI files downloaded for the coral test) as one tidy CSV.
# Sign-corrected so that up = warm. Source loader: ../research/A2/A2_load.R (writes research/A2/out/records.rds).
rec <- readRDS("../research/A2/out/records.rds")$R
keep <- c("DryTortugas_SST", "HonTre_SST", "Fiji_AB_d18O", "Amedee_SrCa", "Rarotonga_SrCa", "Ryukyu_SST", "Freund_NCT", "Wilson_N34_coral")
out <- do.call(rbind, lapply(keep, function(k) { r <- rec[[k]]; data.frame(record = k, site = r$site, variable = r$var, year = r$df$year, value = r$df$x) }))
d <- "../research/A2/data/"
x <- read.table(paste0(d, "freund2019enso.txt"), header = TRUE, sep = "\t", comment.char = "#", na.strings = c("NA", "-999"), quote = "", fill = TRUE)
m <- tapply(x$NWP, x$age_CE, mean); n <- tapply(x$NWP, x$age_CE, length)
nw <- data.frame(record = "Freund_NWP", site = "composite", variable = "Nino Warm Pool index, mean of 4 seasons",
                 year = as.integer(names(m))[n >= 4], value = as.numeric(m)[n >= 4])
out <- rbind(out, nw)
write.csv(out, "data/proxy_annual.csv", row.names = FALSE)
print(aggregate(year ~ record, out, function(y) paste(range(y), collapse = "-")))
