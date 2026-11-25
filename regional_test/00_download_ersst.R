# Download ERSSTv6 monthly SST (NOAA PSL, about 86 MB) into data/. Not stored in the repository.
dir.create("data", showWarnings = FALSE)
if (!file.exists("data/ersst.v6.sst.mnmean.nc"))
  download.file("https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v6/sst.mnmean.nc", "data/ersst.v6.sst.mnmean.nc", mode = "wb")
