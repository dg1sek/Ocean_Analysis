Data needed (not stored in the repository, download once into thermostat_test/data/):
 ERSSTv6        https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v6/sst.mnmean.nc  (as ersst.v6.sst.mnmean.nc; see ../regional_test/00_download_ersst.R)
 OISST v2 1deg  https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2/sst.mnmean.nc  (as oisst.v2.sst.mnmean.nc; ends Jan 2023) and .../noaa.oisst.v2/lsmask.nc (as lsmask.nc; the SST file fills land)
 HadSST 4.2     https://www.metoffice.gov.uk/hadobs/hadsst4/data/data/HadSST.4.2.0.0_actuals_median.nc
 NCEP/NCAR R1   https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss/{dswrf,uswrf,dlwrf,ulwrf,lhtfl,shtfl}.sfc.mon.mean.nc  (saved as ncep_<name>.nc)
