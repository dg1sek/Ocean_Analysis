# ---------------------------------------------------------------------------
# read_nc.R -- read a regular lon/lat netCDF variable into a terra SpatRaster
# through the ncdf4 package. Fallback for installations whose GDAL has no
# netCDF driver (e.g. MacPorts gdal without the +netcdf variant).
# Result: north-up raster, one layer per time step, EPSG:4326, extent spans
# the cell edges; values NA where the file flags missing data.
# ---------------------------------------------------------------------------
rast_nc <- function(f, var = "sst") {
  stopifnot(requireNamespace("ncdf4", quietly = TRUE), requireNamespace("terra", quietly = TRUE))
  nc <- ncdf4::nc_open(f); on.exit(ncdf4::nc_close(nc))
  lon <- as.numeric(ncdf4::ncvar_get(nc, "lon")); lat <- as.numeric(ncdf4::ncvar_get(nc, "lat"))
  a <- ncdf4::ncvar_get(nc, var)                 # lon x lat x time (NA = flagged missing)
  if (length(dim(a)) == 2) a <- array(a, c(dim(a), 1))
  a <- aperm(a, c(2, 1, 3))                      # lat x lon x time
  if (lat[1] < lat[length(lat)]) { a <- a[length(lat):1, , , drop = FALSE]; lat <- rev(lat) }
  dlon <- mean(diff(lon)); dlat <- abs(mean(diff(lat)))
  r <- terra::rast(a)
  terra::ext(r) <- c(min(lon) - dlon / 2, max(lon) + dlon / 2, min(lat) - dlat / 2, max(lat) + dlat / 2)
  terra::crs(r) <- "EPSG:4326"
  r
}
