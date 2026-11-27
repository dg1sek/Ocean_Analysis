# eof_test: dominant modes (EOF/PCA) of the ocean SST field (27 Nov 2026)
Alf's request: principal components of the global ocean as the basis for a coupled-oscillator model.
Run: `Rscript 60_eof.R` from this folder. It needs the ERSSTv6 file `data/ersst.v6.sst.mnmean.nc` (download with `../regional_test/00_download_ersst.R`, then copy or link it into `data/`).
Design: PREREG_eof.md. Results: RESULTS_eof.md. Tools: R/ (Alf's ML sinusoid scan, netCDF reader).
