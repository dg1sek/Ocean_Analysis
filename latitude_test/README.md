# latitude_test: does the 70-80 yr line and the response to eruptions depend on latitude?

Prompted by Alf's mail of 23 Nov 2026 (oscillations of the ocean depend on latitude; singular events can be quantified from the ERSST 2-degree grid).
Ocean-only ERSSTv6, six 20-degree latitude bands from 60S to 60N. Results: `RESULTS_latitude.md`. Design fixed beforehand: `PREREG_latitude.md`.

How to run (R 4.x, RStudio; open this folder, `setwd()` to it): `source("41_latitude_sinus_scan.R")` and `source("40_latitude_events.R")`. The monthly band series are included in `data/`
(`latband_monthly_ersst.csv`), so no download is needed. `39_make_latbands.R` rebuilds that file from the ERSSTv6 grid (86 MB, see `../regional_test/00_download_ersst.R`; packages `terra`, `ncdf4`).
`R/alf_tools.R` holds Alf's analysis steps as functions (same file as in `regional_test`).
