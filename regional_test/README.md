# regional_test: Alf's analysis in three ocean regions

What this folder is: the regional analysis you asked for, run with your own methods, in three regions chosen beforehand
by stated criteria. Everything is base R plus `terra`, `ncdf4` and `Rssa`; every number in `RESULTS.md` can be
regenerated with the scripts below.

## Contents

| File | What it does |
|---|---|
| `REGIONS.md` | the 12 candidate regions, the criteria, the scores and the choice of three (WPWP, GOM, MED) |
| `RESULTS.md` | the results with all numbers (English, can be used as a manuscript section) |
| `10_region_screening.R` | screening of the candidate regions (written before any signal was looked at) |
| `12_regional_analysis.R <ID>` | your steps for one region: maximum-likelihood sin/cos fit with AIC, SSA, FFT with 4 components, knot chain, sunspot regression, coral check |
| `13_gleissberg_checks.R` | is the 80-yr line the Atlantic multidecadal oscillation? compares with an AMV proxy and the sunspot record |
| `11_export_proxies.R` | how `data/proxy_annual.csv` was made from the NOAA coral files of the coral test (only runs in the original working tree; the CSV is included) |
| `R/` | helper functions (`alf_tools.R` holds the analysis steps as short functions) |
| `data/` | proxy sites, annual proxy records, SILSO sunspot numbers (small). ERSSTv6 is downloaded by `00_download_ersst.R` |
| `out/` | all results: tables (`.csv`), logs (`*_log.txt`), figures (`*_fig1_series_fits.png`) |

## How to run (R 4.x)

```r
install.packages(c("terra", "ncdf4", "Rssa"))
setwd("regional_test")          # run everything from this folder
source("00_download_ersst.R")   # 86 MB, once
```
then in a terminal: `Rscript 10_region_screening.R`, `Rscript 12_regional_analysis.R MED` (also `GOM`, `WPWP`), `Rscript 13_gleissberg_checks.R`.
About one minute per region.

## Result in four sentences
1. In the Mediterranean and the Gulf of Mexico the maximum-likelihood fit finds a line at 74-80 yr (amplitude 0.2-0.3 K, above what red noise gives); the warm pool has none.
2. SSA of the annual series finds the same pair (79-80 yr) with the same amplitude, so the two methods agree.
3. The spline-at-solar-minima-minus-LOESS curve can only contain periods of roughly 18-40 yr (its gain at 85 yr is 0.01), which is why the 85-yr cycle is missing there, and why 22-24 yr appears there for random knots too.
4. No region shows a significant link to the sunspot number; the 80-yr line peaks together with an Atlantic multidecadal index and is not confirmed by the Dry Tortugas coral.

Data: NOAA ERSSTv6 (PSL), SILSO sunspot number V2.0, NOAA NCEI Paleoclimatology coral and reconstruction records (see `data/proxy_sites.csv`).
