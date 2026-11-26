# thermostat_test: is there a warm-water ceiling ("tropical thermostat") in the SST data?

Status 2026-11-26: steps 0 and 1 of a six-step plan are done.

| Step | What | Status |
|---|---|---|
| 0 | literature: what does the thermostat idea claim, who disputes it | done (summary in the email to Alf) |
| 1 | ERSSTv6 per-cell check: does warming slow where the water is already warmest? | done, `10_step1_percell.R`, `out/` |
| 2 | repeat with independent SST products (HadSST4 unfilled, OISST 1982+) | proposed |
| 3 | surface-flux data (CERES cloud effect, ERA5/OAFlux latent heat): sign of the feedback above a threshold | proposed |
| 4 | small regional energy-balance model with a threshold loss term T*; is T* identifiable? | proposed |
| 5 | put it into the long-term (forcing-driven) model in place of the logistic ceiling | proposed |
| 6 | report | |

`PREREG_step1.md` fixes the design and the decision rule before any number was computed. `out/step1_log.txt` has all numbers,
`out/step1_*.csv` the tables, the two PNGs the figures.

Run (R 4.x, terra + ncdf4): `Rscript 10_step1_percell.R` from this folder. It needs `ersst.v6.sst.mnmean.nc`
(see `../regional_test/00_download_ersst.R`). About two minutes.

## Result of step 1 (ERSSTv6, tropics |lat| <= 30, 4231 complete cells)
Pre-registered verdict: **signature not seen.**
- Warm cells do not warm more slowly: slope of trend against starting SST (cells >= 27 C) is +0.013 K/decade per K (95 % 0.008 to 0.024) in 1950-2025, positive also in 1980-2025 and 1950-2019.
- The top of the SST distribution warms as fast as or faster than the middle: trend(99 % quantile) - trend(median) = +0.016 K/decade (0.008 to 0.029).
- The area warmer than 30 C grew from about 0.2 % (1950) to 5.6 % (2025) of the tropical ocean, and warmer than 29 C from 5 % to 22 %.
Caveats: ERSSTv6 is a smoothed 2-degree reconstruction; a ceiling of dynamic origin cannot be separated from a thermostat in this test; the warmest bins hold few cells.
