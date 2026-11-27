# volcano_test: what do volcanic eruptions tell us about the inertia of the ocean, and can the Sun produce an 80-year line?

Everything here runs from this folder with the data in `data/` (about 1.3 MB; no downloads). Base R only, plus package `ncdf4` for the solar script.

| File | What it does |
|---|---|
| `20_volcano_events.R` | event study of six eruptions (Krakatau 1883, Santa Maria 1902, Novarupta 1912, Agung 1963, El Chichon 1982, Pinatubo 1991): composite SST response (global and four regions) with a placebo test, and a one-box response model (memory time tau, gain) fitted to the AR6 stratospheric aerosol record |
| `31_ar6_factor_check.R` | where the factor "-20 W/m2 per unit aerosol optical depth" comes from, and how much the conversion matters |
| `30_solar_longperiod.R` | strength of a 60-120 yr line in solar forcing (NOAA TSI 1610-2025, AR6 solar ERF) and the SST amplitude it could produce (energy argument) |
| `PREREG_volcano.md` | design and decision rules, written before the numbers (with the amendment made after the first run) |
| `RESULTS_volcano.md` | all results with numbers and limits |
| `out/` | tables, log files, figure `volcano_fig_composites.png` |
| `data/` | input series (see below) |

## How to run (R 4.x, RStudio)
1. Open the folder `volcano_test` in RStudio (or `setwd()` to it).
2. `source("20_volcano_events.R")` (about 1 minute), then `source("31_ar6_factor_check.R")`.
3. For the solar script install the package once: `install.packages("ncdf4")`, then `source("30_solar_longperiod.R")`.
Results are written to `out/`.

## Data in `data/`
- `global_series.csv`: global mean monthly SST anomaly (ERSSTv6, several variants; column `ersst` is used)
- `regional_monthly_ersst.csv`: monthly regional means of the regions of `regional_test` (columns MED, GOM, WPWP, ...)
- `nino34.long.anom.data`: HadISST Nino3.4 anomaly (NOAA PSL)
- `ar6_volcanic_sAOD_monthly_1850-2019.csv`, `ar6_volcanic_erf.csv`, `ar6_solar_erf.csv`: IPCC AR6 Chapter 7 forcing series (github.com/chrisroadmap/ar6)
- `tsi_v03r00_yearly_s1610_e2025_c20260305.nc`: NOAA Climate Data Record of yearly total solar irradiance
