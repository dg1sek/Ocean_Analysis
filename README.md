# Ocean_Analysis (Siegfried's fork): independent checks of the SST paper

**The paper and the main code live in Alf Loeffler's repository: <https://github.com/aloes2512/Ocean_Analysis>.**
This fork only holds supporting analyses that were written to test and extend that paper. Nothing here replaces
the paper's own code; the fork is kept in step with Alf's `master`, and the extra folders can be merged there if Alf wants them.

## What the paper is about

Alf's manuscript, *Modelling Solar-Locked Natural Variability in Global Ocean SST Anomalies (1850–2026)* (`Ocean_Analysis.qmd`),
asks whether the global sea surface temperature (ERSSTv6) can be split into

- a slow anthropogenic trend (a Richards/logistic growth curve), and
- natural oscillations tied to the sun (a spline with knots at the sunspot minima, filtered, then singular spectrum
  analysis and ESPRIT on the oscillating pairs).

This fork asks the next question: **do the oscillations survive independent tests?** The tests use different
data (coral and reconstruction records, other regions), different methods (maximum-likelihood sinusoid fits with
autocorrelated errors, surrogate-data significance, SSA, FFT) and a statistical null (red-noise surrogates).
The work was done jointly with Claude (Anthropic) and is meant as a transparent second pair of eyes, not as a verdict.

## What is in this fork

| Folder | Question | R code does |
|---|---|---|
| [`coral_test/`](coral_test/) | Is the ~20-yr line in global SST paced by an external clock (lunar nodal 18.6 yr, solar Hale cycle ~22 yr) or is it internal variability? A thermometer record of 176 yr cannot separate these; coral records of several centuries can. | Loads coral Sr/Ca and δ18O records (NOAA NCEI Paleoclimatology), fits the ~20-yr line in sliding windows and checks whether period and phase stay fixed. Figures: spectra, wavelet, sliding period/phase. |
| [`regional_test/`](regional_test/) | Does a *regional* analysis (Mediterranean and others, where sun angle and noise differ) confirm or weaken the global findings? | Screens 12 candidate regions by criteria fixed beforehand, then runs Alf's methods (ML sinusoid scan with AIC, SSA, FFT with 4 components, spline-at-solar-minima chain, sunspot regression, proxy check) for three regions: Mediterranean, Gulf of Mexico, western Pacific warm pool. |
| [`thermostat_test/`](thermostat_test/) | Is there a warm-water ceiling ("tropical thermostat") in tropical SST? Six-step plan; steps 0-3 done (no thermostat signature in ERSSTv6, OISST v2, HadSST4 or in NCEP reanalysis fluxes); step 4 judged pointless. | Pre-registered per-cell test on three SST products (warming rate against starting SST, upper-tail quantiles, area above 27-30 C) plus a heat-flux test (net flux against SST in warm and cool cells). |
| [`latitude_test/`](latitude_test/) | Does the 70-80 yr line and the eruption response depend on latitude (Alf, 23 Nov)? | Ocean-only ERSSTv6 latitude bands: ML sinusoid scan with surrogate null, eruption event study per band. Runs from its own folder. |
| [`volcano_test/`](volcano_test/) | What do volcanic eruptions say about how fast the upper ocean responds, and could the Sun produce the 80-yr line? | Event study of six eruptions with placebo test, one-box response model, AR6 conversion factor check, energy argument for a solar 60-120 yr line. Runs from its own folder with the data included. |
| [`sealevel_test/`](sealevel_test/) | Do tide gauges (Brest, Trieste, Genova, Marseille) independently show the 70-80 yr line? | PSMSL annual sea level: ML/OLS sinusoid scan with AR(1) surrogates (gap-aware), forecast test, phase at 75 yr compared with ERSST. Runs from its own folder with the data included. |
| [`harmonic_test/`](harmonic_test/) | Does the 60-100 yr harmonic forecast better than a plain trend, and how does the SSA window length change the component separation (Alf, 27 Nov)? | Forecast from cut-off years 1950/1975/2000 against AR(1) surrogates; SSA with window lengths N/2 to 40 yr and w-correlation. Runs from its own folder with the data included. |
| [`eof_test/`](eof_test/) | What are the dominant modes of the ocean SST field, and which one carries the 70-80 yr line (Alf, 27 Nov)? | EOF/PCA of ERSSTv6 annual means with North's rule, ML sinusoid scan with AR(1) null on each PC, per-cell line amplitude map, lagged PC cross-correlations. Needs the ERSST file (download script in `regional_test`). |
| [`slope_test/`](slope_test/) | Does the slope of the SST curve fall (Alf, 27 Nov)? | Quadratic trend plus ML sinusoid on global, MED, GOM, WPWP; acceleration with block bootstrap, 30-yr window slopes, size and phase of the sinusoid's own slope. Runs from its own folder with the data included. |

Each folder has its own `README.md` (what was done, how to re-run) and a results file with all numbers
(`coral_test/manuscript_section.md`, `regional_test/RESULTS.md`).

### Where the R code is used

- `regional_test/R/alf_tools.R`: reusable functions: `sin_scan` (ML sinusoid fit over periods 4–150 yr), `null_scan` and `null_scan_pl` (AR(1) and power-law surrogates for the significance of the best period), `ssa_periods`, `fft_top`, `knot_chain`, `solar_lag_fit`.
- `regional_test/R/regional_lib.R`: region boxes, ERSSTv6 loading, area-weighted regional anomalies against the 1991–2020 climatology, top-of-atmosphere insolation.
- `regional_test/10_…13_*.R`: screening, proxy export, the per-region analysis (`Rscript 12_regional_analysis.R MED`), and the check of the ~80-yr line against an Atlantic multidecadal index.
- `coral_test/A2_*.R`: loading (`A2_load`), library (`A2_lib`), run (`A2_run`) and post-processing (`A2_post`) of the coral experiment.
- `scripts/` and `methods_archive/` are Alf's original scripts and are not modified here.

## Headline results so far (details and caveats in the folders)

- A line of 74–80 yr (0.2–0.3 K) is found by the maximum-likelihood fit in the Mediterranean and the Gulf of Mexico, not in the warm pool; SSA of the annual series finds the same line.
- The spline-at-solar-minima-minus-LOESS curve is effectively a band-pass filter (gain 0.01 at 85 yr, about 0.5 at 22–24 yr), so the ~85-yr cycle cannot appear in it, and 22–24 yr periods also appear for random knots.
- No region shows a significant regression on the sunspot number; the 80-yr line coincides with an Atlantic multidecadal index, and one coral record does not confirm it.
- Limits: ERSSTv6 is statistically filled before about 1900, the noise model at 60–100 yr is uncertain, and the Mediterranean and Gulf of Mexico are not independent.

## Running it

R 4.x with `terra`, `ncdf4`, `Rssa`. ERSSTv6 (86 MB) is downloaded by `regional_test/00_download_ersst.R`; small proxy files are included.
Scripts use relative paths and are run from inside their folder.

## Data sources

NOAA PSL ERSSTv6 · SILSO sunspot number V2.0 · NOAA NCEI Paleoclimatology coral and reconstruction records
(sites listed in `regional_test/data/proxy_sites.csv`).

## Contact

Siegfried Loeffler (fork owner) and Alf Loeffler (paper author, [aloes2512/Ocean_Analysis](https://github.com/aloes2512/Ocean_Analysis)).
