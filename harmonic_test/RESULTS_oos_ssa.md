# Results (run 2026-11-27; rules in PREREG_oos_ssa.md, written before the run)

## A. Forecast without knowledge of the future (`out/oos_skill.csv`)
Skill = 1 - RMSE(trend+sinus)^2 / RMSE(trend)^2 over the forecast years; null = 300 surrogates without harmonic.

| series | cut-off | period found on training (yr) | skill | null 95% | p |
|---|---|---|---|---|---|
| GLOBAL | 1950 | 61 | 0.57 | 0.65 | 0.057 |
| GLOBAL | 1975 | 66 | -0.09 | 0.51 | 0.38 |
| GLOBAL | 2000 | 67 | -0.63 | 0.46 | 0.73 |
| MED | 1950 | 69 | 0.66 | 0.53 | 0.030 |
| MED | 1975 | 79 | -1.29 | 0.50 | 0.80 |
| MED | 2000 | 72 | 0.67 | 0.39 | 0.000 |
| GOM | 1950 | 60 | 0.42 | 0.54 | 0.083 |
| GOM | 1975 | 68 | -0.93 | 0.45 | 0.73 |
| GOM | 2000 | 77 | 0.53 | 0.38 | 0.017 |

**Corrected 2026-11-28:** the first run (reported to Alf in mail 18) used a partial year 2026 (January-October) as the 177th annual value although the text said 1850-2025. All scripts now use complete years 1850-2025 (n = 176). Effect: the verdict for GOM changes from "met" to "not met" (1950: skill 0.42 below the null level 0.54); MED stays "met", GLOBAL "not met".

Pre-registered rule (skill above the null 95% level at >= 2 of 3 cut-offs): **met for MED only; not met for GOM and GLOBAL**.
Reading: the verdict is mixed, not clear. At cut-off 1975 the harmonic forecast was much worse than the plain trend in MED and GOM (skill -1.2 and -0.9), at 1950 and 2000 it was better (clearly so for MED). The three forecasts overlap (the 1950 forecast contains the later windows), so they are not independent evidence. 9 tests were made.
Caveat that favours the harmonic: the surrogates assume a quadratic trend is the truth. If the true trend is not quadratic, the sinusoid can absorb the misfit and look like forecast skill. The result therefore does not show that the line is physical. A check with a different trend model (for example a forcing-based trend as in research/W2) has not been done and would have to be registered first.

## B. SSA window length (`out/ssa_windows.csv`), maximum |w-correlation| among components 1-3
| L (yr) | GLOBAL | MED |
|---|---|---|
| 88 (N/2) | 0.07 | 0.70 |
| 80 | 0.08 | 0.50 |
| 60 | 0.31 | 0.68 |
| 58 (N/3) | 0.32 | 0.67 |
| 44 (N/4) | 0.34 | 0.46 |
| 40 | 0.35 | 0.36 |

Alf's expectation is confirmed: with L = N/3 the components 1-3 are strongly correlated for MED (0.67). Separation below the threshold of 0.3 is reached only globally with L = 80 and 88; for MED no window gets below 0.3 (best: L = 40 with 0.36, L = 80 with 0.50).
The first component carries 83-84 % (MED) and 92-94 % (global) in all windows; components 2 and 3 carry only 3-8 %.
**Not usable:** the column `pair_60_100` in the output is an artefact. For L >= 60 my period estimator gives the window length itself for the slow trend components (88,88,88,88 and 80,80,80,80), so "TRUE" does not mean an oscillatory pair of 60-100 yr was found. No oscillatory pair of that length is shown by this diagnostic.
A pair near 20-22 yr appears in most MED windows (L = 88, 60, 44, 40, 80); the periods come from the FFT of the eigenvectors and are coarse. I draw no conclusion from it.
