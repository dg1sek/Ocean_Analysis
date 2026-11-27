# Results (run 2026-11-27; rules in PREREG_oos_ssa.md, written before the run)

## A. Forecast without knowledge of the future (`out/oos_skill.csv`)
Skill = 1 - RMSE(trend+sinus)^2 / RMSE(trend)^2 over the forecast years; null = 300 surrogates without harmonic.

| series | cut-off | period found on training (yr) | skill | null 95% | p |
|---|---|---|---|---|---|
| GLOBAL | 1950 | 61 | 0.55 | 0.56 | 0.053 |
| GLOBAL | 1975 | 66 | 0.02 | 0.53 | 0.33 |
| GLOBAL | 2000 | 67 | -0.39 | 0.44 | 0.65 |
| MED | 1950 | 69 | 0.61 | 0.58 | 0.047 |
| MED | 1975 | 79 | -1.23 | 0.46 | 0.78 |
| MED | 2000 | 72 | 0.64 | 0.46 | 0.003 |
| GOM | 1950 | 60 | 0.43 | 0.35 | 0.043 |
| GOM | 1975 | 68 | -0.93 | 0.45 | 0.77 |
| GOM | 2000 | 77 | 0.57 | 0.37 | 0.010 |

Pre-registered rule (skill above the null 95% level at >= 2 of 3 cut-offs): **met for MED and GOM, not met for GLOBAL**.
Reading: the verdict is mixed, not clear. At cut-off 1975 the harmonic forecast was much worse than the plain trend in MED and GOM (skill -1.2 and -0.9), at 1950 and 2000 it was clearly better. The three forecasts overlap (the 1950 forecast contains the later windows), so they are not independent evidence. 9 tests were made.
Caveat that favours the harmonic: the surrogates assume a quadratic trend is the truth. If the true trend is not quadratic, the sinusoid can absorb the misfit and look like forecast skill. The result therefore does not show that the line is physical. A check with a different trend model (for example a forcing-based trend as in research/W2) has not been done and would have to be registered first.

## B. SSA window length (`out/ssa_windows.csv`), maximum |w-correlation| among components 1-3
| L (yr) | GLOBAL | MED |
|---|---|---|
| 88 (N/2) | 0.07 | 0.65 |
| 80 | 0.08 | 0.41 |
| 60 | 0.31 | 0.66 |
| 58 (N/3) | 0.32 | 0.65 |
| 44 (N/4) | 0.33 | 0.44 |
| 40 | 0.35 | 0.35 |

Alf's expectation is confirmed: with L = N/3 the components 1-3 are strongly correlated for MED (0.65). Separation below the threshold of 0.3 is reached only globally with L = 80 and 88; for MED no window gets below 0.3 (best: L = 40 with 0.35, L = 80 with 0.41).
The first component carries 83-84 % (MED) and 92-94 % (global) in all windows; components 2 and 3 carry only 3-8 %.
**Not usable:** the column `pair_60_100` in the output is an artefact. For L >= 60 my period estimator gives the window length itself for the slow trend components (88,88,88,88 and 80,80,80,80), so "TRUE" does not mean an oscillatory pair of 60-100 yr was found. No oscillatory pair of that length is shown by this diagnostic.
A pair near 20-22 yr appears in most MED windows (L = 88, 60, 44, 40, 80); the periods come from the FFT of the eigenvectors and are coarse. I draw no conclusion from it.
