# Alf's analysis in three regions (ERSSTv6, annual means 1850-2025, anomalies against 1991-2020)

Regions: MED (Mediterranean), GOM (Gulf of Mexico), WPWP (western Pacific warm pool). Selection: `REGIONS.md`.
All numbers come from `out/*_log.txt` and the CSV files next to them; scripts `10`-`13`. Seeds fixed; 400 noise simulations per test, 150 random knot sets.

## What was done (Alf's own steps)
1. **Maximum-likelihood sinusoid** a*sin(wt) + b*cos(wt) + smooth trend, AR(1) errors, scanned over periods 4-150 yr, AIC gain = n log(RSS0/RSS1) - 6 (a, b and the period are fitted). Trend: polynomial of order 1, 2 (main) or 3.
   Significance: how large is the *best* AIC gain over all periods in noise alone? (AR(1) surrogates, and power-law surrogates whose spectrum is fitted to 2-30 yr only.)
2. **SSA** of the annual series (window L = N/3, also 30/58/70/88) and ESPRIT periods; eigenvalues against 150 AR(1) surrogates.
3. **FFT**: sum of the 4 strongest Fourier components.
4. **His knot chain** on the monthly series: cubic spline with knots at the solar minima minus LOESS (span 0.25), compared with 150 random knot sets with the same spacing; SSA (L = 700) of that curve.
5. **Solar tests**: regression of the annual anomaly on lagged annual sunspot number (lags 0-6 yr), also with the Hale-signed sunspot number (22 yr).
6. **Independent proxy**: the same scan on the coral/reconstruction record inside the region where one exists (GOM: Dry Tortugas coral 1734-2008; WPWP: Freund warm-pool index 1617-2007).

## Result 1: the 80-yr line is real in the data of MED and GOM, absent in WPWP

| Region | best period (yr), trend order 2 | amplitude (K) | AIC gain | p vs AR(1) | p vs power law | 85-yr fit, amplitude +- se (order 1 / 2 / 3) |
|---|---|---|---|---|---|---|
| MED | 73.9 | 0.30 | 13.5 | 0.007 | 0.017 | 0.33 +- 0.08 / 0.25 +- 0.05 / 0.29 +- 0.03 K |
| GOM | 80.1 | 0.23 | 17.1 | 0.002 | 0.005 | 0.28 +- 0.05 / 0.22 +- 0.03 / 0.24 +- 0.03 K |
| WPWP | 6.8 (noise) | 0.06 | 2.9 | 0.82 | 0.81 | 0.15 +- 0.08 / 0.04 +- 0.03 / 0.05 +- 0.03 K |

- Order 1 / 3 give 76.6 / 80.9 yr (MED, p = 0.065 / 0.005) and 78.7 / 83.9 yr (GOM, p = 0.010 / 0.005). The amplitude at a fixed 85 yr moves by up to 30 % with the trend choice (about two cycles in 176 yr).
- Noise alone gives a best AIC gain of 9-11 at the 95 % level (`*_log.txt` section A), so the MED and GOM peaks are above, WPWP is far below.
- No other period is supported: fixed 40, 24, 22 and 11 yr give AIC gains of -3.8 to -0.7 (WPWP 24 yr: -0.9, 22 yr: -1.9, order 2).
- FFT, 4 strongest components (quadratic-detrended): MED 88, 58.7, 29.3, 44 yr (49 % of the variance of the sum); GOM 88, 58.7, 8.4, 35.2 yr (42 %); WPWP 58.7, 25.1, 6.8, 19.6 yr (26 %).

## Result 2: SSA sees the same line when it is applied to the series itself

SSA of the quadratic-detrended annual series: components 1 and 2 form one pair, above the AR(1) 95 % limit (MED 29.4 and 24.7 % vs 22.0 and 17.1 %; GOM 26.5 and 18.3 % vs 14.5 and 12.2 %; WPWP: no component above its limit).
ESPRIT period of that pair: MED 78.8 yr (54 % of variance), GOM 80.1 yr (45 %), WPWP 57 yr (17 %, not significant).
The pair does not depend on the window: MED 80.0 / 78.8 / 80.2 / 82.4 yr and GOM 77.7 / 80.1 / 79.6 / 80.3 yr for L = 30 / 58 / 70 / 88 (`*_ssa_window.csv`).
Its rms amplitude equals the ML sinusoid (MED 0.209 vs 0.214 K, GOM 0.151 vs 0.161 K at L = 58). So SSA and maximum likelihood give the same amplitude here; no adjustment of a window or "K" parameter is needed for this series.

## Result 3: why the knot chain cannot show an 85-yr cycle

His "harmonic part" is the spline at the solar minima minus LOESS(0.25). Applied to unit sinusoids it passes (`*_construction_gain.csv`):

| period (yr) | 120 | 85 | 60 | 40 | 30 | 24 | 22 | 18 | 11 | 5 |
|---|---|---|---|---|---|---|---|---|---|---|
| gain | 0.00 | 0.01 | 0.03 | 0.13 | 0.36 | 0.56 | 0.49 | 0.15 | 0.05 | 0.02 |

It is a band-pass filter centred on 22-24 yr. Slow oscillations are removed with the LOESS trend, so an 85-yr line cannot appear in it, and SSA of that curve can only return periods between about 18 and 40 yr.
Consistently, in all three regions the dominant period of the curve is 22.8 / 23.5 / 24.4 yr (MED / GOM / WPWP) with his minima and 22.8 / 24.0 / 24.4 yr (median) with random knots (5-95 %: 19.5-25.8, 19.1-24.8, 21.3-25.8). The rms amplitude with his knots is within the random range in MED and GOM (shares of random sets with larger rms: 0.22 and 0.83) and just above it in WPWP (0.05). SSA of the curve gives pairs at 22-26 yr with modulus about 1 for his knots and for the random sets alike.

## Result 4: no solar signal found in any region

- Sunspot regression (AR(1) errors): smallest p over lags 0-6 is 0.208 (MED, lag 2), 0.163 (GOM, lag 1), 0.273 (WPWP, lag 2); all 1.0 after Bonferroni over 7 lags.
- Hale-signed regressor: smallest p 0.140 (MED), 0.104 (GOM), 0.044 (WPWP, lag 4; 0.053 and 0.056 at lags 3 and 2). One uncorrected p below 0.05 among 3 regions x 7 lags x 2 regressors is what chance produces; we count it as a hint to follow, not as evidence.
- The sunspot record itself (1749-2025) has no 60-120 yr sinusoid: best period 97 yr, AIC gain -3.0 (`13_gleissberg_checks.R`).

## Result 5: the 80-yr line looks Atlantic, not solar

- AMV proxy from the same grid (0-60N, 80W-0E minus global mean): best period 80 yr, amplitude 0.13 K, AIC gain 6.9.
- Most recent maximum of an 85-yr sinusoid: AMV 1941, MED 1946, GOM 1950; WPWP amplitude 0.04 K (`gleissberg_phase_table.csv`). Both Atlantic-side regions swing with the AMV proxy; the Pacific region does not.
- Independent coral (GOM, Dry Tortugas, 1734-2008): its own best period in 4-137 yr is 4.3 yr (p = 0.46). A fixed 80-yr sinusoid has amplitude 0.42 sd and AIC gain 8.5 (not calibrated against noise) but its phase over 1850-2008 is 111 deg against 337 deg for ERSSTv6 (about half a cycle apart), and the high-passed correlation of coral and ERSSTv6 is 0.09. It does not confirm the line. It is one site against a regional mean.
- WPWP reconstruction (Freund warm-pool index): best period 5.8 yr, p = 0.22; correlation with ERSSTv6 -0.12. It carries no usable low-frequency signal here.

## Limits
- ERSSTv6 is statistically filled before about 1900; regional series are less reliable than the global one. HadSST4 (not filled) gridded data would be the check; not done.
- Two cycles in 176 yr. The p-values depend on the noise model at periods of 60-100 yr, which neither AR(1) nor a 2-30 yr power law constrains. Regional piControl runs would be the better null; we only have global means from three CMIP6 models.
- MED and GOM are correlated; they are close to one Atlantic result, not two.
- Three regions were analysed in full; the analysis was chosen before the results, but the scan over periods, the three trend choices and the three regions are separate looks.
- Reproduction check on the MacBook (R 4.6.1, terra 1.9.60, ncdf4 1.24): scripts 01-05 of the earlier global work reproduce the Windows output (global series to 1e-15, 11,074 valid ERSSTv6 cells in every month).
