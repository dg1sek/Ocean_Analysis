# Results: rebuild of Alf's Richards + harmonic fit on ERSST Mediterranean (run 2026-11-28; rules in PREREG_richards.md)

Data: ERSSTv6 2-degree annual means, box 6W-36E, 30-46N (same box and same NOAA product as Alf's `Medi_anomaly.R`, `ext(-6, 36, 30, 46)`; his weights are cell area times ocean fraction, with a land-fraction mask below 10 %, ours the native ERSST ocean cells). 1850-2025, n = 176. Tables: `out/richards_*.csv`; figure `out/richards_fig.png`.

## R1 / free K
Free K (bound 0.2-8): the minimum RSS is at P = 75.5 yr, **K = 8, i.e. on its upper bound**; nu = 2.8, t0 = 2062, B = 0.109, harmonic amplitude 0.25 K. Alf reports 75.8 yr and K = 8 from monthly NOAA data to Oct 2026.
Periods with AIC within 2 of the minimum: **73.4 - 77.8 yr**.

## R2 K profile (K fixed, everything else refitted; 7 free parameters)
| K fixed | P (yr) | amplitude (K) | nu | t0 | AIC - minimum |
|---|---|---|---|---|---|
| 2 | 75.8 | 0.255 | 2.36 | 2030 | 0.95 |
| 3 | 75.7 | 0.254 | 2.07 | 2039 | 0.54 |
| 4 | 75.6 | 0.254 | 1.93 | 2046 | 0.35 |
| 6 | 75.6 | 0.253 | 1.80 | 2056 | 0.17 |
| 8 | 75.6 | 0.253 | 1.74 | 2063 | 0.08 |
| 12 | 75.5 | 0.253 | 1.68 | 2073 | 0.00 |
| 20 | 75.0 | 0.236 | 1.36 | 2076 (bound) | 2.64 |
| 50 | 71.2 | 0.214 | 0.96 | 2076 (bound) | 19.65 |
**K is not determined by the data**: every K from 2 to 12 K fits equally well (AIC differences below 1). The larger values lose only because t0 is capped at 2076 in Alf's bounds, not because the data exclude them. The period stays at 75-76 yr whatever K is.

## R3 against quadratic trend + sinusoid (same data, OLS, 6 parameters)
Quadratic + sinusoid (P = 74.1): RSS 7.48, AIC -543.9. Richards + harmonic (K free, 8 parameters): RSS 5.95, AIC -580.2. The Richards form fits **much better** (AIC difference 36 with two more parameters), i.e. the flat-then-steep shape is closer to the data than a quadratic. Consequence: the quadratic baseline in the forecast test (RESULTS_oos_ssa.md) is probably too weak and favours the harmonic; the forecast test should be repeated with a Richards trend (not yet done, would be pre-registered first).

## R4 residuals of the best model
AR(1) coefficient 0.29. Largest periodogram peak between 10 and 40 yr at 22.0 yr, p = 0.067 against 300 AR(1) surrogates: **not significant at 5 %**. This does not test the amplitude modulation Alf describes, only whether the residual has a spectral line.

## Reading
1. Alf's results are reproduced on our data: period 75.5 vs 75.8 yr, K at the bound.
2. K = 8 is the upper limit of the search range and not an estimate; this is the same statement as "K cannot be derived from the data".
3. Over the data range the Richards curve is convex and without sign of saturation (K far beyond the data).
