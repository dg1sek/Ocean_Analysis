# Pre-registration: the slope of the SST curve (written 2026-11-27, before any result)

Alf's question (27 Nov): has the slope been examined? A logistic ceiling shows up first as a falling slope.

Series: annual means of ERSSTv6 for GLOBAL (global_series.csv, column `ersst`), MED, GOM, WPWP (regional_monthly_ersst.csv).
Periods: 1850-2025 and 1900-2025 (before about 1900 ERSST is statistically filled).

Model (Alf's own): quadratic trend + one sinusoid, period scanned 60-100 yr (step 1), ML with AR(1) errors (alf_tools.R `sin_scan`).
Time centred on 1988: y = a + b x + c x^2 + sinusoid, x = (t-1988)/100.

Quantities
 Q1  acceleration c (K per century squared): slope of the trend at time t is b + 2 c x.
     Deceleration claimed only if the 95% interval of c lies fully below 0; acceleration if fully above 0; otherwise "undetermined".
 Q2  30-yr window slopes (OLS, K/decade) for 1936-65, 1966-95, 1996-2025, for the raw series and for the series minus the fitted sinusoid.
     Compared: 1996-2025 minus 1966-95. Same rule as Q1.
 Q3  Size of the sinusoid's own slope: 2 pi A / P (K/decade) against the trend slope in 1996-2025. Descriptive.
Uncertainty: moving-block bootstrap of the model residuals (5-yr blocks, 500 replicates), whole procedure (period scan, fit, windows) repeated per replicate.
No further tests are added after seeing the result; any extra check is labelled as such.
