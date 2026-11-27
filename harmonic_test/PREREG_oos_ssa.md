# Pre-registration: (A) out-of-sample test of the harmonic line, (B) SSA window lengths -- written 2026-11-27, before running anything

Origin: Alf's mails of 27 Nov (19:17: "Fortsetzung der harmonischen Veränderungen ... Was ist richtig?"; first mail: window K = floor(N/3), other lengths not tried).
Data: annual ERSSTv6 means 1850-2025 for GLOBAL, MED, GOM (as in slope_test).

## A. Forecast without knowledge of the future
Models, both with AR(1) errors (alf_tools ML): M0 = quadratic trend; M1 = quadratic trend + one sinusoid, period scanned 60-100 yr (step 1) on the TRAINING years only.
Cut-off years Y = 1950, 1975, 2000 (training 1850..Y). Forecast Y+1..2025 from the deterministic part (trend + sinusoid), no AR(1) correction.
Score: skill = 1 - RMSE(M1)^2 / RMSE(M0)^2 over the forecast years (positive = the harmonic improves the forecast).
Null: 300 surrogates per series and cut-off, generated from the M0 fit plus AR(1) noise (no harmonic exists in them); the identical procedure (period scan on training, forecast) gives the null distribution of skill.
Decision: the harmonic "has forecast value" for a series if the observed skill exceeds the 95th percentile of its null at >= 2 of the 3 cut-offs. Otherwise "no forecast value shown".
Caveats fixed in advance: ERSST before about 1900 is a reconstruction; 3 cut-offs are not independent; the 1950 cut-off tests the plateau of 1940-70 only through the forecast, not through its fit.

## B. SSA window length (series MED and GLOBAL, annual, 1850-2025, N = 176)
Windows L = N/2 (88), N/3 (58), N/4 (44), and 40, 60, 80 yr. Rssa basic SSA.
Reported per L: share of the first 8 eigen-triples, dominant period of each (from the eigenvector), weighted correlation (w-correlation) among components 1-8, and whether a pair with period 60-100 yr appears.
Alf's criterion: a window length is better if components 1-3 are less strongly w-correlated. Threshold for "separated": |w-corr| < 0.3 between the trend group and an oscillatory pair.
No significance test (no red-noise Monte Carlo) in this step; it is a diagnostic of separability, nothing else.
