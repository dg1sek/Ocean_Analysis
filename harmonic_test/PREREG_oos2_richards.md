# Pre-registration, part 2 of the forecast test: Richards trend instead of a quadratic -- written 2026-11-28 before running

Reason: the rebuild of Alf's fit (RESULTS_richards.md, R3) showed that a Richards curve fits the Mediterranean much better than a quadratic (AIC 36 lower). The quadratic baseline of the first forecast test (PREREG_oos_ssa.md, part A) is therefore probably too weak and favours the harmonic.

Same design as part A, except for the models (annual ERSSTv6 means, complete years 1850-2025; GLOBAL, MED, GOM; cut-offs 1950, 1975, 2000; forecast Y+1..2025; skill = 1 - RMSE(M1)^2/RMSE(M0)^2; decision: skill above the null 95% level at >= 2 of 3 cut-offs):
 M0 = Richards trend, y = A + (K-A)/(1 + nu exp(-B (t-t0)))^(1/nu); K fixed at 8 (the value Alf's bounds give; the K profile in RESULTS_richards.md showed the fit does not depend on K between 2 and 12); B from the logistic step (nu = 1) on the training years; A, nu, t0 free (bounds as in Alf's script).
 M1 = M0 + alpha cos(wt) + beta sin(wt); period scanned 60-100 yr in 1-yr steps on the TRAINING years only, joint refit (A, nu, t0, alpha, beta) for each period.
 Null: 100 surrogates per series and cut-off: M0 fitted on the training years, extended over all years, plus AR(1) noise (phi, sd from the training residuals); the same procedure (M0 and M1 fit on training, forecast) on each surrogate. Surrogates in which a fit does not converge are dropped and counted.
 If M0 or M1 cannot be fitted on the real training data for a cut-off, that cut-off counts as "not fittable" and is reported as such (no fallback to another model).
Caveats fixed in advance: the three forecasts overlap; 9 tests; ERSST before 1900 is a reconstruction; the 1950 training window contains almost no steep part, so a Richards fit there is weakly constrained (K is fixed for this reason).
