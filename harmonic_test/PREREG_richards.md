# Pre-registration: rebuild of Alf's Richards + harmonic fit on ERSST Mediterranean, and a profile over K -- written 2026-11-28 before running

Source of the design: Alf's mail of 28 Nov 10:59 and his script `Richards_81.5harmonic.R` (github.com/aloes2512/Ocean_Basins). His data are Berkeley Earth SST, ours ERSSTv6, so this is a rebuild, not a re-run.
Series: MED annual means, ERSSTv6, 1850-2025 (n = 176). (His fit uses monthly values; annual means used here to keep AR(1) simple.)

Model (his formula): y = A + (K-A)/(1 + nu exp(-B (t-t0)))^(1/nu) + alpha cos(wt) + beta sin(wt), w = 2 pi/P.
Steps as he describes: (1) simple logistic (nu = 1) fixes B; (2) Richards with B fixed; (3) add harmonic, period scanned 60-100 yr in 0.1 yr steps, choose minimum RSS (allowed: same parameter count).
Bounds as in his script: A in [-3, 1], nu in [0.05, 10], t0 in [1850, 2076]; K in [0.2, 8] when K is free.

Questions and rules
 R1 Period: result and the interval of periods with AIC within 2 of the minimum (AIC = n log(RSS/n) + 2k).
 R2 K profile: K fixed at 2, 3, 4, 6, 8, 12, 20, 50 (all other steps repeated for each K, k = 7 free parameters incl. B and P). Report RSS and AIC difference to the best K. "K bounded by the data" only if AIC(K) exceeds the minimum by more than 2 on both sides of the best value within this grid; otherwise "K not determined". Also report the free-K fit and whether it sits on its bound (8).
 R3 Same data, quadratic trend + sinusoid (ML, AR(1)): AIC difference to the best Richards model (only a comparison of how well each fits; parameter counts stated).
 R4 Residuals of the best model: AR(1) coefficient; periodogram peak between 10 and 40 yr and its p against 300 AR(1) surrogates of the maximum peak (a test of Alf's "carrier" idea at the level of spectral peaks, not of its amplitude modulation).
No further tests are added after seeing the result.
