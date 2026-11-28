# Pre-registration: simulation with coupled stochastic oscillators (Alf's 28 Nov point, Siegfried's remark on record length) -- written 2026-11-28 before running

Question: if the true system is a set of coupled, damped, noise-driven oscillators with unrelated periods (no stable 75-yr line), how often does our detection procedure nevertheless report a 60-90 yr line, and how often does the "harmonic forecast" rule pass? That is the false-positive rate of what we did with the real data.

Truth model (n = 176 annual values, 1850-2025): y = quadratic trend (the OLS quadratic of the ERSST Mediterranean series) + sum of K = 3 stochastic oscillators + AR(1) background (phi = 0.3).
 Oscillator k: AR(2) process x_t = 2 r cos(w) x_(t-1) - r^2 x_(t-2) + e_t with period P_k = 2 pi / w and r = exp(-w/(2 Q_k)), Q_k uniform in [2, 10] (decay time from about 0.3 to 1.6 periods... i.e. clearly noisy, not a stable sinusoid).
 Variance: the oscillators share a fraction f uniform in [0.3, 0.7] of the total residual variance (equal shares), the AR(1) background the rest; total sd equals the residual sd of the real Mediterranean series about its quadratic.
Scenarios (400 series each): S0 pure AR(1) (check, expect about 5 % false positives); S1a periods log-uniform in [10, 150] but none inside 55-110 yr (re-drawn if so): variability exists but no oscillator near 75 yr; S1b periods log-uniform in [10, 150], unrestricted (some series will have an oscillator in the band by chance); S2 as S1b but one oscillator forced to P in [65, 85] (power check).
Procedures applied to each series exactly as to the real data:
 D1 Line: ML sinusoid scan 60-100 yr (alf_tools::sin_scan, quadratic trend, AR(1) errors), p against 100 surrogates from AR(1) fitted to the quadratic residuals ("null AR(1)") and, separately, from AR(2) ("null AR(2)"). "Flagged" if the best period is within 60-90 yr and p < 0.05.
 D2 Forecast rule: cut-offs 1950, 1975, 2000, quadratic vs quadratic + sinusoid (period from the training years), skill = 1 - RMSE1^2/RMSE0^2; the rule passes if skill exceeds the fixed thresholds taken from the real Mediterranean surrogates (0.53, 0.50, 0.39) at >= 2 of 3 cut-offs.
Reported: flagged fraction per scenario and null (the false-positive rate), the fraction passing D2, and for the real Mediterranean series the p-values under AR(1) and AR(2) nulls (300 surrogates).
Reading rule fixed in advance: if the false-positive rate under S1a/S1b exceeds 10 % for the AR(1) null, the AR(1) p-values of the real data are called optimistic; if it also exceeds 10 % for the AR(2) null, the line cannot be called distinguishable from what such a system generates by itself.
No further scenarios will be added after seeing the result.
