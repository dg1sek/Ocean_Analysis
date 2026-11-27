# Pre-registration: dominant modes of the ocean SST field (EOF/PCA) -- written 2026-11-27 before any result

Why: Alf (27 Nov) asks for principal components / eigenmodes of the global ocean as the basis for a coupled-oscillator model.

Data: ERSSTv6 2-degree grid, monthly anomalies vs 1991-2020, annual means, cells without missing years, years 1900-2025 (main) and 1850-2025 (check; filled before ~1900).
Weights: sqrt(cos latitude). Method: SVD of the (year x cell) matrix.
Variants: V1 = anomalies as they are (first mode will contain the warming); V2 = global (area-weighted) mean removed from every year first.

Questions and rules
 E1 How many modes are separable? Variance fractions of PC1-8; North's rule: modes i and i+1 separable if lambda_i - lambda_(i+1) > lambda_i * sqrt(2/N*), N* = n (1-phi^2)/(1+phi^2) from the PC's AR(1). Only separable modes are interpreted.
 E2 Which mode carries the 60-100 yr line? ML sinusoid scan (alf_tools::sin_scan, trend order 2, 60-100 yr) on each of PC1-5; significant if the best AIC gain exceeds the 95th percentile of 300 AR(1) surrogates (null_scan).
 E3 Where does the line live? Per-cell amplitude of a sinusoid with period P* = the period found in E2 (if none, 76 yr, the Mediterranean optimum), with quadratic trend, map.
 E4 Link to known modes: correlation of PC1-5 with Nino3.4 (annual) at lag 0; descriptive.
 E5 Coupling hint: lagged cross-correlation between PC pairs (lags -15..+15 yr). Descriptive only; PCs are orthogonal at lag 0 by construction, no significance claimed.
No further tests added after seeing results; extras are labelled.
