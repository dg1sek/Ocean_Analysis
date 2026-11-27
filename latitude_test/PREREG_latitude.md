# Latitude bands of the ocean (ERSSTv6): pre-registered design, written 2026-11-27 before any result
Motivation: Alf (mail of 23 Nov 2026, 11:46): oscillations of the ocean depend on latitude because the imbalance between absorbed solar and emitted power depends on latitude; "singular events" can be quantified from the ERSST 2-degree grid.
Bands (ocean cells only, area-weighted means of monthly anomalies, base 1991-2020): 60-40S, 40-20S, 20S-0, 0-20N, 20-40N, 40-60N. Fixed now.
A. Event study for the six eruptions, identical settings to volcano_test (trend by LOESS, ENSO regression, seasonal-cycle removal, placebo with matching calendar months, one-box memory time and gain), per band.
B. Maximum-likelihood sinusoid scan (periods 4-150 yr, quadratic trend, AR(1) errors) on the annual band series 1850-2025, with the AR(1)-surrogate null (300 simulations): best period, amplitude, AIC gain, p.
Reading: a latitude dependence is claimed only where intervals/p-values separate bands; all six bands are reported. Limits: band means mix basins (Atlantic and Pacific signals can cancel); 60-40S and 40-60N have few ocean cells in some longitudes and ERSSTv6 early data are filled.
