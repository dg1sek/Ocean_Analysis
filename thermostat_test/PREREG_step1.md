# Tropical thermostat, step 1: pre-registered design (written 2026-11-26 before any result was computed)

Question: does ERSSTv6 show the signature of a warm-water ceiling (fewer degrees of warming where the water is already warmest)?

Data: ERSSTv6 monthly SST, cells with |lat| <= 30, ocean cells with no missing month in the analysis period.
Periods (all three reported): P1 1950-2025 (main), P2 1980-2025, P3 1950-2019 (excludes the 2020-25 El Nino/La Nina swings).
Trends: ordinary least squares on annual means, K per decade. Uncertainty: moving-block bootstrap over years (block 5 yr, 1000 resamples,
all cells resampled together so spatial correlation is kept). Area weights cos(lat).

Test A (trend vs warmth). Predictor = each cell's mean annual SST over the first 30 yr of the period (P1, P3: 1950-79; P2: 1980-99).
 A1 area-weighted mean trend per 1 K predictor bin (table + figure).
 A2 slope of cell trend against predictor for cells with predictor >= 27 C (K/decade per K) with bootstrap 95 % interval.
Test B (upper tail). Per year, pooled area-weighted monthly SST distribution of the tropical ocean: quantiles 50, 90, 99 %.
 B1 trend(q99) - trend(q50) in K/decade with bootstrap 95 % interval.
 B2 area fraction above 27, 28, 29, 30 C: trend per decade (relative to its mean) -- reported, not scored.

Decision rule (fixed now):
 SUPPORTS a thermostat: A2 slope < 0 with 95 % interval excluding 0 AND B1 < 0 with 95 % interval excluding 0, in P1 and in P2.
 REFUTES / not seen: both A2 and B1 intervals include 0 or are positive in P1.
 INCONCLUSIVE: anything else (for example only one of the tests, or P1 and P2 disagree). Reported as such.
Known limits stated beforehand: ERSSTv6 is a statistical reconstruction (spatial smoothing may flatten extremes); a lower trend in warm cells
can also come from ocean dynamics (trade winds, ENSO-like pattern) without any thermostat, so step 1 can show the signature but cannot
separate the mechanisms; the 28-30 C bins are covered by few cells; no flux data are used in step 1.
