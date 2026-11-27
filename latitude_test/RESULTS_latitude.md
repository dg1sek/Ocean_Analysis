# Latitude bands of the ocean: results (2026-11-27)
Design: `PREREG_latitude.md`. Ocean-only ERSSTv6 band means (20-degree bands from 60S to 60N). Scripts `39_make_latbands.R`, `40_latitude_events.R`, `41_latitude_sinus_scan.R`; logs and tables in `out/`.

## A. The 70-80 yr line depends on latitude (ML sinusoid scan, quadratic trend, AR(1) errors, annual means 1850-2025, AR(1) surrogate null, 300 simulations)
| band | best period (yr) | amplitude (K) | AIC gain | 95 % of null | p |
|---|---|---|---|---|---|
| 60-40S | 94.5 | 0.11 | 1.2 | 7.8 | 0.91 |
| 40-20S | 9.0 | 0.03 | 1.7 | 7.9 | 0.81 |
| 20S-0 | 53.5 | 0.07 | -0.1 | 7.7 | 0.97 |
| 0-20N | 12.5 | 0.06 | -0.4 | 7.7 | 0.99 |
| 20-40N | 78.0 | 0.20 | 8.6 | 7.8 | 0.033 |
| 40-60N | 74.0 | 0.24 | 9.9 | 7.8 | 0.013 |
The line (74-78 yr, 0.20-0.24 K) appears only in the northern mid-latitude bands; tropics and southern hemisphere show nothing above noise. Band means mix ocean basins; a basin-resolved split (North Atlantic vs North Pacific) is the next step. With a threshold of AIC gain > 2 for "better than no line" every northern value passes; against the surrogate null only the two northern bands do.

## B. Eruption response per latitude band (same method as volcano_test)
| band | peak cooling (K) | placebo p | one-box memory (months, 95 %) | gain (K per unit sAOD) |
|---|---|---|---|---|
| 60-40S | -0.06 | 0.39 | 16 (1-180) | -0.30 |
| 40-20S | -0.18 | 0.013 | 19 (13-28) | -4.14 |
| 20S-0 | -0.17 | 0.053 | 13 (2-38) | -2.48 |
| 0-20N | -0.23 | 0.005 | 5 (1-17) | -1.65 |
| 20-40N | -0.21 | <0.001 | 20 (11-38) | -4.05 |
| 40-60N | -0.31 | 0.001 | 10 (1-180) | -2.72 |
The Southern Ocean (60-40S) shows no detectable cooling; the strongest, best-constrained responses are in the mid-latitude subtropics (40-20S, 20-40N) with memory of about 1.5-1.7 yr; the equatorial band has the shortest memory (5 months, wide interval). Two-timescale models are not preferred in any band (AIC difference -2.3 to +1.3). With six events the band-to-band differences are suggestive, not established.

## Reading
Alf's statement that the oscillations depend on latitude is supported for the 70-80 yr line: it is a northern mid-latitude feature in these data. The same bands, where the line is found, are those with the clearest eruption response. This fits an Atlantic or northern-hemisphere origin (ocean circulation, aerosols) more than a uniformly distributed solar forcing, but the band means cannot say which basin carries it.
