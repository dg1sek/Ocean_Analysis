# Results: shell chronology as a test of the 65-80 yr line (run 2026-11-28; rules in PREREG_bivalve.md)

Data: Butler et al. 2013, North Icelandic Shelf (Grimsey, 66.5 N 18.2 W), 29 shells, 5,705 annual increments, AD 649-2005; own RCS chronology (not pre-whitened). At least 3 shells from AD 1175. Figures: `out/bivalve_fig.png`; tables: `out/bivalve_*.csv`.

## B1 Best period between 40 and 120 yr
| series | years | n | best period (yr) | AIC gain | null 95% | p |
|---|---|---|---|---|---|---|
| own RCS, main | 1175-2005 | 831 | **54** | 6.1 | 7.2 | 0.073 |
| own RCS, first half | 1175-1590 | 416 | 53 | 3.7 | 6.0 | 0.14 |
| own RCS, second half | 1590-2005 | 416 | 40 | 5.2 | 6.6 | 0.077 |
| own RCS, full | 649-2005 | 1357 | 115 | 3.2 | 7.6 | 0.35 |
| published NE.ars (pre-whitened) | 649-2005 | 1357 | 40 | 4.4 | 7.7 | 0.20 |
| published RCS.ars (pre-whitened) | 649-2005 | 1357 | 46 | 0.8 | 7.6 | 0.76 |
In the 65-80 yr band the AIC gain of the main series is at best -2.9 (68 yr; at 75 yr: -5), i.e. the line is not there.
Rule: the 65-80 yr line is "present" only with a best period of 60-90 yr and p < 0.05. **Result: not shown.** No version has p < 0.05 at all; the best period (54 yr in two of three own windows) is shorter than the SST line and not significant.

## B3 Link to SST
Shell index against the ERSSTv6 annual mean in a box around the site (62-70 N, 28-10 W), 1854-2005: r = 0.13 (effective n = 96, p = 0.21). The authors report r about 0.2 for regional SST as well. The chronology is a weak thermometer.

## Reading
1. This 1,357-year shell chronology does not show a 65-80 yr line. It shows no significant periodicity in the 40-120 yr range.
2. That is not evidence against the line in the Mediterranean: the site is near the polar front, the shell growth responds only weakly to SST, and one site cannot represent the whole North Atlantic. The literature claim of a 65-80 yr period comes from other shells (southern Barents Sea, Mette et al. 2021); I did not find a public data file for that record and have not tested it.
3. Consequently the pre-industrial test remains open. It would need a better proxy of Atlantic SST (for example the Barents Sea shells, or an AMO reconstruction).
4. Limits: RCS leaks low-frequency signal into or out of the chronology; low sample depth before AD 1175; a single site.

---
# Part 2: Barents Sea shells (Mette et al. 2021) and a tree-ring AMO reconstruction (Gray et al. 2004) -- run 2026-11-28
Tables: `out/proxies_C1_line.csv`, `out/proxies_C2_phase.csv`; figure `out/proxies_fig.png`; script `96_proxies.R`.

## C1 Best period between 40 and 120 yr (ML sinusoid, AR(1) null, 300 surrogates)
| series | years | n | best period (yr) | AIC gain | null 95% | p | best AIC gain inside 65-80 yr |
|---|---|---|---|---|---|---|---|
| Barents growth, full | 1426-2012 | 587 | **55** | 10.8 | 6.2 | 0.010 | 1.2 |
| Barents growth, first half | 1426-1719 | 294 | 56 | 3.1 | 5.3 | 0.16 | 0.0 |
| Barents growth, second half | 1720-2012 | 293 | 120 | 6.5 | 4.4 | 0.010 | -0.8 |
| Barents d18O, full | 1539-2014 | 476 | **101** | 7.3 | 5.6 | 0.030 | -2.9 |
| Barents d18O, first half | 1539-1776 | 238 | 114 | 2.1 | 4.6 | 0.16 | -3.1 |
| Barents d18O, second half | 1777-2014 | 238 | 101 | 8.3 | 5.2 | 0.007 | -2.2 |
| Gray AMO (annual), full | 1567-1990 | 424 | **94** | 7.2 | 5.1 | 0.017 | -0.2 |
| Gray, pre-instrumental | 1567-1855 | 289 | 99 | 4.2 | 5.1 | 0.087 | -3.5 |
| Gray, instrumental | 1856-1990 | 135 | 70 | 5.0 | 4.1 | 0.027 | 5.0 |
Pre-registered rule (best period 60-90 yr with p < 0.05 in the full series and in a subperiod): **not met in any proxy; the 65-80 yr line is "not shown" before the instrumental era.** Its literature claim (65-80 yr, Mette et al. 2021, from their own spectral method) is not reproduced by this scan: the gain inside 65-80 yr is at most 1.2 (Barents growth) and negative for the isotopes.
The only place where a period of 65-80 yr appears is Gray's instrumental part (70 yr), which is the calibration period of that reconstruction and therefore not independent.

## C2 Phase at a fixed 75 yr (only descriptive)
Amplitudes at 75 yr are tiny in the shell series (0.05-0.08 sd against 0.67 for the Mediterranean SST), so their phases (offset to the Mediterranean maximum: Barents growth -24 yr, d18O +29 yr, intervals of 30-60 yr) say nothing. Gray: offset -0.2 yr pre-instrumental (interval 1987-2055, wide) and -2.5 yr full; the full-record amplitude is 0.33 sd.

## Exploratory observation, not part of the rules and made after seeing the results
Three records give a significant best period near **100 yr** that is not at the edge of the scanned range: Barents d18O (101 yr, p = 0.03; 101 yr in the second half, p = 0.007), Gray AMO (94 yr, p = 0.017; 99 yr pre-instrumental, p = 0.087 only). The Barents growth chronology gives 55 yr (p = 0.01) in the full series and 120 yr in its second half. The sea-level gauges also had best periods of 90-100 yr but at the edge of their scan (60-100), so they do not count here. The common feature of the proxies is a centennial rather than a 70-80 yr period. Warnings: the chronology standardisation of tree rings and shells limits variance at periods beyond roughly 100-150 yr and can shape a peak near its edge; 4-5 cycles in 400-500 yr; no correction for the number of tests. This needs its own pre-registered test before anyone reads anything into it.

## Reading
1. Before 1850 the proxies available to us do not show a 65-80 yr periodicity. That weakens the idea of a stable natural cycle of that length; it does not rule it out (proxies are weak thermometers, standardisation limits low frequencies).
2. The Mediterranean SST line (74 yr, 1850-2025) is so far seen only in the instrumental period.
