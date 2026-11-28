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
