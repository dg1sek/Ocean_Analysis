# Pre-registration: shell chronology (Arctica islandica, North Icelandic Shelf) as a test of the 65-80 yr line -- written 2026-11-28 before any result

Data: Butler et al. (2013), Palaeogeogr. Palaeoclimatol. Palaeoecol. 373, 141-151; NOAA NCEI "North Icelandic Shelf Arctica islandica Shell Measurement 1350 Year Record" (butler2013-raw.txt: annual increment widths in micrometres of 29 absolutely dated shells, 649-2005 AD; butler2013-chron.txt: the authors' published chronologies NE.ars and RCS.ars, which are autoregressively modelled (pre-whitened) versions).
Why this series: annual, marine, the longest of its kind (1,357 years), near the polar front where the Atlantic multidecadal signal is expected. The authors themselves report only weak correlations with regional SST (r about 0.2), so it is an imperfect thermometer.

Own chronology (the published versions are pre-whitened, which suppresses low frequencies, so they are not used for the main test): regional curve standardisation. Mean increment by biological age over all shells, smoothed with an 88-yr spline (as the authors), each increment divided by that curve, then the arithmetic mean over the shells present in each calendar year. Main period: AD 1175-2005 (at least 3 shells); check: AD 649-2005 (1-2 shells before 1175).

Questions and rules
 B1 Line: ML sinusoid scan (alf_tools::sin_scan, AR(1) errors, linear trend) over 40-120 yr in 1-yr steps; best period, AIC gain, p against 300 AR(1) surrogates. The 65-80 yr line is "present" if the best period lies within 60-90 yr with p < 0.05 for 1175-2005 AND the same band holds the best period, or p < 0.05, in at least one of two non-overlapping halves (1175-1590, 1590-2005). Otherwise "not shown".
 B2 Same scan on the published NE.ars and RCS.ars (comparison only, no rule).
 B3 Link to SST: correlation of the chronology with the ERSSTv6 annual mean in a box around Grimsey (62-70 N, 28-10 W, available years, complete months only) over 1854-2005, with the effective sample size corrected for AR(1) of both series. Reported only.
No tests are added after seeing the results.
Limits fixed in advance: RCS assumes one growth curve for all shells and can leak low-frequency signal into or out of the chronology; low sample depth; growth responds to food supply and water mass, not only temperature.

## Part 2 (2026-11-28, written before this part was run): Barents Sea shells and a tree-ring AMO reconstruction
Added after part 1 gave "not shown" for the North Icelandic shells.
Data: (a) Mette et al. 2021, J. Geophys. Res. Oceans, 10.1029/2020JC017074 (NOAA NCEI mette2021.txt): southern Barents Sea Arctica islandica, master shell growth chronology (MShellCrn, 1449-2012; standard ARSTAN chronology, not pre-whitened) and shell d18O (1539-2014); the reconstructed SST in the file is a transfer function of the d18O and not analysed separately. (b) Gray et al. 2004, GRL 31, L12205 (NOAA NCEI amo-gray2004): tree-ring reconstruction of North Atlantic SST anomaly, annual column `ssta`, 1567-1990 (the 10-yr smoothed `amo` column is not used because smoothing creates periodicity). Gray's reconstruction is calibrated on instrumental SST after 1856, so the line may be built in there; therefore the pre-instrumental part 1567-1855 is tested separately and is the one that counts.
Same method as part 1 (alf_tools sin_scan, AR(1) errors, linear trend, 40-120 yr, 300 AR(1) surrogates).
 C1 Rule: the line is "present" in a proxy if the best period lies within 60-90 yr with p < 0.05 for the full series AND the same holds in at least one subperiod (Barents: either half; Gray: 1567-1855). Otherwise "not shown".
 C2 Phase at a fixed 75-yr period (year of the last maximum with a 95% moving-block-bootstrap interval; quadratic trend removed), compared modulo 75 yr with the ERSST Mediterranean maximum (2022). Descriptive only.
No tests are added after seeing the results.
