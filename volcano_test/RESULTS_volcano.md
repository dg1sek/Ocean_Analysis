# Volcano event study: results (2026-11-27)
Design and rules: `PREREG_volcano.md` (with the amendment). Script `20_volcano_events.R`, log `out/volcano_log.txt`, tables `out/volcano_*.csv`, figure `out/volcano_fig_composites.png`.
Six eruptions (Krakatau 1883, Santa Maria 1902, Novarupta 1912, Agung 1963, El Chichon 1982, Pinatubo 1991); monthly ERSSTv6 anomalies 1870-2019, slow variation and ENSO removed.

## A. Composite response (superposed epoch)
| series | peak cooling (K) | month of minimum | placebo median peak (K) | placebo p | recovery after min (months) |
|---|---|---|---|---|---|
| GLOBAL | -0.14 | 11 | -0.04 | 0.002 | 23 |
| WPWP | -0.30 | 27 | -0.13 | 0.015 | (noise) |
| MED | -0.47 | 26 | -0.24 | 0.011 | (noise) |
| GOM | -0.43 | 11 | -0.24 | 0.045 | (noise) |
| NATL_SUBTROP | -0.25 | 12 | -0.12 | 0.050 | (noise) |
Global: all events cool except Novarupta (+0.01): Krakatau -0.17, Santa Maria -0.24, Agung -0.24, El Chichon -0.15, Pinatubo -0.28 K. The global composite is clean.
The regional composites are dominated by noise (single regional events swing by +-1 K); the month of the minimum and the "recovery" are not meaningful there, and a residual ~12-month wave remains after the seasonal-cycle correction. Regional peaks are deeper than the placebo median but only marginally significant. Weak evidence only.

## B. One-box response model (exponential memory, input: AR6 monthly stratospheric AOD)
| series | tau (months), 95 % | gain (K per unit sAOD) | gain per W/m2 if ERF = -20 x sAOD | AIC gain vs no response |
|---|---|---|---|---|
| GLOBAL | 14 (7-28) | -2.25 | 0.11 | 37 |
| WPWP | 9 (4-18) | -3.71 | 0.19 | 26 |
| GOM | 34 (12-108) | -6.19 | 0.31 | 14 |
| NATL_SUBTROP | 70 (30-180, at grid edge) | -6.30 | 0.32 | 15 |
| MED | 86 (1-180, unconstrained) | -6.95 | 0.35 | 10 |
A two-timescale variant (fast + slow) is not preferred (AIC difference -2.1 to +2.0): the data do not identify a slow component.

## Reading
- The upper ocean (global mean, warm pool) answers a stratospheric-aerosol pulse with a memory of about one year (global 14 months, 7-28; warm pool 9 months, 4-18) and cools by about 0.14 K (six-event mean). It has recovered within about two years.
- The Atlantic regions give longer and larger responses (34-86 months) but their intervals are wide; with six noisy events I would not read a physical difference between regions into this.
- These events constrain the FAST response only. The slow, deep-ocean inertia (decades to centuries), which sets the long-term gain, is not visible in pulses that last two years.
- Consequence for a solar cycle with a one-year memory (rough calculation, not computed from data): for the 11-yr cycle the response would lag by about one year and keep about 83 % of its amplitude; for an 80-yr cycle it would follow with a lag of about one year and no damping. The size of the fast-response gain (about 0.1 K per W/m2) is far below the equilibrium sensitivity (of order 0.8 K per W/m2), so the amplitude of an 80-yr response is set by the slow part, which this test cannot measure.

## Limits
Six events; El Chichon and Pinatubo coincide with El Nino years (ENSO removal imperfect); Krakatau in the partly filled ERSSTv6 period; pre-1979 aerosol amounts uncertain; the conversion -20 W/m2 per unit sAOD is quoted from memory (AR6 Chapter 7) and not verified.

## C. The conversion factor "-20 W/m2 per unit sAOD" (`31_ar6_factor_check.R`)
- In the AR6 files, volcanic ERF is EXACTLY -20 x (sAOD - baseline) (regression of the annual AR6 ERF on the annual AR6 sAOD, 170 years: slope -20.000, R2 = 1.000000). The factor is a construction rule of AR6, not a measurement.
- Literature (web search, abstracts and summaries only): AR5 used -25 W/m2 per unit sAOD (from GISS ModelE Pinatubo simulations); Schmidt et al. (2018) report the ERF-AOD relation 13-21 % weaker than AR5 for large eruptions (rapid adjustments), which brings -25 down to about -20; Marshall et al. (2020) note that spatial redistribution of the aerosol can raise the forcing per global mean AOD; other work finds large variation of the forcing efficiency between eruptions (source parameters). So -20 is a reasonable central value with a spread of roughly -13 to -25 and a spread between single eruptions.
- Consequence: our fitted gains are in K per unit sAOD and do not depend on the factor. Only the conversion to K per W/m2 does: global 0.09-0.17 (factor -25 ... -13), 0.11 at -20; warm pool 0.15-0.29; Mediterranean 0.28-0.53 (table `out/volcano_gain_vs_factor.csv`). The conclusion "fast response of order 0.1 K per W/m2, far below the equilibrium sensitivity" holds over the whole range.

## D. Can the Sun produce an 80-year line? (`30_solar_longperiod.R`)
Amplitude of the best sinusoid between 60 and 120 yr:
| series | period (yr) | amplitude |
|---|---|---|
| NOAA TSI 1610-2025 | 110 | 0.11 W/m2 TSI |
| NOAA TSI 1749-2025 | 90 | 0.13 W/m2 TSI |
| NOAA TSI 1850-2025 | 94 | 0.11 W/m2 TSI |
| AR6 solar ERF 1750-2019 | 98 | 0.017 W/m2 ERF |
TSI amplitude x 0.175 (albedo 0.3, geometry) = 0.02 W/m2 ERF, consistent with the AR6 series.
SST amplitude this could drive: 0.003 K (fast response 0.11 K per W/m2), 0.01 K (transient 0.46), 0.02 K (equilibrium 0.76). Observed 74-80 yr line: 0.23-0.30 K (Mediterranean, Gulf of Mexico).
Energy view (heat capacity, following Alf's remark of 23 Nov): a forcing of amplitude F acting on a mixed layer of depth h without any radiative loss gives T = F/(rho cp h w). For P = 80 yr and h = 50-200 m, the solar ERF of 0.017 W/m2 gives 0.008-0.033 K, 9-36 times less than 0.3 K. Heating 75 m of the global ocean by 0.3 K needs 3.3e22 J, the solar line delivers at most 0.8e22 J in a half cycle even if everything is absorbed and nothing lost.
Limits: NOAA and AR6 are low-variability reconstructions; reconstructions with a much larger long-term solar variability exist (to be checked, from memory e.g. Shapiro et al. 2011) and would shrink the gap; the energy argument is for a global layer, regional amplitudes can also be fed by ocean heat transport (a self-sustained or dynamically driven oscillation needs no external clock). The test shows that the Sun in these reconstructions is too weak by an order of magnitude, not that no solar influence exists.
