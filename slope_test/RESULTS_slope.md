# Results: the slope of the SST curve (run 2026-11-27; design in PREREG_slope.md)

ERSSTv6 annual means, model = quadratic trend + one sinusoid (60-100 yr, ML, AR(1)), 500 block-bootstrap replicates.
Full table: `out/slope_results.csv`. Sinusoid phase: `out/sinus_phase.csv` (51_sinus_phase.R, an extra check made after seeing the first result).

## Q1 Does the trend slope fall? (acceleration c, K per century^2)
| series | 1850-2025: c [95%] | 1900-2025: c [95%] | trend slope 1996-2025 (K/decade) |
|---|---|---|---|
| GLOBAL | 0.67 [0.58, 0.77] | 0.67 [0.44, 1.06] | 0.14 |
| MED | 0.75 [0.52, 0.98] | 1.83 [1.37, 2.43] | 0.16 - 0.27 |
| GOM | 0.65 [0.49, 0.81] | 1.41 [0.65, 1.88] | 0.11 - 0.21 |
| WPWP | 0.98 [0.85, 1.09] | 1.05 [0.66, 1.51] | 0.17 - 0.18 |

c > 0 in all four series and both periods: the rule gives **"accelerates", no deceleration anywhere**.
Caveat: the bootstrap intervals are probably too narrow (AR(1) phi about 0.7, 5-yr blocks); the sign is robust, the width is not.

## Q2 30-yr window slopes (K/decade)
Raw series, 1936-65 / 1966-95 / 1996-2025: GLOBAL -0.01 / 0.11 / 0.18; MED 0.02 / 0.02 / 0.46; GOM -0.02 / 0.03 / 0.32; WPWP -0.04 / 0.10 / 0.21.
Difference 1996-2025 minus 1966-95 is positive in all, with the interval above zero for GLOBAL, MED, GOM (WPWP undetermined).
After subtracting the fitted sinusoid: undetermined in all four for 1850-2025; "accelerates" in all four for 1900-2025.
The window slope of MED (0.46) is far above its trend slope (0.16): the sinusoid adds to the slope at present.

## Q3 Size of the sinusoid's own slope
Amplitude of the slope = 2 pi A / P: MED 0.26, GOM 0.18 K/decade (GLOBAL 0.05-0.06, WPWP 0.07-0.08, those lines are weak, dAIC < 4).
For MED and GOM this is **as large as the trend slope itself** (0.11-0.27 K/decade). So the 70-80 yr line decides whether
the regional curve looks steep or flat in a given 30-yr window; it is not a small correction.

## Extra check (made after seeing Q1-Q3): where is the sinusoid now?
| series | fit period | last maximum | sinusoid slope 2025 (K/decade) |
|---|---|---|---|
| MED, 1850-2025 | 74 yr | 2022 | -0.07 |
| MED, 1900-2025 | 85 yr | 2031 | +0.12 |
| GOM, 1850-2025 | 80 yr | 2029 | +0.06 |
| GOM, 1900-2025 | 100 yr | 2048 | +0.19 |
The phase moves by 10-20 yr between the two fit periods, so "now" is uncertain.
If the line is real and the 1850 fit is right, MED is at its maximum now and its slope should fall by about 0.26 K/decade
amplitude over the next decades, which would cancel most of the trend slope (0.16): a plateau or slight decline in MED between about 2025 and 2060.
This is a **testable prediction, not a result**; the 1900 fit says the maximum is still 10 years away.

## Reading
1. No sign of a slowing trend in global, MED, GOM or WPWP: the opposite (acceleration), as expected under a rising forcing.
2. The slope of the regional series is dominated by the sinusoid, not by the trend. An apparent "flattening" in MED would be the line, not a ceiling.
3. Limits: one sinusoid plus quadratic is Alf's model, not a test of it; ERSST is statistically filled before about 1900; not compared with forcing (W2 does that globally).

**Correction 2026-11-28:** `51_sinus_phase.R` and the figure script `52_fig_slope.R` had used a partial year 2026 (January-October) as an extra annual value; the table above is the corrected one (complete years only). Q1-Q3 (`50_slope.R`) were never affected because that script cuts at 2025. The qualitative reading is unchanged: for the 1850 fit the MED sinusoid is at its maximum about now (2022), for the 1900 fit about ten years later (2031).
