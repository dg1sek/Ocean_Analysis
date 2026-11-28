# Results: tide gauges as an independent test of the 70-80 yr line (run 2026-11-28; rules in PREREG_sealevel.md incl. the amendment)

Data: PSMSL RLR annual means (Brest 1807-2023, Trieste 1875-2024, Genova 1884-2025, Marseille 1850-2024, Stockholm 1889-2024 as control; all with gaps). Series standardised. Quadratic trend absorbs land motion and the general sea-level rise. ERSST Mediterranean annual mean as comparison.

## First pass (pre-registered rule: longest continuous stretch; `out/*_stretch.csv`)
Stretches were short (Brest 1861-1943, n = 83; Genova 1928-1996; Marseille 1941-1986; Trieste 1901-2024; composite 1928-2024), the tests inconclusive. Kept for the record.

## P1b: is there a line at 60-100 yr? (gap-aware, all valid years; `out/sealevel_P1b_line_gapaware.csv`)
| series | n | best period (yr) | AIC gain | null 95% | p |
|---|---|---|---|---|---|
| Brest | 187 | 95 | 10.2 | 10.0 | 0.050 |
| Trieste | 131 | 100 | 1.5 | 6.4 | 0.20 |
| Genova | 109 | 100 | 2.5 | 8.0 | 0.22 |
| Marseille | 95 | 90 | 11.5 | 15.6 | 0.10 |
| Mediterranean composite | 121 | 100 | 5.6 | 8.9 | 0.09 |
| Stockholm (control) | 136 | 100 | -4.1 | 2.8 | 0.64 |
| **ERSST Mediterranean SST** (same method) | 176 | **74** | 109.3 | 41.4 | **< 0.003** |
Rule: p < 0.01 clear, p < 0.05 suggestive. **No sea-level series shows a clear line**; Brest is exactly at the suggestive limit (p = 0.050), the Mediterranean gauges are not significant. The best periods sit at the upper end of the scanned range (90-100 yr), not near 75.

## P2b: forecast without knowledge of the future (`out/sealevel_P2b_forecast_gapaware.csv`)
No series meets the rule (skill above the null 95% level at >= 2 usable cut-offs). Brest: -0.37 (1900), 0.31 (1950), -1.69 (1975); Mediterranean gauges and composite: only the 1975 cut-off usable, strongly negative (-0.6 to -18). Same caveat as before (quadratic baseline).

## P3b: phase at a fixed period of 75 yr (year of the last maximum, 95% interval; `out/sealevel_P3b_phase_gapaware.csv`)
| series | last maximum | 95% interval | amplitude / sd |
|---|---|---|---|
| ERSST Mediterranean SST | 2022.2 | 2019.3 - 2024.8 | 0.67 |
| Brest | 2013.9 | 2004.7 - 2021.6 | 0.21 |
| Trieste | 2021.4 | 2005.6 - 2039.6 | 0.11 |
| Mediterranean composite | 2024.3 | 2013.4 - 2032.1 | 0.18 |
| Genova | 1952.4 (+75 = 2027) | 1942.4 - 1961.0 | 0.16 |
| Marseille | 1951.2 (+75 = 2026) | 1941.2 - 1959.8 | 0.25 |
| Stockholm (control) | 2008.1 | 1974.6 - 2040.2 | 0.03 |
Modulo 75 yr the maxima of the four Mediterranean/Atlantic gauges (2014-2027) are in the same part of the cycle as the SST maximum (2022), with intervals of 10-30 yr; the control (Stockholm, amplitude 0.03) is not informative. The amplitudes are small (0.1-0.25 sd against 0.67 for SST), so this is a consistency of phase, not a detection. Period fixed at 75 yr beforehand, so the test cannot show that 75 is the right period.

## Reading
1. Tide gauges do not independently confirm the line: no clear line in any gauge, no forecast value.
2. The phases are compatible with the SST line but are weak signals; a signal this small cannot be told apart from noise with 100-200 years of data.
3. A fair statement: sea level neither supports nor contradicts the 75-yr line. Sea level also includes mass changes and land motion, which makes it a poor thermometer for this particular question.
