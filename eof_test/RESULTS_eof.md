# Results: dominant modes of the ocean SST field (run 2026-11-27; design in PREREG_eof.md)

ERSSTv6, 2-degree grid, annual means, 11,074 ocean cells without gaps, area-weighted. V1 = anomalies as they are, V2 = global mean removed.
Tables: `out/eof_modes.csv`, `out/eof_lagcc.csv`; PCs: `out/eof_pcs_*.csv`. Figures: `out/eof_fig.png` (patterns), `out/eof_pcs_fig.png` (time series), `out/eof_line_map.png`.

## E1 How many modes? (1900-2025)
Variance, PC1-5: V1 48.5 / 13.1 / 4.2 / 3.2 / 2.9 %; V2 23.5 / 9.9 / 7.2 / 5.5 / 5.1 %.
North's rule: in V1 only PC1 and PC2 are separable from their neighbours (PC3 and PC4 are not); in V2 only PC1 and PC5. **The data support two clear modes, then a continuum.**
- V1 PC1: warming, same sign nearly everywhere (strongest in the Arctic edge, North Atlantic, Southern Ocean).
- V1 PC2 (13%): El Nino / La Nina pattern, correlation with Nino3.4 |r| = 0.84 (sign of a PC is arbitrary). In V2 it is PC1 (23%).
- V1 PC3 (4%): North Atlantic warm against the Southern Ocean cool, a basin-scale Atlantic pattern; V1 PC4: a Pacific pattern (north-central Pacific plus eastern equatorial band) of the Pacific-decadal type. Not separable from each other, so read with care.

## E2 Which mode carries the 60-100 yr line?
1900-2025: **none significantly** (best AIC gains 0.6 to 2.6 against 95% null levels 2-4.5). One exception, V2 PC2 at 78 yr (p = 0.043), is one of 10 tests.
1850-2025: several PCs pass, but at **different periods** (V1: 75, 100, 82, 69 yr; V2: 100, 63, 85, 68 yr). No single period is shared by the modes.
Caveat: ERSSTv6 before about 1900 is a statistical reconstruction built from a small number of spatial modes, so an EOF analysis of those years partly reproduces the reconstruction method. I do not interpret the 1850 variant as independent evidence.

## E3 Where does the line live? (sinusoid of 76 yr, no PC was significant, so the default was used)
Strongest in the **North Atlantic (Gulf Stream and subpolar region), the central North Pacific at 30-45 N, the Mediterranean, and the Barents Sea / Arctic margin**, plus a band in the Southern Ocean. Weak in the tropics, especially the Indian Ocean and the western Pacific. This agrees with the latitude result (line only north of 20 N) and shows that the Pacific, not only the Atlantic, carries amplitude at that period. Amplitude is not significance: this map does not test the line against noise.

## E4/E5 Links and coupling (descriptive)
Lagged cross-correlations between PC1-4 stay below |r| = 0.4 at lags of 1-15 yr, except V1 PC2 with PC4 (r = 0.63 at lag 1 yr), i.e. ENSO and the Pacific decadal pattern, which are known to be related. No coherent long-lag relation between a "global" mode and an Atlantic mode appears. No significance is claimed for these numbers.

## Reading for Alf's plan (coupled oscillators)
1. The 126 yr of reliable data separate only about two modes (warming and ENSO). Fitting a coupled model with more modes, coupling strengths and damping will be under-determined: nothing here constrains it.
2. The 70-80 yr line does not appear as a mode of its own in the whole-ocean decomposition. It appears as regional amplitude in the northern basins. A basin-wise decomposition (North Atlantic, North Pacific separately) is the more promising next step, as Alf suggests.
3. A mode-based analysis of the tropical Pacific with the Bjerknes feedback needs subsurface data (thermocline depth, winds), which ERSST does not contain.
