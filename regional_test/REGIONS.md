# Candidate regions for the regional analysis

Selected before any solar or 20-yr signal was looked at (criteria in `10_region_screening.R`). Region boxes: `R/regional_lib.R`.

| Rank | Region | Band | Centre lat | Mean insolation W/m2 | Seasonal range W/m2 | Independent proxy (>=270 yr) | sd annual (K) | early/late sd | r with Nino3.4 | Score (max 10) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **GOM** Gulf of Mexico | subtropical 15-30 | 24.0 | 381 | 196 | Dry Tortugas (Pulaski Reef) | 0.25 | 1.18 | 0.08 | 8 |
| 2 | **WPWP** Western Pacific warm pool | tropical 0-15 | 0.0 | 415 | 53 | Warm pool index (Freund 2019) | 0.16 | 1.04 | 0.27 | 8 |
| 3 | **NWPAC** North-west Pacific (Kuroshio, Ryukyu) | subtropical 15-30 | 27.5 | 371 | 227 | Ryukyu (Kikai Island) | 0.26 | 0.99 | 0 | 7 |
| 4 | **SCS** South China Sea | tropical 0-15 | 13.5 | 404 | 106 | Hon Tre Island (Vietnam) | 0.24 | 1.43 | 0.23 | 7 |
| 5 | **SWPAC** South-west Pacific / SPCZ (Fiji, New Caledonia, Rarotonga) | subtropical 15-30 | -17.5 | 399 | 188 | Fiji Savusavu Bay AB; New Caledonia Amedee; Rarotonga | 0.18 | 1.07 | 0.42 | 7 |
| 6 | **TIND** Tropical Indian Ocean | tropical 0-15 | -5.0 | 414 | 83 | none downloaded | 0.18 | 1.07 | 0.48 | 6 |
| 7 | **MED** Mediterranean Sea | mid/high 30-65 | 38.0 | 334 | 312 | none downloaded | 0.29 | 0.84 | 0.13 | 5 |
| 8 | **NATL_SUBTROP** Subtropical North Atlantic (Sargasso) | mid/high 30-65 | 32.5 | 354 | 268 | none downloaded | 0.24 | 0.72 | 0.02 | 5 |
| 9 | **NINO34** Nino3.4 (central equatorial Pacific) | tropical 0-15 | 0.0 | 415 | 53 | Nino3.4 reconstruction (Wilson 2010 / Freund 2019) | 0.58 | 0.98 | - | 5 |
| 10 | **REDSEA** Red Sea | subtropical 15-30 | 20.5 | 390 | 165 | none downloaded | 0.27 | 1.26 | 0.04 | 5 |
| 11 | **NATL_SUBPOLAR** Subpolar North Atlantic | mid/high 30-65 | 57.5 | 245 | 439 | none downloaded | 0.34 | 0.76 | 0.08 | 3 |
| 12 | **NPAC** North Pacific mid-latitudes | mid/high 30-65 | 37.5 | 335 | 308 | none downloaded | 0.39 | 0.73 | 0.45 | 3 |

**Score** = proxy (0-3) + early-period reliability (0-2) + low ENSO contamination (0-2) + low noise (0-2) + named by Alf (0-1).
Insolation is reported, not scored. NINO34 is scored 0 on ENSO by construction. Proxy records: `data/proxy_sites.csv`
(Rarotonga, 269 yr, falls just under the 270-yr limit used for the count; Red Sea, tropical Indian Ocean and subtropical
North Atlantic corals exist at NOAA NCEI but have not been downloaded).

## Choice for the full analysis: best region in each latitude band

| Band | Region | Score | Why |
|---|---|---|---|
| tropical 0-15 | **WPWP** western Pacific warm pool | 8 | highest annual insolation (415 W/m2), low noise (0.16 K), one independent reconstruction (Freund 2019 warm-pool index), little ENSO leakage (r = 0.27) |
| subtropical 15-30 | **GOM** Gulf of Mexico | 8 | named by Alf, independent coral (Dry Tortugas, 1734-2008), almost no ENSO leakage (r = 0.08) |
| mid-latitude 30-65 | **MED** Mediterranean | 5 (tie with the subtropical North Atlantic) | Alf's own region (he can compare with his own result); largest seasonal insolation range (312 W/m2) of the three; no independent long proxy downloaded |

Runners-up with a coral record (follow-up candidates): NWPAC (Ryukyu, score 7), SCS (Hon Tre, 7), SWPAC (Fiji, New Caledonia, Rarotonga, 7).
Selection did not use any 20-yr, 80-yr or solar signal of the regions.
