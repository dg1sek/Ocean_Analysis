# Tropical thermostat, steps 2 and 3: pre-registered design (written 2026-11-27 before any result was computed)

## Step 2: independent SST products
Same tests and the same decision rule as PREREG_step1.md (A2 slope of trend vs starting SST for cells >= 27 C; B1 trend(q99) - trend(median)),
applied unchanged to three products: ERSSTv6 (2 deg, filled), OISST v2 (1 deg, satellite era, from 1982), HadSST 4.2 "actuals" (5 deg, NOT filled).
Common period P4 = 1982-2025 (baseline 1982-2001) for all three. ERSSTv6 and HadSST4 also P1 (1950-2025) and P2 (1980-2025) where cells are complete.
Cells: |lat| <= 30; only cells with a complete annual series (annual mean needs 12 valid months; every year in the period valid) enter A and B, so that coverage changes cannot fake a trend.
Number of cells is reported for each product.
Reading rule: agreement = same sign of A2 and B1 in all three products on P4. If the products disagree, the verdict for ERSSTv6 is
downgraded to "product dependent". The HadSST4 grid is coarse (5 deg), so its warmest bin may hold very few cells; bins with < 5 cells are not drawn.

## Step 3a: SST relative to the tropical mean
Argued before computing: a relative-SST threshold (convection depends on SST minus tropical mean) is observationally equivalent to uniform warming,
because subtracting the tropical mean from each cell removes exactly the part that a relative thermostat would leave unchanged. The test of
non-uniform warming (does the warm end warm slower than the tropical mean?) is test A, which is invariant to the shift (trend minus mean trend). So 3a
adds no new test; it is reported as: relative upper tail trend = B1 value (trend(q99) - trend(mean)), computed from the same tables.

## Step 3b: flux feedback (cloud, evaporation)
Not run in this round if no usable open data can be fetched without a login. Candidate sources (to be checked): NOAA PSL reanalysis fluxes (latent heat, net shortwave),
CERES EBAF (needs NASA Earthdata login), ERA5 (needs CDS key). Test, fixed now: monthly anomalies of the net surface heat flux into the ocean and of latent heat flux
against SST anomaly in the warm pool (cells with climatological SST >= 28 C), slope in W m-2 per K with block-bootstrap interval, versus cells < 26 C.
A thermostat predicts a more negative slope in warm cells. Reanalysis fluxes are model products; this caveat is part of the result.

## Amendment (2026-11-27, after the first run of ERSSTv6 and HadSST4, before any OISST result existed)
OISST v2 (1 deg monthly file at NOAA PSL) ends in January 2023, so P4 (1982-2025) cannot be run on it. A common period P5 = 1982-2022 (baseline 1982-2001)
is added and is the period used for the cross-product agreement rule; the rule itself is unchanged. The ERSSTv6 and HadSST4 results for P1, P2, P4 stay as first computed.
