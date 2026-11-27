# Thermostat test, steps 2 and 3: results (2026-11-27)
Design and decision rules: `PREREG_step1.md`, `PREREG_step2_3.md` (fixed before the numbers; one amendment recorded there: OISST v2 ends Jan 2023, common period 1982-2022).

## Step 2: three SST products (`11_step2_products.R`, `out/step2_*.csv`, `out/step2_log.txt`)
Common period 1982-2022 (P5), tropical cells |lat| <= 30 with complete annual series.
| Product | cells | A2: slope of trend vs starting SST (K/dec per K, cells >= 27 C) | B1: trend(q99) - trend(median) (K/dec) |
|---|---|---|---|
| ERSSTv6 (2 deg, filled) | 4231 | +0.019 (0.002 to 0.030) | +0.044 (0.021 to 0.057) |
| OISST v2 (1 deg, satellite) | 16351 | +0.021 (0.002 to 0.031) | +0.049 (0.019 to 0.069) |
| HadSST 4.2 actuals (5 deg, unfilled) | 313 | +0.006 (-0.019 to 0.027) | +0.028 (-0.033 to 0.074) |
All point estimates have the same sign (warm end warms as fast as or faster than the rest). ERSSTv6 and OISST exclude zero on the positive side, HadSST4 is too coarse
(313 cells) to say anything either way. No product shows the thermostat signature. Area above 30 C (monthly, tropics) rose 0.6 % -> 4.4 % (ERSSTv6), 1.3 % -> 4.2 % (OISST), 1.0 % -> 6.3 % (HadSST4) over 1982-2022.
ERSSTv6 P1/P2/P4 reproduce step 1 (same cells now also restricted to complete cells for the quantile test).

## Step 3a: SST relative to the tropical mean
No separate test (argued beforehand): a relative-SST threshold is observationally equivalent to uniform warming; non-uniform warming is tested by A2, which is invariant to the shift. Result of A2/B1 above stands for it.

## Step 3b: surface-flux feedback (`12_step3_flux.R`, `out/step3_flux_slopes.csv`, `out/step3_log.txt`)
NCEP/NCAR Reanalysis 1 fluxes (a model product), ERSSTv6 SST, 1982-2025 monthly anomalies, detrended, 1850/1480/1338 cells (cool < 26 C / mid / warm >= 28 C). Slope of the heat flux into the ocean on SST anomaly, W/m2 per K, block-bootstrap 95 % intervals:
| flux | warm | cool | warm - cool |
|---|---|---|---|
| net Q | -1.2 (-3.2 to 0.3) | -1.2 (-3.1 to 1.0) | +0.01 (-2.5 to 1.9) |
| shortwave (net) | -5.5 (-6.1 to -5.0) | +0.2 (-0.4 to 1.0) | -5.7 (-6.8 to -4.9) |
| latent (as heating) | +2.5 (0.8 to 3.8) | -1.9 (-3.0 to -0.4) | +4.4 (2.5 to 5.6) |
| longwave (net) | +1.7 (1.3 to 2.0) | -0.1 (-0.3 to 0.1) | +1.8 (1.3 to 2.2) |
Reading: in warm cells a warmer surface goes with 5.5 W/m2 less absorbed sunlight per K (cloud shading, the Ramanathan-Collins mechanism), but this is offset by more longwave
trapping (+1.7) and by less evaporation (+2.5), so the NET damping is the same in warm and cool cells (difference 0.0, interval -2.5 to +1.9).
Pre-registered rule: thermostat supported if the net-flux difference is negative with the interval excluding zero. It is not. The shortwave component alone does show the shading effect.
Caveats: Reanalysis fluxes depend on the model's cloud and surface scheme. The slope on monthly anomalies mixes feedback and forcing (an El Nino changes winds and clouds and
the SST together), so it is not a clean feedback estimate; a lagged or tendency-based version, and CERES/OAFlux, would be the next refinement.

## Overall
Steps 1-3: no sign, in three SST products and in reanalysis fluxes, that the warmest tropical water warms more slowly or damps its anomalies more strongly than cooler water.
One component of the thermostat (cloud shading in warm cells) is visible, and cancelled in the net. This does not prove the mechanism absent at a different scale or in observations of higher quality.
