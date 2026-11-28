# Results: simulation with coupled stochastic oscillators (run 2026-11-28; rules in PREREG_sim.md, `97_simulation.R`, `out/sim_results.csv`, `out/sim_real_med.csv`, `out/sim_fig.png`)

400 synthetic series per scenario (176 annual values with the Mediterranean quadratic trend, three noise-driven AR(2) oscillators with Q = 2-10, AR(1) background, residual sd 0.287 K as in the real series). The same detection procedure as for the real data; 100 AR(1) or AR(2) surrogates per series.

| scenario | what is true | flagged "line 60-90 yr, p < 0.05" (AR(1) null) | (AR(2) null) | forecast rule passes | AIC gain >= 13.5 (the real MED value) |
|---|---|---|---|---|---|
| S0 | pure AR(1) noise | 2.5 % | 3.2 % | 0.2 % | 0.0 % |
| S1a | oscillators, none within 55-110 yr | 4.8 % | 2.8 % | 3.0 % | 0.2 % |
| S1b | oscillators, periods random (about 59 % of the series have one within 55-110 yr) | 18.0 % | 13.2 % | 8.5 % | 3.0 % |
| S2 | as S1b, one forced to 65-85 yr | 43.5 % | 31.2 % | 13.8 % | 4.0 % |

Real Mediterranean series (1850-2025): best period 74 yr, AIC gain 13.5, p < 1/300 with both AR(1) and AR(2) null.

## Reading (rules fixed in advance: false-positive rate above 10 % would call the p-values optimistic)
1. Where no oscillator lies near 75 yr (S1a) the procedure reports a line in 5 % (AR(1)) or 3 % (AR(2)) of the series: the p-values are not optimistic in that case. The check S0 is also correct.
2. In S1b the flagged fraction is 18 % / 13 %, but most of these series do contain an oscillator within the band by construction; this is not a false positive in the strict sense. It shows what Alf points out: in a system of coupled oscillators a resonance near 75 yr occurs by chance in more than half of the systems (about 59 % under these assumptions) and then appears as a "line".
3. Power is low: with a 65-85 yr oscillator really present (S2) the procedure flags it in only 31-44 % of the series.
4. The real AIC gain (13.5) is reached in 3-4 % of S1b/S2 and in 0.2 % of S1a: it is large compared with what such oscillators produce, **with these assumed strengths and damping**. A stronger or less damped oscillator in the real system cannot be excluded.
5. The forecast rule passes in 3 % (S1a), 8.5 % (S1b), 13.8 % (S2) of the series. The Mediterranean passed it (quadratic baseline), so the forecast evidence corresponds to a chance of roughly 3-9 % under "no stable line" and is weak.
6. What the simulation cannot decide: whether the Mediterranean line is a stable sinusoid or a quasi-periodic noise-driven resonance. Both give the same picture in 176 years.
