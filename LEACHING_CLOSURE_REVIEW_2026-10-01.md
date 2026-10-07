# Dissolved cation export and the background salt closure, 1 October 2026

Package n. Two consistency corrections in how the 0-30 cm box exports what
its drainage water carries, one re-choice of an existing parameter, and the
scoring moved to the month the soil was sampled. No rate constant,
mineralogy, surface area, extract operator or carbonate sink changes.

## Why

August is when the soil cores were taken (Beerling et al. 2024, Methods;
SI Fig. S7). Against the 0-30 cm depth-weighted observation, package m's
basalt-minus-control pH effect is right to within 0.13 through 2018 and
then 2.4 times too strong in 2019 and 4.1 times in 2020. The alkalinity
budget of m's runs 10 and 11 (`runtime/Base_1D/budget.m`) shows why:

| run 11, August | 2016 | 2017 | 2018 | 2019 | 2020 |
|---|---|---|---|---|---|
| base saturation, % (obs 0-10 / 10-30) | 78.3 (80.7 / 83.3) | 80.0 (83.2 / 79.6) | 85.7 (89.6 / 81.1) | 96.0 (87.8 / 80.6) | 99.3 (83.4 / 86.5) |
| f_H, % | 21.6 | 19.8 | 14.3 | 3.96 | 0.69 |
| pore-water alkalinity, meq L-1 | 0.07 | 0.02 | 0.36 | 0.44 | 3.2 |
| alkalinity released by rock, cumulative eq m-2 | 0 | 1.8 | 5.6 | 10.2 | 12.6 |
| Ca+Mg leached, cumulative eq m-2 (control) | 0.16 (0.16) | 0.53 (0.53) | 0.93 (0.92) | 1.53 (1.50) | 1.95 (1.83) |

The exchanger saturated and the runaway coincides with it. The basalt run
exported 0.12 eq m-2 more Ca+Mg and 0.27 mol more C than the control over
four years; of 12.6 eq m-2 released, 10.6 neutralised exchangeable acidity
(that CO2 returns to the air), 1.7 sat in solution, 0.3 left as
bicarbonate. The field exported most of its release (32 % of the rock's Ca
and 12 % of its Mg by 2020, 24 to 29 eq m-2) and its exchanger took about
a quarter of what passed through (base saturation +6 points, CEC weighted).

The control's pore water was neutral (32.4 against 32.0 meq L-1 in 2020)
but its leachate was not: anions out 2.05 (NO3) + 3.1 (residual anion)
against cations out 1.83 (Ca+Mg) + 0.47 (K) + 0.24 (Na). The gap is the
Ca and Mg leaching coefficient of 0.4.

## 1. Dissolved Ca and Mg leach at their concentration

`Biogeochemistry_Parameter.m`: `aCa = aMg = 1` (were 0.4). In Porporato et
al. (2003) and T&C-BG (Fatichi et al. 2019) `a` is the fraction of a lumped
pool in solution, which is why NO3 and K carry 1. `P(60)` and `P(68)` are the
dissolved pools by construction: `ERW_Daily_Step_2D` writes them back from
the pore water each day and keeps the exchangeable part in `P(61)` and
`P(69)`. The same coefficient scales passive uptake (`Biogeo_uptake`) and the
2-D lateral mobility (`Biogeo_Lateral_Pools`), which become the dissolved
concentration too.

    LEAK_X = min(B_X, B_X * min(1, Lk/VT))            [g m-2 d-1], X in {Ca, Mg}
    passive uptake_X = min(Ccrown B_X, B_X * T/VT)    [g m-2 d-1]

## 2. The leached background anion returns as the salt it left as

The legacy closure (`keyword_add = 1`) tops the residual anion and Na
inventories up each day by what leached or was taken up, as bare ions. The
ERW books its own An and Na fluxes, so that cycle is alkalinity-neutral.
Once Ca and Mg leave at their concentration, the leachate is neutral and 3
eq m-2 of them leave with the anion over the run; the anion returning alone
would be an acid input of about 0.65 eq m-2 yr-1 in the control, the size
of its whole present acidification. `ERW_Leached_Salt_Resupply` therefore
pairs the leached part with its cations through the existing `I_Ca` and
`I_Mg` boundary slots, booked in `dAlk_boundary`:

    2 I_Ca + 2 I_Mg = f_Lw * max(0, An_tot - Na_tot)   [umol charge m-2 d-1]
    I_Ca : I_Mg = Ca : Mg in the pore water (equal when both are zero)

Both drivers (`ERW_Daily_Step_2D`, `MAIN_FRAME_BG_ERW_10mins`) call it where
they set `I_Na` and `I_An`, so the 1-D/2-D equivalence holds. The uptake
part of the closure keeps the legacy booking; it has no cation counterpart
in the ERW and is neutral on its own. The control's base-cation loss is then
driven by nitrate and bicarbonate leaving, as it should be.

## 3. Alpha re-chosen with base saturation scored

Export acts only on alkalinity in solution, and in m there was almost none
until the exchanger was full. The share the exchanger takes is the H-site
buffer, which in this model is the Appelo alpha. The sweep
(`Calibration_Sweep_1D`) now scores the August base saturation (digitised
Fig. 1F, `Data_Beerling.mat` `CEC_Exp`) beside pH, as the model's August mean
against each sampled interval and pooled, and the basalt-minus-control
difference as a fifth Pareto objective
(`cumulative_N2O_discrete_depth_pH_base_sat_v11`). The alpha grid keeps the
five earlier values and adds 4.5 and 6 (336 trajectories, 112 combinations;
the pilot is the last case). The 0-30 cm base saturation is the soil-mass
average of the two intervals, a linear quantity unlike pH.

## Dimensional check

| Variable | Expected SI | In the code |
|---|---|---|
| `B60`, `B68` | kg m-2 | g m-2 |
| `Lk`, `T` | m d-1 | mm d-1, summed over the day |
| `VT` | m | mm, `Zbio*(Obio-Ohy)` |
| `Lk/VT`, `T/VT`, `aCa`, `aMg` | 1 | 1, `Lk/VT` capped at 1 |
| `LEAK_Ca`, uptake | kg m-2 s-1 | g m-2 d-1 |
| `f_Lw` | s-1 | d-1, `Lk_tot/VSUM` capped at 1 |
| `An_tot`, `Na_tot` | eq m-2 | umol charge m-2 |
| `Ca`, `Mg` in the split | mol m-3 | umol L-1 (ratio only) |
| `I_Ca`, `I_Mg` | mol m-2 s-1 | umol m-2 d-1, as `ERW_function` adds them |
| `dAlk_boundary` | eq m-2 per step | `(2 I_Ca + 2 I_Mg) dtw`, umol charge m-2 |
| base saturation | 1 | per cent of CEC, `100 sum(f_CEC(:,1:4))` |

Molar masses 40, 24.3 and 39.1 g mol-1 are used in both directions of the
T&C to ERW write-back, so the round trip is exact.

## Failure points checked

1. The control acidifies faster with 1 alone (its 3.1 eq m-2 of leached
   anion would strip cations). Item 2 closes that; the remaining loss is the
   nitrate-driven one. `budget.m` on n's control is the check.
2. Passive uptake rises 2.5 times. It is capped at the crown share of the
   pool and litter returns all but the grain, so the net is the grain alone.
3. Larger daily deltas into `ERW_function`. The inventory check, the `R_alk`
   cap retry and the 1-D/2-D equivalence test guard this; every pin, the
   regen baseline and the three-day references change for every run.

Edge cases: `Lk/VT > 1` capped (at most the dissolved pool leaves in a day);
`VT = 0` with `Lk > 0` already raises; `VT = 0` with `T > 0` gives `Inf`
inside `min`, resolved by the crown cap; `An_tot <= Na_tot` returns no salt;
`Ca = Mg = 0` splits equally; negative or nonfinite inputs are refused.

## Tests

Two edited files (`MAIN_FRAME_BG_ERW_10mins.m`, `test_lateral_routing.m`) had
mixed line endings in m; they are now uniformly CRLF, the ending they
mostly had. No other line of either file changed.


`test_leached_salt_closure` (new, in the audit list): leaching of a
dissolved pool returns C x Lk and every dissolved pool sees the same
drainage fraction; passive uptake returns C x T; the salt carries exactly
the leached anion charge Na did not, in pore-water proportions, with the
three edge cases above and the input guard; both drivers call the function.
`test_lateral_routing` carries the new mobility fractions.
`test_bouchet_1d_joint_sweep_configuration` carries the seven-value grid and
336 trajectories.

## Sources

Beerling et al. 2024 PNAS 121:e2319436121 (Methods, Fig. 1E and 1F, SI
Fig. S7); Porporato et al. 2003 Adv. Water Resour. 26:45; Fatichi et al.
2019 GMD 12:4517; Lewis et al. 2021 Appl. Geochem. 132:105023; Kanzaki et
al. 2024 GMD 17:4515; Appelo 1994 Water Resour. Res. 30:2793.
