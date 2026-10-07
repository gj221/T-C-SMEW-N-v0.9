# Package t: the amphibole rate law, 2 October 2026

Package t is package r with one change to the model: ferroactinolite (and
the unused `actinolite` row) dissolves at the measured amphibole rate of
Hermanska et al. (2022) instead of the tremolite row of Palandri and
Kharaka (2004). The dissolution factor stays 1. The 1-D sweep goes back to
package o's three dimensions (alpha, KNmaxd, Kbot) at dissolution factor 1.

## What was wrong

In the five-year basalt run at dissolution factor 1 one mineral supplies
almost all of the silicate alkalinity. Dissolved by 31 August 2020, after
four applications of 50 t/ha:

| Mineral | Share of rock | Dissolved, mol m-2 | Share of the mineral applied | Alkalinity, eq m-2 |
|---|---|---|---|---|
| ferroactinolite | 11.6 % | 1.00 | 37 % | 11.5 |
| epidote | 25.6 % | 0.16 | 1.5 % | 0.64 |
| chlorite | 36.3 % | 0.05 | 0.4 % | 0.33 |
| albite, quartz | 24.8 % | 0 | 0 | 0 |
| calcite | 0.54 % | 1.08 | 100 % | 2.16 |

The exchanger holds 10.3 eq m-2 of acidity (exchangeable H and Al, a fifth
of a 50.5 eq m-2 CEC), so 11.5 eq m-2 from the amphibole alone uses it up
in 2019. The field's exchanger gained about 3.0 eq m-2 by August 2020
(2.0 net of the difference the plots already had in 2016).

The amphibole was run on the tremolite row of Palandri and Kharaka (2004,
Table 27): acid log k -8.40, E 18.9 kJ/mol, n 0.70; neutral log k -10.60,
E 94.4 kJ/mol. That row is not a fit to tremolite experiments. Their text
(section 3.4.2) takes "glaucophane, hornblende, riebeckite and tremolite,
from Sverdrup (1990)", the coefficients of a soil-weathering model, and
only derives the activation energy from his values at 8 and 25 C. The one
amphibole they regressed from experiments, anthophyllite (Chen and
Brantley 1998), has neutral log k -14.24 in the same table.

Against direct measurements, as log10 of the rate in mol of mineral per m2
(BET) per second, near pH 6:

| Source | 25 C | 37 C |
|---|---|---|
| Palandri and Kharaka tremolite row (as coded) | -10.6 | -10.0 |
| Hermanska et al. (2022) fit, tremolite = ferroactinolite | -12.4 | -12.1 |
| Diedrich, Schott and Oelkers (2014), measured, pH 5.3 to 6.9 | | -11.9 to -12.6 |

The coded row is about 60 times the fit at 25 C and 100 to 400 times the
measurements at 37 C.

- Diedrich et al. (2014, Mineral. Mag. 78, 1449) measured steady-state
  tremolite dissolution in flow reactors, BET-normalised: log r -11.8 at
  pH 2 and 25 C, -11.9 to -12.6 at pH 5.3 to 6.9 and 37 C. They conclude
  that tremolite releases Mg about three orders of magnitude more slowly
  than forsterite and enstatite.
- Hermanska et al. (2022, Chem. Geol. 597, 120807) give acid log k25 -11.3,
  Ea 50 kJ/mol, n 0.2 and neutral log k25 -13.1, Ea 48 kJ/mol. The PHREEQC
  rates database (`phreeqc_rates.dat`, block `RATE_PARAMETERS_HERMANSKA`)
  lists the same numbers for Tremolite, Ferroactinolite and Riebeckite. The
  fit reproduces Diedrich's pH 2 value (-11.7 against -11.8) and lies
  inside the spread of their near-neutral values.
- Critelli et al. (2014, Chem. Geol. 390, 100) dissolved a metabasalt of
  the Blue Ridge assemblage (actinolite and chlorite). Their abstract, the
  only part read here, reports chlorite 3 and 5 times faster than
  actinolite at pH 2.0 and 3.2, and the two releasing Si at the same rate
  at pH 7.2. The coded constants made the amphibole 80 times faster than
  chlorite at neutral pH.

The "fast-weathering" label of the Blue Ridge amphibole in Lewis et al.
(2021, Table 2) and Beerling et al. (2024) rests on the same Palandri and
Kharaka row.

## The change

`v2/T&C_Code/Parameters/ERW_MOD_PARAM_excel.m`, function `min_const`, rows
`ferroactinolite` and `actinolite`:

| | before | after |
|---|---|---|
| `k_diss_H`, mol m-2 s-1 at pH 0 | 10^-8.4 | 10^-11.3 |
| `E_H`, kJ/mol | 18.9 | 50 |
| `n_H` | 0.7 | 0.2 |
| `k_diss_w`, mol m-2 s-1 | 10^-10.6 | 10^-13.1 |
| `E_w`, kJ/mol | 94.4 | 48 |

Stoichiometry, molar mass, surface area, the rate equation and every other
mineral are untouched. `test_erw_mineral_kinetics` now pins the two rows.
Control and calcite runs hold no amphibole and are r's bit for bit.

## What it does to the release

Offline, from the saved factor-1 run and a low-release run (their pH,
temperature, moisture and surface area fed to the model's own rate
equation, which reproduces the stored rates to four digits). No model was
run. Cumulative to 31 August, dissolution factor 1:

| | 2017 | 2018 | 2019 | 2020 |
|---|---|---|---|---|
| ferroactinolite before, mol m-2 | 0.10 | 0.37 | 0.70 | 1.01 |
| ferroactinolite after, mol m-2 | 0.003 | 0.01 | 0.02 | 0.03 |
| silicate alkalinity before, eq m-2 | 1.3 | 4.5 | 8.7 | 12.6 |
| silicate alkalinity after, eq m-2 | 0.12 | 0.43 to 0.46 | 0.85 to 0.96 | 1.3 to 1.5 |
| with the rock's calcite, before, eq m-2 | 1.8 | 5.6 | 10.3 | 14.8 |
| with the rock's calcite, after, eq m-2 | 0.66 | 1.5 | 2.5 to 2.6 | 3.4 to 3.7 |
| field exchanger gain, eq m-2 | 1.2 | 2.3 | 3.0 | 3.0 |
| field gain net of the 2016 difference, eq m-2 | 0.3 | 1.3 | 2.0 | 2.0 |

The ranges are the factor-1 run's pH (low end) and the control's pH (high
end); the corrected run will sit between. (The 14.8 is gross; 1.9 eq m-2 of
it re-precipitated as pedogenic calcite in 2020, which is why earlier notes
quote 12.8.)

With about 3.5 eq m-2 by August 2020 the exchanger keeps most of its
acidity, no pedogenic carbonate forms in five years, and the pH effect
follows the exchange slope, (0.434/f_H + alpha)/100 pH per point of base
saturation. Estimated basalt-minus-control effect, August:

| | 2017 | 2018 | 2019 | 2020 |
|---|---|---|---|---|
| field | +0.18 | +0.28 | +0.45 | +0.47 |
| before, alpha 1.3 (measured) | +0.13 | +0.41 | +1.06 | +2.63 |
| after, alpha 1.3 (estimate) | +0.05 | +0.11 | +0.18 | +0.24 |
| after, alpha 3.4 (estimate) | +0.07 | +0.16 | +0.27 | +0.37 |
| after, alpha 4.5 (estimate) | +0.09 | +0.20 | +0.32 | +0.44 |
| after, alpha 6 (estimate) | +0.11 | +0.24 | +0.40 | +0.56 |

These are estimates, checked against one stopped local run with nearly the
same release (alpha 4.5, +0.09 and +0.22 in 2017 and 2018). The sweep
gives the measured values. Alpha 1.3 was the earlier sweeps' choice because
a flat pH response hid twice the base in 2017 and 2018; with the release
right the choice should move to 3.4 or above, the range the field's own pH
per point of base saturation (0.06 to 0.08) asks for.

## What it does not fix

The control. Its August pH is 0.15 to 0.5 below the field's at alpha 1.3
and 0.26 to 0.67 below at alpha 4.5, because its base saturation falls from
78.4 to 74.2 % where the field's stays near 79 %. A steeper alpha makes
that offset larger, so the sweep has to weigh the basalt effect against the
control level, and the basalt plot inherits the same offset in absolute pH.
No rock chemistry touches the control.

Its budget (saved control run, alpha 1.3, 1 September to 31 August):

| | 2017 | 2018 | 2019 | 2020 |
|---|---|---|---|---|
| change in total alkalinity, eq m-2 | -0.63 | -0.76 | +0.25 | -0.53 |
| nitrate leached, mol N m-2 | 0.41 | 0.56 | 0.69 | 0.31 |
| nitrate leached, kg N/ha | 57 | 78 | 97 | 43 |

The acid is the nitrate leached, 1.7 eq m-2 of alkalinity lost in four
years. The one measured year, 2017, gave 67 kg N/ha, so the acid load is
not the error. What the model lacks is any base supply to the top 30 cm:
net deposition is zero, the soil's own minerals do not weather, the crop
takes every cation from the same 30 cm its residue returns to, and no lime
remains from the last liming of these plots in 2012 (Kantola et al. 2023).
The 2016 pre-application sample (6.07) is also the lowest of the five
control years (6.07, 6.28, 6.31, 6.13, 6.14), and the run is seeded on it.
This is the next thing to resolve; it is not addressed here.

## Consequences to know before adopting it

- Silicate weathering in the five-year window falls to about one ninth
  (12.5 to 1.4 eq m-2). About 60 % of the basalt's alkalinity in 2017 to
  2020 is then its 0.54 % calcite (2.2 of 3.5 eq m-2).
- The 35-year runs still accumulate rock, so the silicate release grows
  each year and the exchanger still saturates, years later than before.
- The rock-loss estimate of Beerling et al. (2024; 32 % of the rock's Ca
  and 12 % of its Mg lost at 0-10 cm, 10.5 t CO2/ha) is not reproduced and
  was not before. Amphibole and calcite hold 23 % of the rock's Ca; the
  rest is epidote, which no rate database dissolves that fast, and a
  release that large cannot leave a soil at pH 6.6 as bicarbonate.

## Sensitivities checked offline

- All three reactive silicates on Hermanska et al. (2022, 2023): 3.9 to
  4.3 eq m-2 by August 2020 instead of 3.4 to 3.7 (epidote and chlorite
  have lower activation energies there and are faster in a cold soil).
- The amphibole on the Hermanska hornblende law (acid log k -10.7): 4.4 to
  4.7 eq m-2, the amphibole giving 1.2 to 1.4 eq m-2 instead of 0.3 to 0.4.

Neither changes the picture; both are left to the dissolution factor of the
uncertainty ensemble (0.8 and 1.2 around 1), which now brackets a central
value the field supports.

## The sweep in this package

Package o's grid at dissolution factor 1: alpha 1.3, 1.825, 2.35, 2.875,
3.4, 4.5, 6; KNmaxd 0.0048, 0.0016, 0.0032, 0.0064 d-1; Kbot free drainage,
0.07, 0.03, 0.01 mm h-1; runs 10 to 12. 112 combinations, 336 trajectories,
about 840 CPU-hours, about 18 hours at 48 concurrent slots. Case 1 to 3 are
the reference tuple; the pilot is case 336 (alpha 6, KNmaxd 0.0064, Kbot
0.01, calcite). Tags keep the `_fd1` of packages p to s. The control and
calcite trajectories repeat r's physics; only the basalt ones are new.

Read `Front_pH_difference_base_sat_difference` first, then the control's
`pH_season_control_RMSE` on the same rows.

## Files changed from r

`v2/T&C_Code/Parameters/ERW_MOD_PARAM_excel.m` (two rows of `min_const`);
`v2/T&C_Code/Tests/test_erw_mineral_kinetics.m` (pins them);
`runtime/Calibration_Sweep_1D/Sweep_Cases.m`, `Run_Sweep_Task.m`,
`Validate_Sweep.m`, `Summarize_Sweep.m`, `submit_sweep.sh`, `README.md`
(o's grid, counts and smoke cases, with the dissolution factor carried as a
single value); `test_bouchet_1d_joint_sweep_configuration.m`; the root
README, the runbook (section 0t), the dissolution-factor note and the audit
document's reference list. The pins are regenerated.

## Validation

On Bouchet, through the sweep's own `SweepCheck` job (the full audit and
the three-case smoke), then the one-year control gate, whose residual is
o's (-0.0248) because the control holds no rock. `ERW_MOD_PARAM_excel.m` is
inside the equivalence fingerprint, so this needs a fresh root; the
controller builds the baseline as its first job. Locally, as a check of the
edits only, all passing: `checkcode` on the edited files; the 38 audit and
campaign tests that read an edited file or run the rock (the driver
equivalence among them, on a local baseline); a replica of the three-case
smoke; and a ten-day control and basalt pair from 8 November 2016 run in
the r and the t tree. The control pair is identical in every model
variable. In the basalt pair (soil 9 to 12 C, pore pH 5.65) the amphibole
dissolves 24 times more slowly in t, the stored rate equals the new law
evaluated at the run's own state, and epidote, chlorite and calcite
dissolve the same to four digits.
