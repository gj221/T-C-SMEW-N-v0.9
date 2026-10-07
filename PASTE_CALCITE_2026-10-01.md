# Package q: the soil-water paste dissolves the plot's calcite, 1 October 2026

Package q is package p (the dissolution-factor sweep,
`DISSOLUTION_FACTOR_SWEEP_2026-10-01.md`) with one change to the model: the
bulk-pH operator lets the plot's solid calcite dissolve in the soil-water
paste. Nothing else differs from p.

## Why

The 50 t/ha calcite run of package o (run 12, workbook case 7, 5000 g m-2
on each November date) lifts the bulk pH from 6.10 to 7.4-7.8 within three
weeks and then holds it there whatever more is applied. The pore water
does what it should: the first application neutralises the exchanger
(5.3 of 50 mol m-2 dissolve, f_H 20 % to 1 %), the water saturates with
calcite, and the pore pH sits at the calcite-CO2 equilibrium for the soil's
pCO2 and Ca (7.0-7.25 at 1-2 % CO2 with 6 mM Ca). That equilibrium is
thermodynamics and is not changed.

The cap is in the reading. The reported bulk pH is the pH of a 1:2.5
soil-water paste, and `Porewater_to_Bulk_pH` followed SCEPTER (Kanzaki et
al. 2024) in freezing every field solid in the paste "because of the short
duration for soil pH measurements". That holds for silicates. It does not
hold for calcite, which dissolves within the measurement at the model's own
surface area: a soil sample holding 45 mol m-2 of free calcite carries it
into the paste, the pore water is diluted about eightfold, and the calcite
dissolves until the slurry is saturated again. A soil with free carbonate
reads 7.8-8.3 in a water paste for that reason. With the solid frozen the
operator could only report the diluted pore water, 0.4-0.5 above the pore
pH, and the calcite plot looked capped near the basalt plot's end state.

## What changed

`Porewater_to_Bulk_pH.m`: the solid calcite inventory the caller passes
(`bulk_pH_params.solid_CaCO3_umol_m2`, previously forced to zero inside the
operator) now enters the closed slurry. At a trial extract H the amount
dissolved, xi [umol per litre of extract water], raises the Ca and the DIC
inventories by xi each, and

    xi = S                         if the slurry is undersaturated with all of it dissolved
    xi = 0                         if it is already saturated before any dissolves
    gamma_div^2 Ca(xi) CO3(xi) = K_CaCO3(25 C)     otherwise

with K from Plummer and Busenberg (1982) at the 25 C of the paste and
gamma_div from the slurry's ionic strength (the activity basis of
`ERW_function`). Nothing precipitates in the paste. The residual of the
charge balance is the existing closed-DIC residual evaluated at that xi
(`charge_residual_solid_calcite`); without solid calcite the operator takes
the old path bit for bit, so the control, the seed inversion and the
equivalence baseline are unchanged.

`ERW_Daily_Step_2D.m` and `MAIN_FRAME_BG_ERW_10mins.m` pass the inventory:
the applied grains not yet dissolved (`M_min` of the calcite column over its
molar mass) plus the pedogenic pool (`CaCO3`), in umol m-2. The two
statements are the same in both drivers.

`Tests/test_pH_conversion.m`: the block that pinned an inert solid now
checks that an excess of calcite leaves the paste at saturation with part
of the solid remaining and the DIC raised by exactly what dissolved, that a
trace dissolves completely, and that a zero inventory reproduces the
solid-free operator bit for bit.

## Dimensional check

| Variable | Expected SI | In the code |
|---|---|---|
| solid calcite inventory | mol m-2 | umol m-2 (`M_min` g m-2 over `MM_min` g umol-1, plus `CaCO3` umol m-2) |
| extract water | m3 m-2 | L m-2, `n_extract * rho_bulk * Zrw` |
| S, xi | mol m-3 | umol L-1 of extract water |
| K_CaCO3 | mol2 m-6 | umol2 L-2, `10^(logK) * conv_mol^2` |
| gamma_div | 1 | 1, Davies at the slurry's ionic strength |
| saturation ratio | 1 | `Ca * CO3 / (K / gamma_div^2)` |

## Numerics

At a fixed H the saturation ratio rises monotonically with xi (both
inventories grow), so xi is bracketed on [0, S] and bisected to 1e-13 of S;
the undersaturated end of the last bracket is returned, so the remaining
solid is never negative. The calcite dissolved at a trial H adds
xi(2 CO2 + HCO3)/DIC of positive charge to the residual, more at low pH, so
the residual stays monotonic in H and the bracket the existing scan selects
is the root; the strict full-range rescan the dead solid branch used is
removed. Cost: about 45 evaluations of the exchange state per residual
call, on days when solid calcite is present.

## Evidence (package o calcite run, states recomputed under four closures)

| day after 1 Nov 2016 | pore pH | closed, solid frozen (o) | closed, calcite dissolves (q) | open to air, frozen | open to air, calcite dissolves |
|---|---|---|---|---|---|
| 5 (before application) | 5.62 | 6.09 | 6.09 | 6.11 | 6.11 |
| 11 (day after) | 6.34 | 6.65 | 7.08 | 7.31 | 8.14 |
| 25 | 6.91 | 7.28 | 7.79 | 7.66 | 8.14 |
| 45 | 7.03 | 7.39 | 7.85 | 7.77 | 8.12 |
| 75 | 7.39 | 7.78 | 8.32 | 7.78 | 8.14 |
| 119 | 7.01 | 7.41 | 7.95 | 7.66 | 8.13 |

The calcite plot now reads 7.8-8.3 once its exchanger is neutralised. The
acid states do not move (no solid). Opening the paste to atmospheric CO2
would give the classical 8.1 and would move the acid states by only 0.01,
but it changes the convention for every plot and is not part of this
package.

Two properties of the closed slurry are worth knowing. With an acid
exchanger and unreacted calcite in the same sample, the carbonate that
neutralises exchange acidity stays in the slurry as CO2, so the paste rises
less than an open beaker would (pH 6.50 in the unit test's deliberately
unreacted state). And a trace of calcite in an acid paste lowers the pH by
a few thousandths for the same reason. Neither state persists in the model:
the field neutralises the exchanger within weeks and the rock's own calcite
is gone within days.

## What it does not change

The pore-water chemistry, the calcite kinetics and the per-step saturation
limiter, the exchanger, the hydrology and the export operator are p's. The
basalt runs see the paste change only while the rock's own calcite
(0.27 mol m-2 per application) is still solid, a matter of days, and where
pedogenic calcite has formed. That needs the plot's exchanger acidity to be
used up first, so a combination that keeps some acidity through August 2020
scores as in p, and one that runs out reads higher than in p from then on.
(Package r adds the Mg carbonate pool to the same paste:
`PASTE_MG_CARBONATE_2026-10-01.md`. Package u opens the paste to the
laboratory air, the fourth closure of the table above:
`PASTE_OPEN_TO_AIR_2026-10-02.md`.)

## Validation

On Bouchet, through the sweep's `SweepCheck` job (the full audit, which
includes `test_pH_conversion` and `test_calcite_event_smoke`, then the
three-case smoke) and the five-year calcite pilot that precedes the
production array. Locally only the edited test, the calcite event smoke and
a replica of the three-case smoke were executed, with a 120-day calcite run
and a 75-day basalt run as model-level checks; no audit was run locally.
