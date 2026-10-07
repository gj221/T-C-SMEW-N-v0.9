# Package r: the soil-water paste also dissolves the plot's Mg carbonate, 1 October 2026

Package r is package q (`PASTE_CALCITE_2026-10-01.md`,
`DISSOLUTION_FACTOR_SWEEP_2026-10-01.md`) with one further change to the
bulk-pH operator: the plot's solid Mg carbonate, the nesquehonite pool of
`MG_CARBONATE_REVIEW_2026-09-30.md`, enters the soil-water paste and
dissolves there, as the plot's calcite does since q. Nothing else in the
model differs from q.

The notes of this package also say "dissolution factor" where p and q said
"release factor". It is the same quantity (`DissolutionFactor`, the
workbook `diss_f`); no option, column, tag or number changes.

## Why

The reported bulk pH is the pH of a 1:2.5 soil-water paste. Package q let
calcite dissolve in that paste and left the Mg carbonate frozen, which the
Mg carbonate review had already named as a bias: "in a real extract some of
it would redissolve and raise the pH, so once the pool is large the
modelled bulk pH of the basalt soils is a lower bound". That is the wrong
way round for the two solids. The solubility product of nesquehonite
(MgCO3.3H2O) is about 1,600 times calcite's, and it comes to equilibrium
with water within hours (Harrison et al. 2013). A soil sample that holds it
carries it into the paste, the pore water is diluted about ninefold, and
the solid dissolves until the slurry is saturated again or the pool is
spent. A measurement on that soil reads the result, so the operator now
computes it.

The pool is exactly zero in the controls, the calcite runs and every
five-year basalt run (no workbook case seeds it, and it forms only once
calcite has taken the Ca). It forms in the 35-year basalt runs, after
about 2031 at dissolution factor 1 (`MG_CARBONATE_REVIEW_2026-09-30.md`)
and presumably later at a lower factor. Those runs are the only ones this
change moves.

## What changed

`Porewater_to_Bulk_pH.m`: the caller passes the pool as
`bulk_pH_params.solid_MgCO3_umol_m2`. At a trial extract H the slurry takes
x_Ca of calcite and x_Mg of nesquehonite [umol per litre of extract
water]: Ca grows by x_Ca, Mg by x_Mg and the DIC inventory by both. Each
solid obeys the rule q gave calcite, and nothing precipitates in the paste:

    x = S                      the whole solid dissolves and Omega <= 1
    0 < x < S                  part dissolves and Omega = 1
    x = 0                      none dissolves because Omega >= 1 already

with Omega = gamma_div^2 M CO3 / K, K_CaCO3 from Plummer and Busenberg
(1982) and log K_MgCO3 = -5.27 (Harrison et al. 2019), both at the 25 C of
the paste, and gamma_div from the slurry's ionic strength (the activity
basis of `ERW_function`). The two solids share the carbonate ion, so the
two extents are settled together (`charge_residual_solid_carbonates`). The
charge-balance residual is the existing closed-DIC residual at that pair.

A paste without Mg carbonate takes q's code path unchanged, so every run in
which the pool is zero reproduces q bit for bit.

`ERW_Daily_Step_2D.m` and `MAIN_FRAME_BG_ERW_10mins.m` pass the pool
(`MgCO3`, umol m-2) beside the calcite inventory.

`Tests/test_pH_conversion.m` gains the checks listed under Validation.

## Dimensional check

| Variable | Expected SI | In the code |
|---|---|---|
| nesquehonite pool | mol m-2 | umol m-2 (`bg.MgCO3`) |
| extract water | m3 m-2 | L m-2, `n_extract * rho_bulk * Zrw` |
| S_Mg, x_Mg, S_Ca, x_Ca | mol m-3 | umol L-1 of extract water |
| K_MgCO3 | mol2 m-6 | umol2 L-2, `10^(-5.27) * conv_mol^2` |
| gamma_div | 1 | 1, Davies at the slurry's ionic strength |
| saturation ratio | 1 | `Mg * CO3 / (K / gamma_div^2)` |
| Mg inventory of the slurry | mol m-3 | `Mg_d + cec * f_Mg / 2 + x_Mg`, umol L-1 |
| DIC inventory of the slurry | mol m-3 | `DIC_d + x_Ca + x_Mg`, umol L-1 |

## Numerics

At a fixed H either saturation ratio rises with either extent: the cation
and the DIC inventories grow, and the added cation displaces the other
from the exchanger. The calcite extent that meets its own rule therefore
falls as more nesquehonite dissolves, and the nesquehonite ratio taken
along that calcite extent still rises. The pair is one bracketed root for
x_Mg with the calcite extent resolved inside it. Both roots use Brent's
method on the logarithm of the saturation ratio, closed to 1e-13 of the
solid, and return the undersaturated end of the last bracket, so the solid
remaining is never negative and dissolution never leaves the slurry
supersaturated.

The carbonate dissolved at a trial H adds x(2 CO2 + HCO3)/DIC of positive
charge, more at low pH, so the residual rises with H as the solid-free one
does and has one root. Two consequences are used. The narrow scan around
the previous day's root visits its grid outward from the anchor and stops
at the first sign change (`bulk_pH_residual_scan_outward`); the bracket is
the one the whole scan selects, and on ten test states the result is
bit-identical to the whole scan. And the check that looks for a second,
pore-proximal root when the paste stands more than 2.5 pH units above the
pore water is skipped for a paste that dissolved Mg carbonate, where such
a distance is chemistry and there is no second root.

Cost per daily call, measured on a loaded PC: 0.02-0.05 s without solids,
0.8 s with calcite only (q's path, unchanged), 0.01-0.2 s with Mg carbonate
alone or with a pool the paste dissolves completely, 1.3 s when both
solids dissolve in part. A 35-year basalt run that holds the pool for
twenty years spends at most about 2.5 hours on it.

## What it does to the reading

Synthetic pore waters at 15 C, saturated with calcite and nesquehonite,
with the campaign's exchanger and extract (CEC 50.5 eq m-2, alpha 1.3,
1:2.5). They are constructed states, not output of a 35-year run.

| soil pCO2 | pore pH | paste, solids frozen (p) | calcite dissolves (q) | pool 0.1 | 0.5 | 1 | 2 | 3 mol m-2 or more |
|---|---|---|---|---|---|---|---|---|
| 0.3 % | 8.77 | 8.85 | 9.19 | 9.24 | 9.43 | 9.61 | 9.83 | 9.98 |
| 1 % | 8.44 | 8.55 | 8.89 | 8.95 | 9.16 | 9.35 | 9.60 | 9.71 |
| 3 % | 8.15 | 8.24 | 8.57 | 8.63 | 8.87 | 9.09 | 9.36 | 9.43 |

The last five columns are package r with 20 mol m-2 of calcite beside the
nesquehonite pool shown. About 2.4 to 3.1 mol m-2 of nesquehonite
saturates the paste; a larger pool reads the same. Calcite then dissolves
very little (0.006 mol m-2) because the nesquehonite has already raised
the carbonate.

## Limits

- The paste stays closed to the atmosphere, the convention of every plot.
  A slurry at pH 9.5 to 10 would take up CO2 if it stood open, so the
  reading is the first reading of such a sample. (Package u opens it:
  `PASTE_OPEN_TO_AIR_2026-10-02.md`. A paste saturated with nesquehonite
  then reads about 9.3.)
- The model has no aqueous ion pairs (MgCO3 and MgHCO3+), in the pore water
  or in the paste. The saturated reading does not depend on them, because
  saturation fixes the free-ion product; the amount of solid the paste
  needs to reach it is larger with pairs than the 2.4 to 3.1 mol m-2 here.
- The bulk pH is also the pH the biological response functions take
  (`BIOGEOCHEMISTRY_DYNAMIC3`). In a soil that holds the pool the reading
  now stands up to 1.2 units above the frozen one, and the nitrogen and
  carbon responses follow it. At the 1 % CO2 state above, the
  decomposition factor is 0.77 at the frozen reading (8.55), 0.69 at q's
  (8.89) and 0.51 at 9.71; the N2:N2O ratio of denitrification is 9, 12
  and 28. The pore water these soils hold is at pH 8.2 to 8.8. Whether the
  organisms should answer to the paste or to the pore water in such a soil
  is a modelling choice this package does not change: every plot's biology
  reads the paste.
- This change reads the state the model holds; it does not decide whether
  the pool should exist. The pool forms in the 35-year basalt runs because
  the 0-30 cm box keeps nearly all the Mg and alkalinity the rock releases
  (`LEACHING_CLOSURE_REVIEW_2026-10-01.md`,
  `DISSOLUTION_FACTOR_SWEEP_2026-10-01.md`).
- A sample with unreacted Mg carbonate and an acid exchanger is not a state
  the model holds: the pool exists only where the exchanger has no acidity
  left. The operator's answer there is still worth knowing. A small pool
  dissolves completely, its carbonate ends as CO2 trapped in the closed
  slurry, and the paste reads up to 0.5 below the solid-free one (q's note
  describes the same for calcite). A large pool neutralises the sample's
  exchange acidity inside the paste and reads 9.6 or more while the pore
  water is still acid, as in the seeded check under Validation.
- With the legacy `TC_GT` exchange scheme the Gaines-Thomas partition can
  fail when Mg stands far above Ca (seen for a pool of 20 mol m-2 or more
  in an acid test sample). The operator then stops with
  `ERW:noBulkPHRoot`; it does not return a wrong reading. The production
  scheme (`appelo`) has a partition for any inventory.

## Validation

On Bouchet, through the sweep's `SweepCheck` job (the full audit, then the
three-case smoke) and the five-year calcite pilot, as for q.

Locally, as a check of the edits only (MATLAB R2026a): `checkcode` on the
edited files; the 30 audit tests that call the operator or read the text
of a file this package edits, the wording edit included; the driver
equivalence test on a baseline built for it; a replica of the three-case
smoke; and a six-day calcite run and a six-day basalt run across the first
application, each identical to q's in every saved model variable. No
five-year run holds the pool, so the drivers were also exercised on a
throwaway copy that seeds 20 mol m-2 of nesquehonite into the control soil:
over 25 days the pool dissolved to 13 mol m-2, and on every day the stored
bulk pH equalled the operator's reading of that day-end state with that
day's pools to 3e-9. The full audit was not run locally.
`test_pH_conversion` now also checks that a paste without Mg carbonate
reproduces the calcite operator bit for bit, that a trace of nesquehonite
dissolves completely, that an excess leaves the paste at nesquehonite
saturation with solid remaining and the DIC and Mg inventories raised by
exactly what dissolved, that an excess of both carbonates settles at both
saturations at once on a single root, that a small pool dissolves
completely while calcite stays saturated, that calcite left supersaturated
is not touched, and that a tracked call reaches the same root through the
outward scan.

## Sources

* Harrison, A.L., Mavromatis, V., Oelkers, E.H., Bénézeth, P., 2019.
  Solubility of the hydrated Mg-carbonates nesquehonite and dypingite from
  5 to 35 C. Chemical Geology 504, 123-135. doi:10.1016/j.chemgeo.2018.11.003
* Harrison, A.L., Power, I.M., Dipple, G.M., 2013. Accelerated carbonation
  of brucite in mine tailings for carbon sequestration. Environmental
  Science & Technology 47, 126-134. doi:10.1021/es3012854
* Plummer, L.N., Busenberg, E., 1982. The solubilities of calcite,
  aragonite and vaterite in CO2-H2O solutions between 0 and 90 C.
  Geochimica et Cosmochimica Acta 46, 1011-1040.
  doi:10.1016/0016-7037(82)90056-4
* Kanzaki et al. 2024, the closed-DIC laboratory extract of SCEPTER, as
  cited in `Porewater_to_Bulk_pH.m`.
