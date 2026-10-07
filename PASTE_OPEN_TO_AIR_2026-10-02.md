# Package u: the soil-water paste stands open to the laboratory air, 2 October 2026

Package u is package t (`AMPHIBOLE_KINETICS_2026-10-02.md`) with one change
to the model and two to the campaign. The model change is in the bulk-pH
operator: the 1:2.5 soil-water paste is open to the laboratory air, so its
dissolved CO2 is the air's, and the plot's solid carbonates still dissolve
in it (`PASTE_CALCITE_2026-10-01.md`, `PASTE_MG_CARBONATE_2026-10-01.md`).
The campaign changes are the Appelo alpha, 3.4 in place of 1.3, and the
uncertainty ensemble's alpha range, 1.3 to 6 around it; they are at the end
of this note.

## Why

The reported bulk pH is the pH of a 1:2.5 soil-water paste at 25 C. Until
t the operator carried the field pore water's dissolved inorganic carbon
into the paste and closed it to gas, after SCEPTER (Kanzaki et al. 2024).
Three things are against that closure for this measurement.

- The sample is air-dried and sieved before it is wetted (Beerling et al.
  2023 preprint, quoted in `runtime/Base_2D/README.md`). The CO2 its pore
  water held in the field, several to about fifty times the atmosphere's,
  leaves on the bench. It is not there to be carried into the paste.
- The paste is mixed in a vessel that holds air. Its dissolved CO2 moves
  toward the air's, not toward the field's.
- SCEPTER's own code is not consistent on the point. Its released default
  closes the laboratory run and carries the field DIC, but the path that
  does not track DIC holds the paste at atmospheric CO2
  (`scepter.f90`, untracked gases stay at their boundary value), and its
  spin-up scripts add lime fully dissolved. There is no convention to
  inherit, so the operator follows the measurement.

A closed paste on a Kantola soil state held 7 to 98 umol/L of CO2 through
the year against about 15 in a paste at equilibrium with air. For an acid
soil that barely moves the reading. For a soil that holds free carbonate
it decides it: closed, the CO2 that calcite releases while it neutralises
exchange acidity stays in the slurry and the paste read as low as 6.5 with
unreacted calcite in it; open, that CO2 leaves and the paste goes to the
calcite-air equilibrium, the 8.1 to 8.3 a calcareous soil reads in water.

## What changed

`Porewater_to_Bulk_pH.m`. The dissolved CO2 of the paste is fixed at

    CO2_lab = k_H(25 C) * CO2_air_extract

where `CO2_air_extract` is the atmospheric CO2 the callers already passed
(the forcing's value of the day, ppm/22.41, in umol per litre of air) and
k_H = 0.83 is the dimensionless Henry constant of `ERW_function` at its
25 C reference: the dissolved CO2 the model gives any water at equilibrium
with that air at 25 C. At a trial extract H,

    HCO3 = k1 * CO2_lab / H        CO3 = k2 * HCO3 / H

with the constants conditioned at the paste's own ionic strength as
before, and the charge balance is solved for H. The field DIC no longer
enters the solution. It is still handed over, and the operator reports the
carbon the paste exchanged with the air,
`diag.CO2_air_exchange_umol_L = DIC_paste - DIC_field * dilution - carbonate dissolved`,
negative when the paste degasses.

Only the dissolved gas is reset. Everything else the sample held is carried
into the paste as before: the dissolved cations and anions, the exchanger
with its base cations and its acidity, and the solid carbonates. The
alkalinity of the field soil is therefore kept, through the charge balance;
what the paste loses is the excess CO2 that the soil air held against it.
The field's CO2 stays where it acts, in the pore-water pH and in the
weathering and exchange chemistry that run on it.

The carbonate solids follow the rule of q and r unchanged: each dissolves
until the paste is saturated with it or the solid is spent, and nothing
precipitates. Only the bookkeeping differs: a dissolving solid adds its
cation to the conserved inventory, and the carbonate ion is set by H and
the air instead of by a conserved DIC.

Two smaller changes go with it. The check that looks for a second root
when the paste stands more than 2.5 pH units above the pore water is
skipped for a paste that dissolved calcite, as r already skipped it for Mg
carbonate: an acid sample with unreacted calcite reads above 8 by chemistry.
And under the legacy Gaines-Thomas exchange scheme (`TC_GT`) the exchanger
of an acid sample has no partition at pH 8 until some calcite has
dissolved, because the cations the sample holds cannot fill its sites; the
calcite extent solve now treats such a slurry as short of Ca and moves its
bracket up, where it used to stop with `ERW:partitionExchange`. The
production scheme (`appelo`) has a partition for any inventory and never
takes that path.

`Invert_Bulk_pH.m` calls the forward operator, so the seed is inverted
through the open paste; only its comments changed. `ERW_Daily_Step_2D.m`,
`MAIN_FRAME_BG_ERW_10mins.m` and `ERW_MOD_PARAM_excel.m` are unchanged:
they already passed the day's atmospheric CO2 and the field DIC.

`Tests/test_pH_conversion.m` pins the open paste (under Validation).

## Dimensional check

| Variable | Expected SI | In the code |
|---|---|---|
| atmospheric CO2 | mol m-3 of air | umol per litre of air, `Ca`/22.41 with `Ca` in ppm |
| Henry constant, 25 C | 1 (water over air) | 0.83, `ERW_function`'s `k_H` at `T_ref` |
| dissolved CO2 of the paste | mol m-3 | umol L-1, `k_H_lab * CO2_air_extract` |
| HCO3, CO3 | mol m-3 | umol L-1, `k1 CO2/H`, `k2 HCO3/H`, constants in umol L-1 |
| carbon exchanged with the air | mol m-3 of extract water | umol L-1, `DIC - DIC_d - x_Ca - x_Mg` |
| saturation ratio | 1 | `Ca * CO3 / (K / gamma_div^2)`, unchanged |

## Numerics

At a fixed H the air fixes the carbonate ion, so either saturation ratio
rises with the amount dissolved through the cation alone and the extent
solves of q and r keep their brackets. Bicarbonate and carbonate fall with
H and the dissolved cations rise with it, so the charge residual rises
with H and has one root; a cold start on a calcite-bearing sample finds one
sign change over pH 2 to 12 in both exchange schemes. The cost per daily
call is that of r: a few hundredths of a second without solids, about a
second with calcite, about two with both carbonates.

## What it does to the reading

Constructed states of the unit test (a 1:2 extract on 230 kg m-2 of soil,
air at 400 ppm):

| sample | closed (t) | open (u) |
|---|---|---|
| acid soil, pore pH 5.2, 18 % of sites H, no solid | 5.926 | 5.902 |
| the same with the field DIC doubled | 5.908 | 5.902 |
| the same with 50 mol m-2 of calcite, Appelo | 6.243 | 8.203 |
| the same with 50 mol m-2 of calcite, Gaines-Thomas | 6.500 | 8.218 |
| alkaline Mg-rich soil, pore pH 8.4, no solid | 8.996 | 8.710 |
| the same with 20 mol m-2 of calcite | 9.338 | 8.861 |
| the same with 20 mol m-2 of nesquehonite | 10.081 | 9.284 |

Doubling the air's CO2 lowers the calcite paste by 0.198, two thirds of
log10(2), the slope of the calcite-CO2-water system.

Saved five-year soil states read under both closures: the same pore water,
exchanger and solids, every third day of three local runs of earlier
packages, 568 states each, read by this package's operator as the drivers
call it and set beside the reading stored with the run. The first two
columns are August means for 2016 to 2020.

| run | closed (stored) | open (u) | single days, open minus closed |
|---|---|---|---|
| control, alpha 1.3 | 6.049, 5.921, 5.793, 5.977, 5.836 | 6.052, 5.921, 5.797, 5.977, 5.836 | -0.002 to +0.037, median +0.009 |
| basalt at a low release (alpha 4.5, factor 0.25) | 6.015, 5.942, 5.948, 6.432, 6.474 | 6.018, 5.941, 5.955, 6.431, 6.475 | -0.006 to +0.105, median +0.016 |
| basalt before the amphibole fix, which holds pedogenic calcite in 2020 | 6.049, 6.047, 6.200, 7.031, 8.475 | 6.052, 6.046, 6.209, 7.018, 8.119 | -0.70 to +0.38 |

So an acid soil reads the same August mean to within 0.01 and differs by
up to 0.04 to 0.10 on the days its soil air is richest in CO2, when the
closed paste was the lower of the two; in the run that comes to hold
calcite the two differ by -0.7 to +0.4. None of the 1,704 states failed to
solve; a call took 0.007 s in the median and at most 0.2 s.

Twenty-day runs of this package from 8 November 2016, across the first
application (a local check of the edits; the daily bulk pH is 6.100 in
every run until the application):

| run | pore pH at the seed | bulk pH, first day after | bulk pH, day 17 | pore pH, day 17 | solid calcite, day 17 |
|---|---|---|---|---|---|
| basalt, campaign set | 5.618 | 6.160 | 6.157 | 5.689 | none |
| calcite 50 t/ha, campaign set | 5.618 | 8.129 | 8.126 | 6.895 | 45.4 of 50 mol m-2 |
| calcite 5 t/ha, campaign set | 5.618 | 7.816 | 7.812 | 6.892 | 0.40 of 5 mol m-2 |
| basalt, lower set | 5.556 | 6.139 | 6.136 | 5.620 | none |
| calcite 50 t/ha, lower set | 5.556 | 8.125 | 8.123 | 7.589 | 44.8 mol m-2 |
| calcite 50 t/ha, upper set | 5.546 | 8.131 | 8.125 | 6.281 | 46.4 mol m-2 |

At the forcing's 404 ppm and the campaign's 1:2.5 extract the calcite paste
reads 8.13, whichever set. The 5 t/ha run dissolved all but 0.40 of its
5 mol m-2 within seventeen days, into the exchanger's acidity, as the
50 t/ha run dissolved 4.6 of its 50; what is left cannot saturate the
paste, which reads 7.81. No run raised a warning, the lower set's calcite
run included, whose paste stood 2.6 above its pore water on the first day.
A three-day 2-D run of the upper set's 5 t/ha arm (run 9, twelve columns)
read 8.03 in the eight margin columns that received calcite and left the
four main columns at 6.10.

## What it changes downstream

- The seed. The initial pore-water state is the one whose paste reads the
  observed bulk pH, so it is now inverted through the open paste. Every
  trajectory therefore differs from t's from the first step, the control
  included, by the small difference between the two closures at the seed:
  at alpha 3.4 the control's seed of 1 January 2016 has a pore pH of 5.621
  where the closed paste gave 5.628, both reading 6.100.
- The soil biology. The pH response functions of decomposition,
  nitrification and denitrification read the bulk pH, so in a plot that
  holds solid carbonate they now see about 8.1 where they saw 7.8 to 8.5.
- The 35-year basalt runs. Once they hold nesquehonite their paste reads
  about 9.3 where r read 9.4 to 10.0, and it takes about twice as much of
  the pool to saturate an open paste (3.9 against 1.8 mol m-2 in the unit
  test's sample), because the bicarbonate it forms draws CO2 from the air.

## Limits

- Equilibrium is the state the reading tends to, not the state the
  protocol guarantees. The field's protocol is ten seconds of vortexing
  and an immediate reading. Neither the exchange of CO2 with the air nor
  the dissolution of calcite is complete in that time, so for a
  carbonate-bearing sample the operator gives the end state and the
  measured value could lie between the diluted pore water and it. No
  Kantola plot holding free carbonate has a measured pH to test this
  against. For a soil without solid carbonate the choice is worth about
  0.01 or less in the August mean and up to 0.1 on single days.
- The laboratory air is the forcing's atmosphere of the day. In the
  35-year runs it follows the scenario, so a calcite-saturated paste reads
  about 0.04 (SSP1-2.6) to 0.10 (SSP5-8.5) lower in 2050 than in 2016.
  Indoor air is often richer in CO2 than outdoor air, which would lower a
  calcite reading further by two thirds of the log of the ratio.
- The model's air concentration is ppm/22.41 at every temperature, the
  molar volume at 0 C. At 25 C that is 9 % above the true concentration,
  so the paste holds 9 % more CO2 than a rigorous Henry calculation gives
  (15.0 against 13.7 umol/L at 404 ppm). It is kept because it is what the
  model calls equilibrium with the atmosphere everywhere else; it is worth
  0.025 in a calcite paste and nothing in an acid one.
- With the legacy `TC_GT` scheme a sample holding Mg carbonate beside an
  acid exchanger still stops with `ERW:noBulkPHRoot`, as in r. The model
  holds that pool only where the exchanger has no acidity left.
- Silicates stay inert in the paste, as before.

## Validation

On Bouchet, through the `Base1Dcheck` job of the all-12 controller (the
full audit, the pre-application identity and the smoke runs) before any
production run (`RUNBOOK_2026-09-24.md`, section 0u).

Locally, as a check of the edits only (MATLAB R2026a, one process at a
time): `checkcode` on the operator (no new message); the 98 tests of the
audit list, with the driver-equivalence test on a baseline built for it,
none failing; the fourteen tests that read an edited file again on the
final tree; the saved-state comparison and the six twenty-day runs above;
a three-day 2-D run of the upper set's 5 t/ha arm; and a dry run of both
controllers against mocked Slurm commands. The paste change reaches every
run, which is why the whole list was run here; the authoritative run is
Bouchet's. No five-year run was made locally.

`test_pH_conversion` now checks that the paste holds exactly the air's
dissolved CO2; that doubling the field DIC leaves the reading unchanged bit
for bit and moves only the reported air exchange, by the transferred DIC;
that more CO2 in the air lowers an acid paste slightly and a calcite paste
by two thirds of log10(2); that an excess of calcite leaves the paste at
saturation between 7.9 and 8.5 on one root, in both exchange schemes, with
the Ca it dissolved in the conserved inventory and the CO2 its exchange
acidity released given up to the air; that a nesquehonite paste takes
carbon from the air; that the air exchange closes the paste's DIC with the
transferred DIC and both solids; that a non-scalar air CO2 is refused; and
that both drivers hand the day's atmospheric CO2 to the operator. The
checks of q and r that do not concern the carbon balance are kept.

## The campaign's alpha and the ensemble around it

`runtime/Base_1D/Campaign_Parameters.m` carries alpha 3.4, formed as the
sweep forms it (the last point of `linspace(1.3,3.4,5)`), with the
exchange constants slid by the campaign's rule so that the seed keeps its
measured 80 % base saturation: log KH 5.136, 4.036, 8.943, 9.259, 14.179
for Na, K, Ca, Mg, Al. It is the Appelo (1994) proton-loading slope and
SCEPTER's default; the 1.3 it replaces was fitted on 26 September, before
the leaching closure and the amphibole correction. KNmaxd 0.0048 and free
drainage are unchanged. Alpha sets how far the bulk pH moves for a given
change in base saturation, (0.434/f_H + alpha)/100 pH units per point:
0.035 at alpha 1.3 and 0.056 at 3.4 for the seed's f_H of 0.20.

The uncertainty ensemble (`runtime/Uncertainty_1D/README.md`) keeps its
sixteen parameters and its ranges, with alpha at 1.3 in the lower set and 6
in the upper set, the ends of the sweep grid. No range is published for
alpha, so this is an exploratory bracket, as the KNmaxd range is. The
dissolution factor stays at 0.8 and 1.2 around 1.

Both controllers can now submit the 1-D runs alone (`TC_ALL12_DIMS=1d`,
`TC_UC_DIMS=1d`), so the five-year 1-D results can be read before the 2-D
and 35-year campaigns are released; `readout_1d.m` at the package root
prints them.

## Sources

* Appelo, C.A.J., 1994. Cation and proton exchange, pH variations, and
  carbonate reactions in a freshening aquifer. Water Resources Research
  30, 2793-2805. doi:10.1029/94WR01048
* Beerling, D.J. et al., 2023 preprint. Enhanced weathering in the U.S.
  Corn Belt delivers carbon removal with agronomic benefits.
  arXiv:2307.05343 (soil pH method).
* ISO 10390:2021. Soil, treated biowaste and sludge: determination of pH.
* Kanzaki, Y., Chiaravalloti, I., Zhang, S., Planavsky, N.J., Reinhard,
  C.T., 2024. In silico calculation of soil pH by SCEPTER v1.0.
  Geoscientific Model Development 17, 4515-4532.
  doi:10.5194/gmd-17-4515-2024, with the SCEPTER 1.0.2 source
  (`get_soilpH_time.py`, `scepter.f90`).
* Plummer and Busenberg 1982 and Harrison et al. 2013, 2019, as cited in
  `PASTE_CALCITE_2026-10-01.md` and `PASTE_MG_CARBONATE_2026-10-01.md`.
