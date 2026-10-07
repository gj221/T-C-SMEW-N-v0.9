# Post-legacy fixes, 25 September 2026

Source: `.codex_phcharge_20260923`, the package behind
`tcsmewn_v09_kantola_knmaxd_20260925` (runs 10-12, array 27485170). That
package is unchanged. Every fix below was written test first in this copy,
which is packaged as `tcsmewn_v09_20260926j`. The three-set uncertainty
ensemble and the SSP controller's pin path (last section) follow in
`tcsmewn_v09_20260926k`, whose model is byte for byte j's.

## Changes that move results

**Exchange seed.** logKH now moves with alpha by `charge*(alpha-1.3)*0.20`
instead of `charge*(alpha-1.3)`. Selectivity is `10^(-logKH + alpha*charge*f_H)`,
so the new shift holds it at the seed's `f_H = 0.20` and the pH 6.1 inversion
keeps the measured 80 % base saturation at any alpha. At alpha 2.35 the old
shift seeded 68.7 % base saturation with pore Ca and the residual anion about
five times higher. The rule is set in `Campaign_Parameters.m` and
used by the uncertainty task, its summary and the joint sweep cases.
Test: `test_kantola_seed_charge` now also seeds at the campaign alpha.

**Harvest reserve.** At harvest the remaining nutrient reserve is exported off
the field through `AgrHarNut`, as in legacy. The post-legacy `case -2` block in
`Plant_Exports.m` sent it to surface litter, which put the unused part of the
legacy 20 g N m^-2 (maize) and 15 g N m^-2 (soybean) planting reserve into the
soil every year. Tests: `test_crop_harvest_reserve_ledger`,
`test_crop_harvest_scale_reserve_export`.

**N fixation cost.** `Root_Exudation.m` (both branches) now uses the Fisher et
al. (2010) cost the soil side already used (`s = -6.25`, `c = 25.15`). With
`s = -30` on the plant side the soil credited 4.8 times the N the plant paid
for, in the soybean years. Test: `test_biological_fixation_carbon_allocation`.

**Calcite solubility.** `K_CaCO3` follows soil temperature with the Plummer and
Busenberg (1982) polynomial (-8.48 at 25 C) in the seed and in every ERW step,
replacing the fixed `10^-8.35`. Each step now records `Omega_CaCO3`. Tests:
`test_erw_omega_activity_basis`, `test_calcite_event_smoke`.

**Plant tissue nutrients.** Plant nutrients follow legacy T&C again. The
September design kept the reserve at or above zero. Growth drew on the
reserve alone (`Navail = Nreserve`) through one scale,
`Plant_Reserve_Growth_Limit`, that cut growth and translocation together,
and guards rejected any negative reserve. In legacy a negative reserve is tissue
diluted below its nominal content, `rXc = 1 + Xreserve/TXns`, and growth
may dilute the non-structural tissue to 0.65 of nominal:
`Xavail = Xreserve + 0.35*TXns`, clipped at zero, with `rCa_min` and
`rMg_min` (0.65) setting the Ca and Mg floors. Each nutrient in turn cuts
NPP to what its availability supports and the rejected carbon is respired
(`B(8)`). Translocation stops once an availability is exhausted. As
tissue N falls, photosynthesis and maintenance respiration fall with `rNc`,
which is the N dilution response of field crops (Greenwood et al. 1990;
Plenet and Lemaire 1999 for maize). The growth step in
`VEGETATION_DYNAMIC` is legacy again line for line but one: on a day of
negative NPP the reserve pays leaf dark respiration, and that transfer now
replaces the leaf-out translocation to the roots as well (`Tr_r = 0`),
where legacy gave the roots carbon the reserve never paid.

Where the September limiter never bound, the two designs give the same run
bit for bit. The ensemble's all-minimum member, whose soybean used up its
N on 15 September 2016 and stopped at the guards, now dilutes its tissue
instead. The tolerance and absorption of package g went with the
guards.

At harvest the organs leave at the concentration the plant has that day.
A deficit is diluted tissue, so it goes with the grain and the residue and
nothing is left for the reserve export; a surplus stays in the reserve and
leaves through `AgrHarNut`, as before. `Nutrients_Available` and the
harvest reject a deficit larger than the tissue nutrient, the one state in
which the plant would hold negative N, and the N audit checks plant N,
tissue plus reserve, where it checked the reserve alone. Tests:
`test_plant_tissue_dilution`, `test_crop_harvest_reserve_ledger`,
`test_kantola_n_magnitude_audit`, and `test_checkpoint_fingerprint_dependencies`,
which now also requires every source the checkpoint lists to exist.

## Changes that do not move the current runs

**Na and residual anion closure.** Both drivers now run the legacy steady
closure (`keyword_add = 1`), which production already used. The 1-D frame no
longer switches it off for Al-pinned runs, and `test_2d_driver_equivalence`
checks that its baseline used the workbook's closure. With that, the 1-D and
2-D drivers agree bit for bit on the production branch.

**Albite and anorthite constants.** These are now the llnl.dat Al+++ reaction
constants rooted per proton (log K 2.7645 and 26.578). The previous values
belonged to the legacy kaolinite reactions. Albite still does not dissolve,
because ERW receives the whole T&C Si pool as dissolved Si (legacy).
Test: `test_erw_mineral_kinetics`.

**Daily dates in post-processing.** Row j of the daily state holds the BIOGEO
day from 01:00 of `Date(1)+(j-3)`. `Plot_1D_Then_All_12`,
`Summarize_Uncertainty` and `Summarize_Sweep` now date
rows that way (they ran two, two and one day late), and daily WFPS and drainage
average the same hours. Model N2O is therefore scored against the right
chamber days. The first retained row is the lone first hour, dated the day
before the start, and the last complete day of 2020 is 29 December. The
sweep's pre-period gate and pH mean, its 2020 carbon year and the
uncertainty summary's calendar years count complete days only. Tests:
`test_plot_configuration_figures`, `test_bouchet_1d_joint_sweep_configuration`,
`test_uncertainty_campaign`.

**2-D lateral transport (runs 1-9 only).** Lateral conductivity is now the
horizontal `K*aR` (the adapter passes `aR`; it was dropped) and frozen-pore
impedance is the legacy `10^(-5*F_ice)` the vertical solver uses (it was
Omega 6 laterally). Soil air displaced by lateral water is no longer booked by
a piston term in `ERW_Lateral_Transport`, because the next ERW step's `ADV`
books the same change; the lateral step now closes inorganic carbon, and its
unused `CO2_atm` input is gone. Dissolved DIC still moves with the water and
soil-gas CO2 still diffuses between cells. The atmosphere residual that
measured the piston is gone from the driver and from
`Check_Lateral_Conservation`, where a zero exchange would have
failed every 2-D run on accumulated round-off; DIC in the lateral mass
residual carries the carbon check. Tests: `test_lateral_energy_erw`,
`test_lateral_routing`, `test_model_adapter_contracts`,
`test_pool_overdraw_guards`, `test_2d_amendment_conservation`,
`test_bouchet_2d_parallel_configuration`.

**Campaign names.** Campaign folders and files carry broad names; only the
location-specific inputs keep the site or study in theirs. `Base_1D` and
`Base_2D` hold the five-year validation of runs 10 to 12 and 1 to 9,
`SSP_1D` and `SSP_2D` the projections, `Calibration_Sweep_1D` the joint
sweep and `Uncertainty_1D` and `Uncertainty_2D` the ensemble, with the calibrated set in
`Base_1D/Campaign_Parameters.m`. Results go under the same names in
`v2/SIM/Subfield_12`; projection runs and ensembles add horizon and scenario
(`SSP_1D_35yr_ssp245`, `Uncertainty_1D_35yr_ssp245`). Slurm job names follow
(`Sweep*`, `Base1D`, `Base2D`, `SSP1D`, `SSP2D`, `UC1D`, `UC2D`). The
superseded campaigns (the pre-CEC13 sweep, the schedule-fix reruns, the
DenitPH7p8 2-D runs, the moisture sweep and the Appelo literature pair) are
gone with their tests, manifests and the old all-12 controller, which the
validation fingerprint no longer hashes; `phobs.m` moved to `Base_1D`.
`test_bouchet_2d_parallel_configuration`, which both 2-D campaigns run, now
checks `Base_2D`. Test: `test_bouchet_linux_portability` rejects site,
cluster or parameter values in any name under `runtime/`.

**5 t/ha calcite on the validation forcing.** `submit_all_12.sh` takes
`TC_CALCITE_CASE=8`, as the SSP controller does, and runs the calcite
treatments (1-D run 12, 2-D runs 3, 6 and 9) at the agronomic 5 t/ha rate of
workbook case 8 as a campaign of their own (`Base_1D_case8`,
`Base_2D_case8`), beside the 50 t/ha one. The SSP 2-D task refused a
calcite-case campaign's validation element (array element 1, run 1) before
it reached validation, so no 5 t/ha projection could start; both 2-D tasks
now check the run after validation. Test: `test_calcite_case_campaigns`.

**Pools emptied to rounding.** Once its soybean could dilute instead of
stopping, the all-minimum member emptied the soil of mineral N in
mid-September 2016, and a pool that the day's limited sinks drain to zero
can land a few ulps below it. Strict sign checks then stopped the run
twice: the weathering step's NO3 and NH4, which follow T&C's daily states
through the sub-steps (-3.5e-15), and the next day's soil check after
uptake (-8.7e-27 g g^-1 soil), fed by the unchecked daily update
`P = Ptm1 + dP` on a day when immobilization took all the mineral N. The
same arithmetic sat in the uptake cap, `uptake*(headroom/uptake)`, which
can round one ulp above the headroom, and in the hourly lateral exchange of
the 2-D runs. Each of these now ends an emptied pool at zero, within the
band of 16 eps (per gram of soil in T&C) that the donor limits already
treat as nothing, and the uptake cap stops at the headroom exactly. A
larger deficit, a nonfinite value or a negative T&C state coming into the
weathering step is still an error, now raised on the day it happens. Every
run that passed before is unchanged bit for bit, since each change acts
only where a run would have stopped. Tests: `test_erw_input_guards`,
`test_biogeo_pool_exhaustion`.

**Baseline script.** `regen_equiv_1d.sbatch` builds the equivalence
baseline under the campaign identity: a `cpumodel:9575f` node, R2025a, one
computational thread and the thread variables, as the sweep's own baseline
job does. It carried none of those, so the baselines of the 26 September
roots were built on whatever node Slurm chose, with MATLAB free to thread.
`test_2d_driver_equivalence` holds the 2-D driver to the baseline exactly
for water, energy and vegetation state, and against such a baseline it
fails on round-off from the third hour, which is what stopped both all-12
validations of 27 September (jobs 27638811 and its 5 t/ha twin) before
production. Test: `test_all_12_configuration` pins the script.

## Kept as is

The bulk-pH extract leaves undissolved calcite out, following SCEPTER. The
denitrification plateau stays at pH 7.5 (set on 22 September with its
rationale in `Biogeo_Denitrification_PH_Response.m`); the README and audit
now say 7.5, and `test_nitrogen_process_shapes` pins it. ERW still receives the
whole T&C Si pool, and T&C still drains and moves Ca and Mg at 0.4 of pool 60
and 68, both as in legacy. The Blue Ridge mineral formulas stay the
preconfigured legacy stoichiometries.

**Microbial maintenance rates.** mrb and mrf stay at the legacy T&C-SMEW
code's 0.010 and 0.005 g C g^-1 C d^-1, under which `Data_Kantola_BG.mat` was
spun up; the SI's Table S2 lists 0.005 and 0.002 as its originals. In
`BIOGEOCHEMISTRY_DYNAMIC3.m` these rates also set microbial mortality and
enzyme production, so steady-state biomass scales as 1/mr. Starting the
spin-up at the SI values would drive its bacteria toward twice and its fungi
toward 2.5 times their biomass, taking up roughly 4.8 g N m^-2 in the scored
years. The uncertainty ensemble's lower and upper sets still move both to
the Table S2 minima and maxima (0.0012 and 0.03). Test:
`test_uncertainty_campaign`, which pins both sets to Table S2 and builds
their parameters the way the soil model does.

## Calibrated set

Alpha 1.3, KNmaxd 0.0048 and free drainage, combination
`a1p3_kn0p0048_kbfree` of the refit sweep of 26 September (root
`tcsmewn_v09_kantola_postlegacy_20260925`, arrays 27559013 to 27559017). That
sweep ran with the legacy maintenance rates but under the September
plant-reserve design. Both two-objective fronts sit at alpha 1.3, the floor of the grid,
where pH changes little with alpha. KNmaxd 0.0032 ties 0.0048 on the N2O
log-ratio, but only 0.0048 reproduces the total chamber N2O (0.97 of observed
against 0.68), and free drainage scores best on N2O and NO3.
`Campaign_Parameters.m` carries the set, and the uncertainty ensemble's upper
set moves alpha to the grid point 2.35.

Under legacy tissue dilution runs 10 to 12 no longer match the sweep's cases
1 to 3 bit for bit (ensemble base member, 27 September): the September
limiter did bind in them. The shifts are small against the margins that
chose the set. Daily bulk pH moves by at most 0.008, 0.012 and 0.052
(control, basalt, calcite), total N2O by about 1% and 2017 NO3 leaching by
0.3%, and the one-year gate's pre-period residual by 0.001 pH. KNmaxd and
drainage were chosen on total chamber N2O (0.97 of observed at 0.0048
against 0.68 at 0.0032) and on NO3, where free drainage scored best. Alpha
sits at the floor of a pH gradient, about 0.003 in the seasonal and 0.01 in
the treatment-difference RMSE per grid step and the same way in every
combination, which comes from the exchange chemistry that the plant change
shifts in common. The sweep was therefore not rerun: its pilot and
production were cancelled on 27 September, and the ensemble's alpha range
carries the alpha uncertainty.

## Uncertainty ensemble and SSP pin path (package k)

**Three parameter sets.** On 27 September the user replaced the 36-member
one-at-a-time ensemble with the legacy analysis's own design: the original
values, every parameter at the lower end of its range and every parameter
at the upper end. The fourteen legacy parameters (Table S2's twelve
literals with the dissolution factor and pore-saturation exponent) gain
this model's KNmaxd (0.002 to 0.01 d^-1) and alpha (the calibrated 1.3 to the sweep
point 2.35, constants slid by the campaign rule). The sets cover every
treatment in 1-D and 2-D, the 5 t/ha calcite arm included, over five years
on the legacy forcing and 35 on each SSP. The original is not rerun: it is the Base and SSP campaign result of
root j, and the controller refuses to submit unless this checkout's model
code and inputs are byte for byte j's (`same_model_code.sh`). Every 2-D
element runs at its original's worker count (6 at five years, 12 at 35).
A saved set is accepted only with its full provenance, checked per
column. The ensemble runs on the `day` partition, whose per-user limits
(1500 CPUs, 20 TB) are far wider than `week`'s (96 CPUs, 1.5 TB, where
the Base and SSP campaigns run and a 384 GiB element runs four at a
time). A run longer than a day requeues itself when Slurm signals the
limit half an hour ahead and continues from its checkpoint under the
same job id, at most five times. `TC_UC_WITH_ORIGINAL=1` runs the original
set there as well, for when the original's own campaign is held back by
the week limits; an original run in both roots must match bit for bit. The 36-member summariser is removed; the summary and figures for
three sets come with the post-processing package. Tests:
`test_uncertainty_campaign`.

**SSP pin list.** `submit_ssp.sh` looked for its pin list in
`runtime/SSP_1D`, where it has never been, so the controller refused every
campaign until the list was copied there. It now reads
`runtime/Base_1D/ssp_pins.sha256`. Test: `test_ssp_controller_files`, which
resolves every file the controller requires.
