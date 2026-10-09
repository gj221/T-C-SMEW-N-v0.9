# T&C-SMEW-N v0.9

Package ab is package aa with the five-year field campaign
(`runtime/Field_2D`) held to z's fast ERW kernels. Its batch file sets
`TC_ERW_KERNEL=fast` for every element, so nothing needs exporting. Its
task refuses any other kernel and checks that every result records the
fast one (`run_provenance.erw_kernel`), and its validation saves the
one-day grid run and checks it the same way. The fast kernels compute the
original code's numbers bit for bit; nothing else changes.
`RUNBOOK_2026-09-24.md`, section 0ab. The paragraph below describes aa,
which ab contains.

Package aa is package z with the campaign side of the 108 35-year runs
on the 10 m Kantola grids (`FIELD_PROJECTIONS_2026-10-08.md`):
`runtime/Field_SSP_2D` runs the central SSP runs on a field grid and
`runtime/Field_Uncertainty_2D` the lower and upper sets of the
uncertainty ensemble on it, both with Field_2D's gates, `TC_FIELD_WORKERS`
workers, two-hourly checkpoints, a memory floor from y's daily arrays,
z's ERW kernels (`TC_ERW_KERNEL`, checked and recorded per run) and
24-hour pieces on the day partition that requeue themselves. The central
controller adds runs beside running ones, so the pilot (run 2) sizes the
rest while it runs. `Summarize_Uncertainty.m` at the package root gathers
the three sets of an arm from y's daily stores and ERW day summaries (full
records for older results), field-averaged as Domain averages. No
equation, parameter or solver setting changes: `Run_One_Subfield` only
keeps a result or checkpoint name within 255 characters
(`Bounded_File_Tag`), since every field-grid ensemble name was 258 to 272
characters; a name that fits is unchanged. `RUNBOOK_2026-09-24.md`,
section 0aa, is the procedure. The paragraph below describes z, which aa
contains.

Package z is package y with restructured kernels for the ERW chemistry
step (`ERW_KERNELS_2026-10-08.md`). `TC_ERW_KERNEL` selects them: unset or
`matlab` runs the original code, `fast` restructured MATLAB kernels that
do the same floating-point operations in the same order, and `mex` the
same kernels compiled with MATLAB Coder in the deployed root
(`build_erw_mex.sbatch`). All three compute every value bit for bit
alike. On Bouchet, against package x's code, `fast` ran the ten-minute
step about 11 times faster and a one-column year 4.8 times faster; `mex`
ran the step about 21 times faster. MATLAB runs no MEX function on a
thread-pool worker, and the squares of a 2-D run share a thread pool, so
the field runs use `fast`; `mex` serves one-column and serial runs, and
`Run_One_Subfield` refuses it with a pool. Each result records its kernel
in `run_provenance.erw_kernel`. No equation, parameter or solver setting
changes, and a run that sets nothing runs y's code. `compare_erw_kernels.sbatch`
runs the same runs under each kernel, and in x, and checks every saved
value; `compare_x_y.sbatch` can run y's side with `CMP_Y_KERNEL=fast`
(`RUNBOOK_2026-09-24.md`, section 0z). The paragraph below describes y,
which z contains.

Package y is package x with the records of a multi-square run kept daily
and its squares computed in parallel, for the 35-year runs on the 10 m
Kantola grids (`FIELD_RECORDS_DAILY_2026-10-08.md`). No equation,
parameter or solver setting changes, and every computed value is x's, bit
for bit. A one-column run keeps x's full records, so the 1-D campaigns
and tools are unchanged. A run of several squares keeps the last 50 hours
of each hourly array and one day of ten-minute chemistry, and saves daily
stores per square (the daily arrays the readers use, daily means and
totals of the hourly water fluxes, one summary per ERW day), hourly field
means of the tower fluxes and its final state, instead of the whole
workspace (`TC_FULL_RECORDS=1` keeps the full records). The daily
vegetation and biogeochemistry, the hourly vertical physics and the daily
ERW chemistry run square by square on the thread pool; lateral transfer
stays serial. The field maps and the bulk-pH summary read the daily
stores. `Load_Field_DEM` loses the backslashes that made x fail the
portability gate of the field campaign. `compare_x_y.sbatch` runs x and y
side by side on the plot (40 days, and 400 days across two applications)
and on the basalt grid and checks every output (`RUNBOOK_2026-09-24.md`,
section 0y). The paragraph below describes w;
package x added v1.0's DEM topography, lateral transfer and field grids to
it (`FIELD_GRIDS_2026-10-07.md`, `runtime/Field_2D/README.md`).

Package w is package v with one guard relaxed (`ROUNDOFF_NEGATIVE_2026-10-05.md`):
`ERW_Ionic_Strength` counts a concentration within 1e-10 umol-conv/L below
zero as zero instead of stopping the run. A 35-year 2-D uncertainty run
stopped on -4.5e-16, a cancellation residue. Every run that got past the
check gives the same numbers in v and w. One audit test file differs from
v's (`test_erw_activity_coefficient`), and the comparison utility at the
package root (`compare_u_v.m`) now also compares two results that both
keep daily particle rows, as v and w do, and compares the fsolve options
by their settings. `RUNBOOK_2026-09-24.md`, section 0w, is the procedure.
The paragraph below describes v, which w contains.

Package v is package u with two changes to how a run holds and saves its
state, and none to the model (`PARTICLE_RECORDS_DAILY_2026-10-04.md`).
The particle-class histories of an amended column (`psd`, `d`, `lamb`,
`SSA`, `delta_d`) keep their last three days whole and, of every earlier
day, its last row: the model reads only the previous row and each
application's first row, so every computed value is u's, and the results
hold those histories at the initial rows and the last row of every ERW
day (`weath.particle_rows`) instead of every ten-minute row. Held whole
they took 16.8 GiB per amended column over 35 years, about half of a
whole-field basalt run's memory. The SSP 2-D campaign checkpoints every
two hours instead of every thirty minutes, and the 35-year 2-D runs ask
256 G instead of 384 G, so six run together on the week partition. One
audit test is rewritten (`test_erw_particle_history`) and one pinned line
changes (`test_uncertainty_campaign`). `RUNBOOK_2026-09-24.md`, section
0v, is the procedure for this root. The paragraph below describes u,
which v contains.

Package u is package t with three changes. The soil-water paste of the
bulk-pH operator stands open to the laboratory air
(`PASTE_OPEN_TO_AIR_2026-10-02.md`): its dissolved CO2 is the air's at
25 C instead of the CO2 the field pore water carried, and the plot's solid
carbonates still dissolve in it, so a soil holding free calcite reads the
calcite-air equilibrium of its paste, about 8.1, and a soil without solid
carbonate reads within about 0.01 of the closed paste in its August mean
(up to 0.1 apart on single days). The
campaign's Appelo alpha is 3.4, the Appelo (1994) slope SCEPTER carries as
its default, in place of the 1.3 fitted on 26 September, and the lower and
upper sets of the uncertainty ensemble carry alpha 1.3 and 6 around it
(`runtime/Base_1D/Campaign_Parameters.m`,
`runtime/Uncertainty_1D/README.md`). Both campaign controllers can submit
their 1-D runs first (`TC_ALL12_DIMS=1d`, `TC_UC_DIMS=1d`), and
`readout_1d.m` at the package root prints the August read-out of those
runs. Every trajectory differs from t's: the seed is inverted through the
open paste, and the soil biology reads the open reading. Two audit test
files differ from t's (`test_pH_conversion`, which pins the open paste, and
`test_uncertainty_campaign`), with the all-12 configuration test.
`RUNBOOK_2026-09-24.md`, section 0u, is the procedure for this root. The
paragraph below describes t, which u contains.

Package t is package r with one change to the model and the 1-D sweep
returned to dissolution factor 1 (`AMPHIBOLE_KINETICS_2026-10-02.md`). The
basalt's amphibole (ferroactinolite, 11.6 % of the rock) ran on the
tremolite row of Palandri and Kharaka (2004), which is taken from Sverdrup
(1990) rather than fitted to tremolite experiments and is about 60 times
the measured rate at pH 6. That one mineral supplied 11.5 of the 12.5 eq
m-2 of silicate alkalinity by August 2020 and used up the exchanger's
10.3 eq m-2 of acidity in 2019. It now carries the amphibole rate law of
Hermanska et al. (2022). Re-evaluated offline from saved runs, the rock
then releases about 0.7, 1.5, 2.5 and 3.5 eq m-2 by August 2017 to 2020 at
dissolution factor 1, against a field exchanger gain of 1.2, 2.3, 3.0 and
3.0. Control and calcite runs hold no amphibole and are r's bit for bit.
The sweep is package o's again: seven alphas, four KNmaxd, four bottom
boundaries, 112 combinations, 336 trajectories, with `DissolutionFactor` 1
in every case. Two audit test files differ from r's
(`test_erw_mineral_kinetics`, which pins the new constants, and
`test_bouchet_1d_joint_sweep_configuration`). Packages p, q, r and s are
superseded for anything that holds basalt. The paragraph below describes
r, which t contains.

Package r is package q with the plot's solid Mg carbonate (the nesquehonite
pool) dissolving in the soil-water paste beside its calcite
(`PASTE_MG_CARBONATE_2026-10-01.md`). The pool is zero in every five-year
run, so the sweep, the controls, the calcite runs and the five-year basalt
runs are q's bit for bit; the 35-year basalt runs read a higher bulk pH
once the pool has formed. Two audit test files differ from q's:
`test_pH_conversion` (new checks) and
`test_bouchet_1d_joint_sweep_configuration` (the wording of one message).
The notes now say "dissolution factor" where p and q said "release
factor"; the option, the workbook column and every number are unchanged.
The paragraph below describes q, which r contains.

Package q is package p with one change to the model: the bulk-pH operator
lets the plot's solid calcite dissolve in the soil-water paste, as it does
in a real pH measurement, instead of freezing it with the silicates
(`PASTE_CALCITE_2026-10-01.md`). The pore-water chemistry is untouched; a
plot holding free calcite now reads 7.8 to 8.3 where it read 7.4 to 7.8,
and plots without solid calcite are unchanged bit for bit. Two audit tests
differ from o's (`test_bouchet_1d_joint_sweep_configuration`,
`test_pH_conversion`). The paragraph below describes p, which q contains.

Package p is package o with the dissolution factor (Run_One_Subfield's
`DissolutionFactor`, the workbook `diss_f`) added as the fourth dimension of
the 1-D calibration sweep and the other three grids trimmed to keep the
sweep's cost: 150 combinations, 450 trajectories
(`DISSOLUTION_FACTOR_SWEEP_2026-10-01.md`). The model, every campaign wrapper
and 93 of the 94 audit tests are o's unchanged; the inherited export
operator is kept as inherited by decision. Case 1 to 3 of the sweep remain
the reference tuple, so the one-year gate and the all-12 campaign label
(`a1p3_kn0p0048_kbfree_fd1`) carry the same physics as o. The paragraphs
below describe o, which p contains.

Package o is package n with its Bouchet validation repaired and made
diagnosable (`VALIDATION_FIXES_2026-10-01.md`): the equivalence baseline
now fingerprints only the four input files it reads, so the scenario
forcings copied into a root cannot stale it; a refused baseline names the
files that moved; both 1-D validators print a timestamped banner at
every stage; and the one audit test n's salt closure broke
(`test_erw_lateral_cation_daily_delta`, whose last-day expectation
predated a boundary Ca and Mg supply) now accounts for it. The model is
n's, unchanged, and the audit suite and smoke runs were executed on it
before packaging. Everything below describes
that model.

Package n lets the 0-30 cm box export what its drainage water carries.
Dissolved Ca and Mg now leach and are taken up at their concentration
(`aCa = aMg = 1`, as NO3 and K already were; the legacy 0.4 belonged to a
lumped pool), and the leached background anion returns as the salt it left
as instead of a bare anion. In package m the leachate carried 2.5 times
less Ca+Mg charge than anion charge, 84 % of the alkalinity the rock
released stayed in the box, the exchanger saturated by 2019 and the basalt
pH ran away; the field reached pH 6.9 at 83 to 90 % base saturation. The
sweep now scores that base saturation beside pH, over August, when the
soil was sampled, with the alpha grid extended to 6; the calibrated alpha
is to be re-chosen from it before production. Evidence, equations and
checks are in `LEACHING_CLOSURE_REVIEW_2026-10-01.md`. The Mg carbonate
sink of m (`MG_CARBONATE_REVIEW_2026-09-30.md`) and the particle-history
storage of l (`PARTICLE_HISTORY_REVIEW_2026-09-28.md`) are kept. Every run
changes, the control included, so this build needs its own Bouchet R2025a
equivalence baseline, validation and sweep before production; use the
runbook here.

## Repository layout

```
v2/run/          launchers: ERW_run, Run_Matrix_12, Run_One_Subfield,
                 Run_All_12_Parallel, Run_Subfield_Worker,
                 Terminate_Process_Tree
v2/plotting/     figure scripts: Plot_1D_Then_All_12, Plot_Grid_Setup,
                 Plot_Biogeo_Pools
v2/figures/      every PNG and FIG the plotting scripts write
v2/T&C_Code/     model source, Inputs, Parameters and Tests
v2/SIM/          simulation results (not tracked)
v2/logs/         worker logs and status files (not tracked)
v2/Experimental observations/   digitised field data used by the benchmark
```

## Site scope

This v0.9 checkout targets the Energy Farm (Kantola) simulations only.
Bruce Field support is separate from this campaign and lives in the full
v0.9 checkout; neither the Bruce wrapper nor its generated forcing files
are carried here.

Every plotting entry point is a runnable script: `cd` into `v2/plotting` and
run it by name, with no arguments. Each one puts the folders it needs on the
path, so it works from any working directory.

## Energy Farm (Kantola) simulations

The twelve simulations are defined in `v2/run/Run_Matrix_12.m`. Each group
contains control, basalt and calcite treatments.

| Runs | Scenario group | Domain and application |
|---:|---|---|
| 1 to 3 | whole field | 12 cells, amendment on the main field and 5 m margin |
| 4 to 6 | margin with crop | 12 cells, amendment on the cropped 5 m margin |
| 7 to 9 | margin with prairie | 12 cells, amendment on the prairie 5 m margin |
| 10 to 12 | vertical 1-D | one 3600 m2 column, amendment across the plot, no lateral transport |

The retained 12-cell geometry has four 25 m by 25 m field cells, four
50 m by 5 m side margins and four 5 m by 5 m corners. The prairie uses
big bluestem as the high C4 PFT and purple coneflower as the low C3 PFT.

All runs use the control Kantola T&C parameters, forcing and initial
state. The active biogeochemical depth is 300 mm in every run, matching the
workbook value and the spin-up that produced Data_Kantola_BG.mat. Treatment rows only supply
the amendment properties. Untreated cells use the control chemistry row.
The configuration workbook contains only the five Kantola cases used by
the single-column comparison and the twelve-run matrix.

Aqueous ERW speciation and mineral reactions use pore-water pH. T&C microbial
rates, including nitrification and denitrification, use the reconstructed 1:2.5
DI-water bulk-extract pH, following the Beerling et al. preprint soil-pH method.
The Kantola seed uses the measured pretreatment bulk
pH as an inversion boundary while retaining the prescribed exchange inventory.
There is no additive pH offset.

The field charge balance includes dissolved Al hydroxide ions. `Alk_tot` is the
conserved strong-ion inventory; pore-water `Alk` is carbonate alkalinity.
Primary calcite contributes to the carbonate, not silicate, alkalinity source
diagnostic.

ERW cation partitioning follows the selected pH exchange scheme. `TC_GT` keeps
the legacy T&C fixed-pool rates and disables only its kinetic solution-exchange
rates. `appelo` disables both legacy rate pairs, leaving the conserved
SCEPTER/Appelo equilibrium as the only K, Ca and Mg partition operator.
Standalone T&C runs keep all legacy rates.

### Nitrogen transformation contract

The active ammonium and nitrate pools are `B(31)` and `B(32)`. They are
converted from g N m^-2 to g N per g dry soil inside
`BIOGEOCHEMISTRY_DYNAMIC3` and converted back to g N m^-2 d^-1 for output.
The corrected gross denitrification rate is:

```text
DOC_v = (B(12) + B(13)) * rho_soil                         [g C m^-3]
f_DOC = DOC_v / (DOC_v + Kc), Kc = 17.5 g C m^-3
DENIT = k_den * f_T * f_water * f_den_pH * f_DOC * B(32)  [g N gsoil^-1 d^-1]
```

This retains the first-order nitrate reservoir loss in Dickinson et al. (2002)
and legacy T&C-SMEW. The Li et al. (1992) nitrate half-saturation parameter is
volumetric, `Kn = 82.5 g N m^-3`, but remains inactive because TC-SMEW does not
carry the denitrifier biomass, yield and maintenance states required by the
DNDC double-Monod rate. The rejected `k*N^2/(N+Kn)` form is not a total
Michaelis-Menten rate.

Gross denitrification uses a clipped linear bulk-pH response from pH 4 to a
plateau at pH 7.5. The shape is adapted from Blanc-Betes et al. (2021), which
used a pH 7.0 breakpoint. 7.5 is the top of the range Olaya-Abril et al.
(2021) tested and inside the pH 7.0 to 8.0 optimum band of Simek and Cooper
(2002); it is not a universal optimum. The N2:N2O pH
factor is Blanc-Betes Eq. 9 with no floor of its own; the published 0.1 minimum
applies to the complete N2:N2O product, as in their Eq. 6. The nitrate,
respiration and WFPS terms are the legacy TC-SMEW DayCent forms and change only
the product split, so the composite ratio is a hybrid rather than the published
Blanc-Betes ratio. N2 plus denitrification-derived N2O
therefore equals the gross nitrate sink.

Ammonia volatilization remains the legacy empirical term
`VOL = kv*B(31)`, with `kv = 8.6e-5 d^-1`. It has no direct pH multiplier or
explicit NH3/TAN equilibrium. pH affects it indirectly because the pH response
of nitrification changes the competing NH4 demand and the NH4 pool available
on later days. Treat `VOL` as an NH4-based legacy loss diagnostic, not a full
surface NH3 transfer model.

Nitrate leakage remains an inventory-limited mass-flow term driven by drainage.
It has no direct pH multiplier. pH affects leakage indirectly by changing
nitrification, denitrification and the nitrate pool available to drain.

The detailed equations, dimensions and citations are in
[`v2/POST_LEGACY_EQUATION_AUDIT.md`](v2/POST_LEGACY_EQUATION_AUDIT.md).

### Phosphorus from weathering

Mineral P released by weathering returns once per ERW window to the T&C
mineral-P state `P(:,43)`, in g P m^-2, through
`ERW_Return_Weathered_Phosphorus`. The 12-cell and vertical 1-D layouts share
that one transfer. The window budget is checked at run time, so a replayed or
inconsistent window raises an error instead of creating P.

No feedstock row used by the twelve runs contains a P-bearing mineral, so
basalt-derived P supply is exactly zero in the current campaign. The Blue Ridge
metabasalt mineralogy follows Lewis et al. (2021),
doi:10.1016/j.apgeochem.2021.105023, which did not detect apatite or any
K-bearing mineral in that rock. Bulk P and K that XRD did not resolve into a
named phase are not released, because release is computed from named mineral
stoichiometry. Do not add an apatite fraction without a feedstock measurement
that assigns P to a phase.

The archived Kantola background file contains areal pool totals and no
depth metadata. It is reused unchanged in every case. This keeps the runs
comparable, and the 300 mm active depth is the depth it was spun up at.

The 12-cell runs exchange subsurface water with a finite-volume
Buckingham-Darcy calculation. Water advects the dissolved pools used by
`Biogeo_Leakage`: DOC, DON, DOP, NH4, NO3, mineral P, solution K, Ca, Mg
and Si. ERW-only Na, Al, residual anion and DIC use the same active-depth
water flux. Soil-gas CO2 diffuses across shared faces. Soil heat moves by
Fourier conduction and with liquid-water sensible heat. The outer plot
boundary is closed laterally. Existing vertical drainage and solute
leaching remain active. Frozen hydraulic conductivity uses the Hansson
et al. (2004) ice-impedance factor. Fluxes blocked by frozen pore space
are excluded when the adaptive lateral time step is calculated.

### Energy Farm Bouchet R2025a runs

`RUNBOOK_2026-09-24.md` holds the Bouchet deployment and the order of the
campaigns. The twelve-run validation goes through
`runtime/Base_1D/submit_all_12.sh`, which needs the R2025a equivalence
baseline in place (`regen_equiv_1d.sbatch`). It submits the 1-D and 2-D
validations (6 CPUs each, 40 GB and 32 GB) and both production jobs held,
then releases them; production waits for both validations, after which
runs 10 to 12 and the 2-D array of runs 1 to 9 become eligible together.
Every job pins `cpumodel:9575f` and records its hardware, and production
needs a fresh R2025a validation stamp, stopping if source, inputs, tests
or the numerical platform change after validation. Calling a campaign's
batch file directly bypasses these checks.

From PowerShell at the repository root:

```powershell
matlab -batch "cd('v2/run'); Run_All_12_Parallel('Runs',10:12,'MaxConcurrent',3,'MaxRuntimeHours',168)"
matlab -batch "cd('v2/run'); Run_All_12_Parallel('Runs',1:9,'MaxConcurrent',2,'ParallelWorkers',6,'MaxRuntimeHours',168)"
matlab -batch "cd('v2/plotting'); Plot_1D_Then_All_12"
```

The plotting command first draws the vertical 1-D weathering and nitrogen
figure (`vertical_1d_weathering_nitrogen.png` and `.fig`). It then writes one
PNG and FIG pair per summary metric and one pair per time-series metric to
`v2/figures`. The weathering panels report solution Ca, Mg and Si in g
element m^-2, and N2 plus total N2O as cumulative gaseous N loss. The
reported denitrification N2O share removes the saved nitrification N2O flux
before comparing denitrification N2O with N2. Filenames begin with
`Subfield_12_summary_` or
`Subfield_12_time_dynamics_`. Each figure has one axis and contains control,
basalt and calcite for the 1-D, 2-D whole-field, 2-D cropped-buffer and 2-D
native-vegetation-buffer configurations. Treatment controls colour and
configuration controls line style. Cumulative nitrate leaching and the soil
nitrate pool remain separate metrics. Set `FigureDirectory` beforehand to
write the figures somewhere else. The all-12 stage is skipped with a warning
listing the missing runs until the 2-D results finish. NEE is respiration
minus crown-weighted NPP; positive values are carbon release to the
atmosphere. Carbon is normalised by the full 3600 m^2 domain area.
Daily rows are dated by the BIOGEO day they hold: row j runs from 01:00 of
Date(1)+(j-3) to 00:00 of the next day, and row 2 is the lone first hour.
Daily WFPS averages the same hours, and the initial retained state has no
preceding interval and is plotted as `NaN`.

The completed 1-D runs can be plotted before the 2-D runs finish, from
`v2/plotting`:

```matlab
Runs = 10:12; Plot_1D_Then_All_12
```

This writes separate PNG and FIG files for each metric. With only runs 10 to
12 available, each comparison figure contains 1-D control, basalt and calcite.
When runs 1 to 9 finish, the three 2-D groups can be plotted without
reloading the 1-D results:

```matlab
Runs = 1:9; Plot_1D_Then_All_12
```

This writes separate `Subfield_2D_summary_*` and
`Subfield_2D_time_dynamics_*` PNG and FIG files. Each contains control,
basalt and calcite for the three available 2-D configurations. Running the
script without setting `Runs` plots the 1-D group first and then compares
all twelve treatment and configuration combinations.

From an open MATLAB session at the repository root:

```matlab
cd v2/run
Run_All_12_Parallel('DryRun',true)
Run_All_12_Parallel('Runs',10:12,'MaxConcurrent',3,'MaxRuntimeHours',168)
Run_All_12_Parallel('Runs',1:9,'MaxConcurrent',2,'ParallelWorkers',6,'MaxRuntimeHours',168)
cd ../plotting
Plot_1D_Then_All_12
```

The first command completes the three vertical 1-D comparisons before the
2-D runs begin. Stop between the two commands if the second stage needs to
wait. `Ctrl+C` stops active child runs; completed result files are retained
and skipped when the same command is run again. On Windows the launcher
stops the complete worker process tree, so a MATLAB worker is not left behind.
The chemistry reruns started with `run_phcharge(10)`, `run_phcharge(11)` or
`run_phcharge(12)` write a rolling checkpoint every 30 minutes to
`results/checkpoints`. A matching checkpoint resumes on the next launch. The
checkpoint is removed after the final result is saved.

Results go to `v2/SIM/Subfield_12`. Logs and worker status files go to
`v2/logs/Subfield_12`. Existing results are skipped unless
`'Overwrite',true` is passed. Use overwrite after a solver change so all
comparisons use the same code. A one-day check uses:

```matlab
cd v2/run
Run_All_12_Parallel('MaxConcurrent',2,'NDays',1,'Overwrite',true)
addpath(fullfile('..','T&C_Code','Tests')); test_subfield_12_smoke_results(1)
cd ../plotting
NDays = 1; Plot_1D_Then_All_12
```

The September 2019 mineral-N regression uses a small saved fixture:

```matlab
addpath(fullfile('v2','T&C_Code','Tests')); test_no3_immobilization_replay
```

`ParallelWorkers` only applies to the independent daily ERW chemistry in
each 2-D run. Keep `MaxConcurrent * ParallelWorkers` within the available
CPU allocation. Daily biology, hourly physics and lateral routing remain
serial. Full 12-column runs also hold several GB of state.

### Current campaigns (2026-09-25)

The calibrated set and every campaign that reads it live under `runtime/`;
`RUNBOOK_2026-09-24.md` at the checkout root orders them for Bouchet:

* `Base_1D` and `Base_2D`: the 5-year validation, all twelve runs,
  through `runtime/Base_1D/submit_all_12.sh` (optional
  `TC_CALCITE_CASE=8` adds the 5 t/ha calcite runs as a campaign of their
  own). `Campaign_Parameters.m` is the only file that spells the
  calibrated numbers, and `phobs.m` scores August pH against both
  observed depth intervals and `budget.m` prints the alkalinity budget
  and exchanger state of the control and basalt runs.
* `SSP_1D` and `SSP_2D`: the ssp126/245/585 projections to 2050 through
  `runtime/SSP_1D/submit_ssp.sh` (`TC_SSP`, `TC_HORIZON`, optional
  `TC_CALCITE_CASE=8`).
* `Calibration_Sweep_1D`: the joint sweep that calibrates alpha and
  KNmaxd, with the Kbot alternatives beside them.
* `Uncertainty_1D` and `Uncertainty_2D`: the parameter-uncertainty
  ensemble, the T&C-SMEW field paper's sensitivity analysis in the coupled
  model. Beside the original set, reused from the Base and SSP results,
  all sixteen parameters run at the lower and at the upper end of their
  ranges, in 1-D and 2-D, with the 5 t/ha arm, over both horizons
  (`runtime/Uncertainty_1D/submit_uncertainty.sh`). `Run_One_Subfield`
  carries its three knobs (`BiogeoOverrides`, `DissolutionFactor`,
  `SaturationExponent`), bit-identical when unset.
* `Worker_Scaling_Test`: sizes the 35-year 2-D jobs.

Figures: `v2/plotting/Plot_1D_Then_All_12.m` (`Scenario` for projection
results). Its `UncertaintyFile` bands and `v2/plotting/Plot_Uncertainty.m`
read the retired 36-member summary; the three-set ensemble's summary and
figures come with its post-processing package.
