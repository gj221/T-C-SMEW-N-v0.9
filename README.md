# T&C-SMEW-N v0.9

## Repository layout

```
v2/run/          launchers: ERW_run, Run_Matrix_12, Run_One_Subfield,
                 Run_All_12_Parallel, Run_Subfield_Worker,
                 Terminate_Process_Tree
v2/plotting/     figure scripts: Plot_1D_Then_All_12, Plot_Grid_Setup,
                 Plot_Biogeo_Pools, Plot_Nflux_Diagnostics,
                 Plot_Kantola_Benchmark, and the Cflux_post_processing helper
v2/figures/      every PNG and FIG the plotting scripts write
v2/T&C_Code/     model source, Inputs, Parameters and Tests
v2/SIM/          simulation results (not tracked)
v2/logs/         worker logs and status files (not tracked)
v2/Experimental observations/   digitised field data used by the benchmark
```

## Site scope

This v0.9 checkout supports the Energy Farm (Kantola) simulations and Bruce
Field forcing preparation. The Bouchet twelve-run campaign below targets the
Energy Farm only. Bruce Field support is separate from this campaign.

### Bruce Field ERA5 forcing (full checkout only)

This workflow remains in the full v0.9 checkout. The Energy Farm Bouchet
archive omits its generated forcing files and site-specific wrapper.

`Preprocess_Bruce_Field_ERA5` reads the split ERA5-Land time-series NetCDF
files and writes `Data_Bruce_Field_Spinup_2014_2023.mat` and
`Data_Bruce_Field_Transient_2024_2025.mat` in `v2/T&C_Code/Inputs`.
The period labels use UTC boundaries. `Date` is fixed Eastern Standard Time
(GMT-5), as required by the radiation partitioning calculation.

From the repository root:

```matlab
addpath(fullfile('v2','T&C_Code','Inputs','Preprocess'))
Preprocess_Bruce_Field_ERA5
```

The ERA5 files contain the grid point at 36.7 N, 78.5 W. The preprocessing
uses a 127.4689 m USGS 3DEP elevation at that point. Negative precipitation
and shortwave values within the observed NetCDF packing noise are set to zero.

Every plotting entry point is a runnable script: `cd` into `v2/plotting` and
run it by name, with no arguments. Each one puts the folders it needs on the
path, so it works from any working directory. `Cflux_post_processing` stays a
function because the benchmark script calls it as a helper.

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
rates, including nitrification and denitrification, use the reconstructed 1:2
DI-water bulk-extract pH. The Kantola seed uses the measured pretreatment bulk
pH as an inversion boundary while retaining the prescribed exchange inventory.
There is no additive pH offset.

Kantola ERW runs require `appelo` because their conserved initial exchange
inventory is built with the same Appelo constants used by the daily solver. A
missing scheme or `TC_GT` is rejected before Kantola initialization. `TC_GT`
remains available for legacy-compatible non-Kantola ERW configurations, where
it keeps the legacy T&C fixed-pool rates and disables only the kinetic
solution-exchange rates. Standalone T&C runs keep all legacy rates.

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
plateau at pH 7.8. The shape is adapted from Blanc-Betes et al. (2021), which
used a pH 7.0 breakpoint. The 7.8 substitution is a treatment sensitivity
supported by Simek et al. (2002), not a universal optimum. The N2:N2O pH
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

Do not reuse the first Bouchet runs 10 to 12 for treatment comparisons. The
amendment schedule used the allocation-vector index rather than the timestamped
state index. This placed each event one 10-minute step early and dropped an
event at the exact run start.

The corrected 1-D runbook is in
`runtime/Kantola_5yr_1D_Runs10to12_ScheduleFix_R2025a_Bouchet`. It runs the
three treatments in parallel on one physical node, with one CPU and one
numerical thread per treatment. The 2-D runbook is in
`runtime/Kantola_5yr_2D_Runs01to09_DenitPH7p8_Alpha1p3_Bouchet`. It uses six
thread workers for the independent daily ERW column calculations.

Both runbooks pin `cpumodel:9575f`, record the assigned hardware and require a
fresh R2025a validation stamp. Use the coordinated controller below for the
twelve-run campaign. Calling either batch file directly bypasses its account,
environment, clean-state and dependency checks. Production stops if source,
inputs, tests or the numerical platform changed after validation.

Build the release from the reviewed working-tree manifest. The manifest
includes `v2/POST_LEGACY_EQUATION_AUDIT.md` with the source and unit rationale
for this correction. It excludes old simulation and runtime output by
construction:

```bash
set -euo pipefail
if ! tar --version | grep -Fq 'GNU tar'; then
  echo 'GNU tar is required to build this archive.' >&2
  exit 1
fi
tar --format=ustar -czf tc-smewn-v09-energy-farm-p-return-bouchet.tar.gz \
  --verbatim-files-from --numeric-owner --owner=0 --group=0 \
  --mode='u=rwX,go=rX' --no-acls --no-xattrs --no-selinux \
  --exclude='*.asv' --exclude='*.log' --exclude='*.status' \
  --exclude='slurm_*.out' --exclude='slurm_*.err' --exclude='.DS_Store' \
  --exclude='v2/T&C_Code/Inputs/Data_Bruce_Field_*.mat' \
  --exclude='v2/T&C_Code/Inputs/Preprocess/Preprocess_Bruce_Field_ERA5.m' \
  -T runtime/bouchet_all_12_package_manifest.txt
sha256sum tc-smewn-v09-energy-farm-p-return-bouchet.tar.gz > tc-smewn-v09-energy-farm-p-return-bouchet.tar.gz.sha256.txt
```

The Energy Farm archive omits the generated Bruce forcing files and the
Bruce-specific wrapper. They remain in the full v0.9 checkout.

Create a new, empty Bouchet repository root with the release archive and its
checksum sidecar in the parent directory. Do not extract it over the earlier
`tc-smewn-v09` tree because that tree contains results with an older campaign
fingerprint. Confirm the account reported by `slurm_checkup.sh`, then run this
block from the new root:

```bash
set -euo pipefail
(cd .. && sha256sum -c tc-smewn-v09-energy-farm-p-return-bouchet.tar.gz.sha256.txt)
tar --no-same-owner -xzf ../tc-smewn-v09-energy-farm-p-return-bouchet.tar.gz
test -r 'v2/T&C_Code/Inputs/Data_Kantola.mat'
test -r 'v2/T&C_Code/Inputs/Data_Kantola_BG.mat'
test -r 'v2/T&C_Code/Tests/run_post_legacy_equation_audit_tests.m'
test -r 'v2/POST_LEGACY_EQUATION_AUDIT.md'
MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=pi_np363 bash runtime/submit_bouchet_all_12.sh
```

The controller submits all five jobs held, builds the complete dependency
graph, then releases them. The R2025a baseline runs first. After it succeeds,
the 1-D and 2-D validation allocations become eligible independently and may
run at the same time. The 1-D validation uses four MATLAB process workers
inside one 6-CPU, 40 GB allocation for four 315-day identity trajectories:
control, basalt, calcite and repeated calcite. The 2-D validation uses array
element 1 in a separate 6-CPU, 32 GB allocation.

Both production allocations wait for both validations. After both pass, the
nine-element 2-D array and the three-task 1-D allocation become eligible
together. The array permits up to nine active elements, and the 1-D allocation
starts three exclusive one-CPU MATLAB steps after Slurm assigns its node. All
twelve simulations can therefore run concurrently when Bouchet has the
resources, but eligibility does not guarantee simultaneous start times. No
Priority Tier QOS is requested.

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
Daily WFPS is calculated from exact 24-hour blocks. The initial retained state
has no preceding 24-hour interval and is plotted as `NaN`; each later state is
paired with the mean for the day that ends at that state.

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
Active runs restart from the beginning because transient checkpoints are not
written.

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
