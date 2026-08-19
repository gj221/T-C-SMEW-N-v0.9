# T&C-SMEW-N v0.9

The twelve simulations are defined in `v2/Run_Matrix_12.m`. Each group
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
state. The active biogeochemical depth is 2500 mm in every run, set by
the deepest configured prairie root depth. Treatment rows only supply
the amendment properties. Untreated cells use the control chemistry row.
The configuration workbook contains only the five Kantola cases used by
the single-column comparison and the twelve-run matrix.

The archived Kantola background file contains areal pool totals and no
depth metadata. It is reused unchanged in every case. This keeps the runs
comparable, but it is not a new 2500 mm equilibrium spin-up.

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

From PowerShell at the repository root:

```powershell
matlab -batch "cd('v2'); Run_All_12_Parallel('Runs',10:12,'MaxConcurrent',3,'MaxRuntimeHours',168)"
matlab -batch "cd('v2'); Run_All_12_Parallel('Runs',1:9,'MaxConcurrent',2,'MaxRuntimeHours',168)"
matlab -batch "cd('v2'); Plot_1D_Then_All_12"
```

The plotting command first draws the vertical 1-D weathering and nitrogen
figure (`vertical_1d_weathering_nitrogen.png` and `.fig`) and then writes
`Subfield_12_summary.png`, `Subfield_12_summary.fig`,
`Subfield_12_time_dynamics.png` and `Subfield_12_time_dynamics.fig` to
`v2/SIM/Subfield_12`. The all-12 stage is skipped with a warning listing
the missing runs until the 2-D results finish. The time plot
uses two columns per scenario group and no more than three stacked rows.
It shows control, basalt and calcite through time for cumulative N2O,
cumulative nitrate leaching, the soil nitrate pool, pore-water pH, NEE
and cumulative crop yield. NEE is respiration minus crown-weighted NPP;
positive values are carbon release to the atmosphere. Crop yield is
harvested fruit/grain carbon in g C m^-2, so prairie PFTs contribute zero
yield. Both carbon metrics use the full 3600 m^2 domain area as the
denominator.

The completed 1-D runs can be plotted before the 2-D runs finish:

```matlab
Runs = 10:12; Plot_1D_Then_All_12
```

This writes updated six-panel PNG and FIG outputs for runs 10 to 12 only.
When runs 1 to 9 finish, the three 2-D groups can be plotted without
reloading the 1-D results:

```matlab
Runs = 1:9; Plot_1D_Then_All_12
```

This writes PNG and FIG copies of `Subfield_2D_summary` and
`Subfield_2D_time_dynamics`. Running the script without setting `Runs`
plots the 1-D group first and then all 12 scenarios together.

From an open MATLAB session in `v2`:

```matlab
Run_All_12_Parallel('DryRun',true)
Run_All_12_Parallel('Runs',10:12,'MaxConcurrent',3,'MaxRuntimeHours',168)
Run_All_12_Parallel('Runs',1:9,'MaxConcurrent',2,'MaxRuntimeHours',168)
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
Run_All_12_Parallel('MaxConcurrent',2,'NDays',1,'Overwrite',true)
addpath(fullfile('T&C_Code','Tests')); test_subfield_12_smoke_results(1)
NDays = 1; Plot_1D_Then_All_12
```

The September 2019 mineral-N regression uses a small saved fixture:

```matlab
addpath(fullfile('T&C_Code','Tests')); test_no3_immobilization_replay
```

Two 2-D workers is the default for a 32 GB machine. The three 1-D runs can
run together. Full 12-column workers hold several GB of state, so increase
the second-stage limit only after checking available memory.
