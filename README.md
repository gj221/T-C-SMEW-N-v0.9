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
leaching remain active.

From PowerShell at the repository root:

```powershell
matlab -batch "cd('v2'); Run_All_12_Parallel('MaxConcurrent',2)"
matlab -batch "cd('v2'); Plot_Subfield_12"
```

From an open MATLAB session in `v2`:

```matlab
Run_All_12_Parallel('DryRun',true)
Run_All_12_Parallel('MaxConcurrent',2)
Plot_Subfield_12
```

Results go to `v2/SIM/Subfield_12`. Logs and worker status files go to
`v2/logs/Subfield_12`. Existing results are skipped unless
`'Overwrite',true` is passed. A one-day check uses:

```matlab
Run_All_12_Parallel('MaxConcurrent',2,'NDays',1,'Overwrite',true)
Plot_Subfield_12('NDays',1)
```

Two workers is the default for a 32 GB machine. A full 12-column worker
holds several GB of state, so use three only when more memory is available.
