# Kantola five-year 2-D runs 1 to 9 on a field grid (Bouchet, R2025a)

Campaign tag: `Field_2D__<grid name>`

Runs 1 to 9 of `Run_Matrix_12.m` (control, basalt and calcite on the `whole_field`, `margin_crop` and `margin_prairie` layouts) on the real Kantola field instead of the 12-column plot: T&C-SMEW-N v1.0's polygon grid, with topography from the USGS 1 m DEM when the grid carries it. Everything else is the calibrated set of `runtime/Base_1D` (`Campaign_Parameters.m`, `Assert_Base_Provenance.m`), so a result here differs from the matching `Base_2D` result only by the grid and the topography.

## The grid

A grid is a MAT file holding `topo`, saved by `Build_Field_Grid_Topology` on a machine with the Mapping Toolbox:

```matlab
topo = Build_Field_Grid_Topology('kantola_field.shp', ...
    'GridCellSizeM',5,'MarginWidthM',5, ...
    'FieldDEM','.../DEM/Energy Farm_IL');      % omit FieldDEM for a flat grid
save('kantola_dem_5m.mat','topo');
```

* **Squares**: v1.0's majority rule (`Build_Projected_Field_Topology`, `'majority'`): a square is simulated when more than half of it lies in the field, on the rotation and offset that minimise the misrepresented area. The width comes from the terrain: `v2/interpolation/Run_Topographic_Grid_Convergence` steps 20, 10, 5, 3.5, 2.5, 1.5 and 1 m and takes the coarsest width from which every finer step changes the terrain statistics by less than 5 %.
* **Margin**: the field's outer band of `MarginWidthM` metres (default 5 m, the plot's) is labelled margin and the rest main (`Build_Margin_Field_Topology`), exactly as the 12-column plot keeps its margin as the outer 5 m of the 60 m plot. All nine runs simulate the same squares; the layouts place rock and the prairie by label. A square wider than the margin cannot resolve it and is refused.
* **Topography** (`FieldDEM`): elevation, slope, aspect, curvature, TWI and sky view per square (`Attach_Topography_To_Topology`), at a DEM spacing of the square width by default. In the 2-D frame (`OPT_DEM`) it drives total-head Darcy flow, overland routing of runoff and its solutes to the outlet cells, and per-square radiation.
* **Results** weigh each square by its overlap with the field (`Field_Aggregation_Areas`), so domain means are field means.

Copy the file to Bouchet **beside** the deployed root (a validated root never gains files), with a name of at most 32 letters, digits and underscores; the name becomes the campaign tag. The task refuses a file inside the root. The file's content enters the validation stamp, so replacing it after validation stops production with `StaleValidation`.

## Resources

The hourly T&C loop is serial over columns, so time and memory grow with the number of squares; the daily ERW step runs on twelve workers. The 12-column plot takes about 6 h and under 64 GB for five years. Set `TC_FIELD_MEM` and `TC_FIELD_TIME` per grid from its cell count (`topo.N_cells`). The week partition lets one user hold 1.5 TB at once, which sets how many elements run together.

## Submission

From the checkout root, after the equivalence baseline has been regenerated for this deployment:

```bash
TC_FIELD_TOPOLOGY=/full/path/beside/the/root/kantola_dem_5m.mat \
TC_FIELD_MEM=320G TC_FIELD_TIME=4-00:00:00 \
MATLAB_MODULE=MATLAB/2025a TC_BOUCHET_ACCOUNT=<pi account> \
bash runtime/Field_2D/submit_field_2d.sh
```

The controller submits one validation element (`day`, twelve CPUs, 96 GB, 4 h: the complete audit, the Linux portability checks and a one-day run of the grid across the first application) and the production array (`week`, runs 1 to 9), which starts only after the validation passes. `TC_FIELD_STAGE=resume` continues after a timeout or node loss: finished runs are kept, interrupted ones resume from their two-hourly checkpoints, and the controller lists any stale run lock to remove first.

## Files

| file | role |
|---|---|
| `submit_field_2d.sh` | controller: preflight, validation element, production array |
| `submit_field_2d.sbatch` | array batch file: validation mode and production |
| `Run_Field_2D_Task.m` | one element: gates, lock, existing-result check, `Run_Subfield_Worker` with `FieldTopologyFile`, post-run checks, provenance (calibrated set and the grid file's SHA-256), lateral conservation |
| `Check_Lateral_Conservation.m` | byte copy of the 2-D campaigns' whole-run check, which also reports the overland export of a topographic run |

## Outputs

* Results: `v2/SIM/Subfield_12/<tag>/TC_2D_runNN_<treatment>__<layout>__field_<shapefile>__cell<w>m__margin<m>[__dem<r>]__alpha3p4__...__1826d__from20160101.mat`.
* Logs, status, locks and stamp: `v2/logs/Subfield_12/<tag>/`.
* Slurm output: `runtime/Field_2D/slurm_<array job>_<element>.out|err`.

Compare with `Base_2D`: the same run on the 12-column plot. A flat grid of the same field (no `FieldDEM`) separates the effect of the grid from that of the topography.
