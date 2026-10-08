# Kantola five-year 2-D runs 1 to 9 on a field grid (Bouchet, R2025a)

Campaign tag: `Field_2D__<grid name>`

Runs 1 to 9 of `Run_Matrix_12.m` (control, basalt and calcite on the `whole_field`, `margin_crop` and `margin_prairie` layouts) on the real Kantola field instead of the 12-column plot: T&C-SMEW-N v1.0's polygon grid, with topography from the USGS 1 m DEM when the grid carries it. Everything else is the calibrated set of `runtime/Base_1D` (`Campaign_Parameters.m`, `Assert_Base_Provenance.m`), so a result here differs from the matching `Base_2D` result only by the grid and the topography.

## The grid

A grid is a MAT file holding `topo`, saved by `Build_Field_Grid_Topology` on a machine with the Mapping Toolbox:

```matlab
topo = Build_Field_Grid_Topology('kantola_basalt_10m.shp', ...
    'GridCellSizeM',20,'MarginWidthM',20, ...
    'FieldDEM','.../DEM/Energy Farm_IL');      % omit FieldDEM for a flat grid
save('kantola_basalt_dem_20m.mat','topo');
```

* **Squares**: v1.0's majority rule (`Build_Projected_Field_Topology`, `'majority'`): a square is simulated when more than half of it lies in the field, on the rotation and offset that minimise the misrepresented area. `v2/interpolation/Run_Topographic_Grid_Convergence` steps 20, 10, 5, 3.5, 2.5, 1.5 and 1 m and reports, per terrain statistic, the coarsest width from which every finer step changes it by less than 5 %.
* **Margin**: the field's outer band of `MarginWidthM` metres (default 5 m, the plot's) is labelled margin and the rest main (`Build_Margin_Field_Topology`), as the 12-column plot keeps its margin as the outer 5 m of the 60 m plot. All nine runs simulate the same squares; the layouts place rock and the prairie by label. A square wider than the margin cannot resolve it and is refused.
* **Topography** (`FieldDEM`): elevation, slope, aspect, curvature, TWI and sky view per square (`Attach_Topography_To_Topology`), at a DEM spacing of the square width by default. In the 2-D frame (`OPT_DEM`) it drives total-head Darcy flow, overland routing of runoff and its solutes to the outlet cells, and per-square radiation.
* **Results** weigh each square by its overlap with the field (`Field_Aggregation_Areas`), so domain means are field means.

Copy the file to Bouchet **beside** the deployed root (a validated root never gains files), with a name of at most 32 letters, digits and underscores; the name becomes the campaign tag. The task refuses a file inside the root. The file's content enters the validation stamp, so replacing it after validation stops production with `StaleValidation`.

## The Kantola grids (2026-10-08)

The two 3.8 ha maize/soybean trial fields of Kantola et al. (2023), each with its eddy-covariance tower (`MaizeBasalt`, `MaizeNoBasalt`):

* **Boundaries**: their 10 m boundaries are kept in `Model case study data/Maps/Energy Farm_IL`, with a README and a preview. They are the field edges measured on NAIP (25 Jul 2017, 1 m), reduced to the Sentinel-2 10 m grid by the majority rule. Sentinel-2 confirms them: 2016-2020 crops of soy, maize, maize, soy, maize inside both. The grids are in NAD83 / UTM 16N, the DEM's frame.
* **Squares and margin**: 20 m squares and a 20 m margin, one lime-spreader pass (about 60 ft), so the margin is the outer ring of squares.
* **Width**: the convergence ladder gives no joint width.
  * On the basalt field, slope and relief stay within 5 % from 20 m down.
  * TWI keeps falling as the 1 m lidar resolves microtopography, on both fields, so it never converges.

| grid file | squares (main + margin) | simulated | field left out | relief | mean slope | outlets |
|---|---|---|---|---|---|---|
| `kantola_basalt_dem_20m` | 81 (49 + 32) | 32,400 m2 | 1,600 m2 | 4.15 m | 2.2 % | 1 |
| `kantola_control_dem_20m` | 90 (56 + 34) | 35,000 m2 | 1,900 m2 | 1.86 m | 0.95 % | 3 |

* **Flat twins**: `kantola_basalt_flat_20m` and `kantola_control_flat_20m` hold the same squares and labels without topography.
* **Field left out**: the strip along one edge that is exactly half a square wide, which the majority rule drops.

## Resources

The hourly T&C loop is serial over columns, so time and memory grow with the number of squares; the daily ERW step runs on twelve workers.

* **Measured, five years:**
  * one column: 2 h 14 min and 3.1 GB;
  * the 12-column plot: about 6 h, peaks of 18-27 GiB.
* **Per column:** the daily ERW chemistry takes about 2.1 h, and the serial hourly step about 0.16 h.
* **An N-square grid:** therefore takes about N x 0.33 h with twelve workers, and about 1.4-2.2 GiB per square.
* **The Kantola grids:** that is 27-30 h and 110-200 GiB. Submit them with `TC_FIELD_MEM=256G TC_FIELD_TIME=3-00:00:00`.
* **Concurrency:** the week partition lets one user hold 1.5 TB at once, so five such elements run together.
* **Overruns:** an element that times out or runs out of memory keeps its two-hourly checkpoint. Resume it with `TC_FIELD_STAGE=resume`, raising the limit first.

## Submission

From the checkout root, after the equivalence baseline has been regenerated for this deployment:

```bash
TC_FIELD_TOPOLOGY=/full/path/beside/the/root/kantola_basalt_dem_20m.mat \
TC_FIELD_MEM=256G TC_FIELD_TIME=3-00:00:00 \
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

## Figures

`v2/plotting` maps any result cell by cell (no display needed; PNG, FIG, CSV and MAT go to `v2/figures/field_maps/<result>`):

```bash
matlab -batch "addpath('v2/plotting'); f = '<result>.mat'; \
  Plot_Field_Layout(f,'Visible','off'); \
  Plot_Field_Maps(f,'Variables',{'N2O','NO3_leaching','NEE','pH_bulk','lateral_water'},'Interval','year','Visible','off'); \
  Plot_Field_Maps(f,'Variables',{'N2O','pH_bulk'},'Interval','season','Reference','<control result>.mat','Visible','off')"
```

`Plot_Field_Layout` draws the squares, the margin band, the vegetation and feedstock of each square and, on a DEM grid, elevation with contours, the overland flow paths, the outlets, slope and wetness index. `Plot_Field_Maps` draws the chosen variables (`Field_Cell_Variables('list')`) per year, growing season, season, month, week, day, N-day block, date range or snapshot date, with the field-mean series above the maps, optionally as the difference from another run on the same grid or as a GIF.
