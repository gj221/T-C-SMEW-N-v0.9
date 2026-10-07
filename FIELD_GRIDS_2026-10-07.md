# Package x: v1.0's DEM topography, lateral transfer and field grids (2026-10-07)

Package x is package w plus T&C-SMEW-N v1.0's topography and field-grid machinery, ported for the Kantola runs. Without a field option every run is package w's, bit for bit.

## What came over from v1.0 (branch feature/topographic-lateral-flow)

* **Topography** (`OPT_DEM`, set by Run_One_Subfield when the grid carries a DEM): `Load_Field_DEM` (USGS GeoTIFF tiles, area-mean resampled to the square width), `Attach_Topography_To_Topology` and `Field_Topographic_Attributes` (elevation, slope, aspect, curvature, TWI, sky view, Quinn multiple-flow-direction weights, outlet cells). In the 2-D frame: total-head Darcy flow (H = psi + z) on slope-normal faces, slope-aware lateral heat, overland routing of runoff and its dissolved load (T&C pools and the ERW Na, Al, An and DIC) to the outlet cells with field-origin tracers, run-on offered to the next hour, and per-square radiation (beam projected on slope and aspect, diffuse on sky view, `Topographic_Radiation_Projection`).
* **Field grids**: v1.0's polygon grid on the field shapefile under the majority rule (a square is simulated when more than half of it lies in the field, on the rotation and offset of least misrepresented area), the width from the terrain (`v2/interpolation/Run_Topographic_Grid_Convergence`), domain results weighted by each square's overlap with the field (`Field_Aggregation_Areas`).
* **Per-cell soils** (`OPT_CELL_INPUTS`, `CellAgronomicFile`): v1.0's `Resolve_Cell_Agronomic_Inputs` and `Configure_Cell_Agronomic_Inputs`, without v1.0's per-cell spin-up states. Not used for Kantola (no soil samples).

Both switches are explicit and default off; the frame refuses a topology carrying elevation with `OPT_DEM` off, and `OPT_DEM` on without one. The port follows the PFAS fork's switch design (commit 48fb515) with the PFAS module left out. Package w's own lateral solver is kept (anisotropy `aR`, ice impedance), with only the elevation term and the face cosine added.

## New in this tree: the margin on a field grid

v1.0 grids whole fields only, every square main. The 12-column plot keeps its 5 m margin as the outer ring of the 60 m plot (2500 m2 main inside 1100 m2 of margin), so a field grid keeps it the same way: `Build_Margin_Field_Topology` labels a square main when more than half of it lies in the field shrunk inward by `MarginWidthM` (default 5 m) and margin otherwise. The squares are v1.0's for every layout; all nine runs simulate the same grid and place rock and the prairie by label. A square wider than the margin is refused.

## Running on Bouchet

Building a grid needs the Mapping Toolbox, so grids are built locally and saved:

```matlab
topo = Build_Field_Grid_Topology(shapefile,'GridCellSizeM',w,'MarginWidthM',5,'FieldDEM',dem);
save('<name>.mat','topo');
```

The file goes beside the deployed root (never inside it). `runtime/Field_2D` runs the five-year runs 1 to 9 with the calibrated Base_1D set on that grid; see its README. Deploy as a new root, regenerate the equivalence baseline (`regen_equiv_1d.sbatch`), then submit `runtime/Field_2D/submit_field_2d.sh` with `TC_FIELD_TOPOLOGY`, `TC_FIELD_MEM` and `TC_FIELD_TIME`. Its validation element runs the complete audit before any production element starts.

## Files

New: `Build_Field_Grid_Topology.m`, `Build_Margin_Field_Topology.m`, 41 functions from v1.0 (topography, routing, polygon grids, per-cell inputs; `v2/interpolation` for the width analysis), `runtime/Field_2D`, tests `test_total_head_darcy`, `test_surface_runoff_routing`, `test_field_switches`, `test_field_grid_majority`, `test_margin_field_topology`, `test_field_aggregation`, `test_cell_agronomic_inputs`, `test_2d_topographic_smoke` (all in the audit list; the margin test needs the Mapping Toolbox and reports a skip without it). Changed: `MAIN_FRAME_BG_ERW_2D.m`, `Run_One_Subfield.m`, `Darcy_Lateral_Topology.m`, `Soil_Lateral_Heat_Transport.m`, `Build_Shapefile_Grid_Setups.m` (margin argument), `Check_Lateral_Conservation.m` (reports overland export).
