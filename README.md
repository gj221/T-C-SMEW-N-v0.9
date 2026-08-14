# T-C-SMEW-N-v0.9

The nine sub-field simulations are defined in `v2/Run_Matrix_9.m`.
The default domain has 12 columns: four 25 m by 25 m field cells, four
50 m by 5 m side margins, and four 5 m by 5 m corners. This keeps the
earlier geometry and avoids the cost of a uniform 5 m grid.

From PowerShell at the repository root:

```powershell
matlab -batch "cd('v2'); Run_All_9_Parallel('MaxConcurrent',3)"
```

From an open MATLAB session in the `v2` directory:

```matlab
Run_All_9_Parallel('DryRun',true)
Run_All_9_Parallel('MaxConcurrent',3)
```

Results go to `v2/SIM/Subfield_9`. Child-process logs and status files
go to `v2/logs/Subfield_9`. Existing results are skipped unless
`'Overwrite',true` is passed. Shortened runs and grid overrides use
separate filenames, so they cannot cause a full default run to be skipped.
