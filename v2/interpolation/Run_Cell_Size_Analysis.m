function [analysis,recommendation] = Run_Cell_Size_Analysis(shape_file, ...
        feature_indices,cell_sizes_m,minimum_cell_fractions,summary_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Score a range of cell widths against the per-cell rule, for each polygon.
%%
%% The retained set keeps only cells with at least the requested fraction of
%% their own area inside the polygon, so no simulated column is mostly
%% outside the field. Because the failing cells are dropped, the cost is a
%% thin unsimulated band along the boundary; the largest width that still
%% captures a given share of the polygon is reported for each threshold.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

script_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(script_directory);
addpath(script_directory,fullfile(v2_directory,'T&C_Code'));

if nargin < 3 || isempty(cell_sizes_m)
    cell_sizes_m = [5 4 3 2.5 2 1.5 1.25 1 0.8 0.6 0.5 0.4 0.3 0.25 0.2];
end
if nargin < 4 || isempty(minimum_cell_fractions)
    minimum_cell_fractions = [0.90 0.95 0.99];
end
if nargin < 5 || isempty(summary_name)
    [~,shape_name] = fileparts(shape_file);
    summary_name = sprintf('Cell_Size_Analysis_%s', ...
        regexprep(shape_name,'[^A-Za-z0-9]+','_'));
end
capture_targets = [0.98 0.99 0.995];

setups = Build_Shapefile_Grid_Setups(shape_file,[],feature_indices);
features = shaperead(shape_file);
analysis = table();
recommendation = table();
for k = 1:numel(setups)
    label = Field_Setup_Label(features,setups(k).source_feature_indices,k);
    % Orient_Field_Polygon_Grid picks its rotation to fit the requested cell
    % width, so a candidate width other than the calibrated one needs a
    % fresh orientation search from the true, unrotated field boundary --
    % exactly what Build_Projected_Field_Topology does for a prescribed
    % width. Reusing the calibrated rotation (Field_Setup_Polygon) would
    % compose a second, unrelated re-orientation on top of it and silently
    % undercount cells at any width the calibration did not choose.
    field_polygon = polyshape( ...
        setups(k).field_boundary_projected_x_m, ...
        setups(k).field_boundary_projected_y_m, ...
        'Simplify',false,'KeepCollinearPoints',true);
    candidates = unique([setups(k).cell_size_m;cell_sizes_m(:)],'stable');
    fprintf(['\n=== %s: %.1f m2, %.1f m perimeter, ' ...
        'calibrated %.3f m with %d cells ===\n'], ...
        label,setups(k).field_area_m2,perimeter(field_polygon), ...
        setups(k).cell_size_m,setups(k).N_cells);
    polygon_analysis = Analyze_Field_Cell_Size( ...
        field_polygon,candidates,minimum_cell_fractions);
    polygon_analysis.Polygon = repmat(k,height(polygon_analysis),1);
    polygon_analysis.Field = ...
        repmat(string(label),height(polygon_analysis),1);
    polygon_analysis.Calibrated = ...
        abs(polygon_analysis.CellSizeM-setups(k).cell_size_m) < 1e-9;
    analysis = [analysis;polygon_analysis]; %#ok<AGROW>

    for threshold = reshape(minimum_cell_fractions,1,[])
        rows = polygon_analysis(abs( ...
            polygon_analysis.MinimumCellFraction-threshold) < 1e-12,:);
        rows = sortrows(rows,'CellSizeM','descend');
        fprintf('  cells kept at >=%.0f%% of their own area inside:\n', ...
            100*threshold);
        for row = 1:height(rows)
            fprintf(['    h=%6.3f m  kept %6d of %6d  captured %7.3f%%  ' ...
                'missed %7.1f m2 (%5.3f m band)  worst kept cell %6.2f%%\n'], ...
                rows.CellSizeM(row),rows.CellsStrict(row), ...
                rows.CellsAny(row),100*rows.CapturedFraction(row), ...
                rows.MissedAreaM2(row),rows.MissedBandM(row), ...
                100*(1-rows.StrictOutsideFraction(row)));
        end
        for target = capture_targets
            feasible = rows(rows.CapturedFraction >= target,:);
            if isempty(feasible)
                fprintf(['    capture >=%.1f%% not reached ' ...
                    'by any tested width\n'],100*target);
                continue
            end
            best = feasible(1,:);
            fprintf(['    capture >=%.1f%%: h=%.3f m, %d cells, ' ...
                'captured %.3f%%\n'],100*target,best.CellSizeM, ...
                best.CellsStrict,100*best.CapturedFraction);
            recommendation = [recommendation; ...
                table(k,string(label),threshold,target,best.CellSizeM, ...
                best.CellsStrict,best.CellsAny,best.CapturedFraction, ...
                best.MissedAreaM2,best.StrictOutsideFraction, ...
                'VariableNames',{'Polygon','Field','MinimumCellFraction', ...
                'CaptureTarget','CellSizeM','CellCount','CellsAny', ...
                'CapturedFraction','MissedAreaM2', ...
                'StrictOutsideFraction'})]; %#ok<AGROW>
        end
    end
end

analysis = movevars(analysis,{'Polygon','Field'},'Before','CellSizeM');
figure_directory = fullfile(v2_directory,'figures');
if ~isfolder(figure_directory); mkdir(figure_directory); end
writetable(analysis,fullfile(figure_directory,[char(summary_name) '.csv']));
if ~isempty(recommendation)
    writetable(recommendation,fullfile(figure_directory, ...
        [char(summary_name) '_recommended.csv']));
end

end
