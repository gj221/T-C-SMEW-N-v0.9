function boundary_polygon = Field_Setup_Polygon(field_setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rebuild the exact field boundary of a setup in its rotated grid frame.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(field_setup) || ~isscalar(field_setup)
    error('Field_Setup_Polygon:BadFieldSetup', ...
        'field_setup must be one grid setup or site boundary.');
end

% A site boundary already holds the union as a polyshape. Rebuilding it from
% NaN-separated vertices would merge its disjoint fields into one ring.
if isfield(field_setup,'union_polygon')
    boundary_polygon = field_setup.union_polygon;
    if isempty(boundary_polygon.Vertices) || area(boundary_polygon) <= 0
        error('Field_Setup_Polygon:EmptyBoundary', ...
            'The setup boundary must enclose positive area.');
    end
    return
end
if ~isfield(field_setup,'field_boundary_x_m') || ...
        ~isfield(field_setup,'field_boundary_y_m')
    error('Field_Setup_Polygon:BadFieldSetup', ...
        'field_setup must carry its field boundary coordinates.');
end

boundary_polygon = polyshape( ...
    field_setup.field_boundary_x_m,field_setup.field_boundary_y_m, ...
    'Simplify',false,'KeepCollinearPoints',true);
if isempty(boundary_polygon.Vertices) || area(boundary_polygon) <= 0
    error('Field_Setup_Polygon:EmptyBoundary', ...
        'The setup boundary must enclose positive area.');
end

end
