function RM = Run_Matrix_9()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run_Matrix_9
%%
%% The nine sub-field simulations: three treatments crossed with three
%% spatial layouts, all on the sub-plot domain with 5 m field margins.
%%
%%                     | whole_field    margin_crop     margin_prairie
%%   ------------------+----------------------------------------------
%%   control           |    1               4                7
%%   basalt            |    2               5                8
%%   calcite           |    3               6                9
%%
%% TREATMENT  what is applied
%%   control   no amendment
%%   basalt    crushed basalt      (workbook case 3, Blueridge)
%%   calcite   agricultural lime   (workbook case 7)
%%
%% LAYOUT  where the amendment goes, and what grows in the margin
%%   whole_field     amendment over main grid AND margin; margin cropped
%%   margin_crop     amendment on the 5 m margin only; margin cropped
%%   margin_prairie  amendment on the 5 m margin only; margin converted
%%                   to prairie (big bluestem + purple coneflower)
%%
%% All nine carry the fixed negative-NO3 pool: under net immobilization
%% the N demand is charged to NH4 and NO3 in proportion to pool size
%% rather than at a fixed 90:10, so the NO3 sink is bounded by the pool
%% (BIOGEOCHEMISTRY_DYNAMIC3.m, inherited from main).
%%
%% NOTE ON A DEGENERATE PAIR
%%   Runs 1 and 4 (control x whole_field, control x margin_crop) are the
%%   same physical configuration: no amendment anywhere and crop in the
%%   margin. They are kept as separate entries so each layout subset has
%%   its own baseline. The launcher runs both because all nine matrix
%%   entries are retained as separate outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Workbook case numbers in ERW_MOD_PARAM_parameters.xlsx (Kantola rows)
CASE_CONTROL = 5;    %% crop_control  -- no rock
CASE_BASALT  = 3;    %% crop_Beerling -- Blueridge basalt
CASE_CALCITE = 7;    %% crop_calcite  -- agricultural lime

treatments = { ...
    'control', CASE_CONTROL; ...
    'basalt',  CASE_BASALT;  ...
    'calcite', CASE_CALCITE};

%% apply_main / apply_margin: which regions receive the amendment
%% margin_pft: '' keeps the base crop parameters in the margin cells
layouts = { ...
    'whole_field',    true,  true,  ''; ...
    'margin_crop',    false, true,  ''; ...
    'margin_prairie', false, true,  'PFT_Prairie_Margin'};

%% --- Domain geometry, shared by all nine ------------------------------
geom.field_side_m    = 60;    %% [m] square plot side, main grid + margin
geom.margin_width_m  = 5;     %% [m] field margin width
geom.main_grid_nx    = 2;     %% main grid columns
geom.main_grid_ny    = 2;     %% main grid rows
%% 2x2 main grid + 4 side strips + 4 corners = 12 columns. Raise
%% main_grid_nx/ny to 4 for a 24-column run when the extra interior
%% resolution is wanted; cost scales with the column count.

RM = struct('id', {}, 'name', {}, 'treatment', {}, 'layout', {}, ...
    'case_number', {}, 'apply_main', {}, 'apply_margin', {}, ...
    'margin_pft', {}, 'geom', {});

k = 0;
for il = 1:size(layouts,1)
    for it = 1:size(treatments,1)
        k = k + 1;
        RM(k).id           = k;
        RM(k).treatment    = treatments{it,1};
        RM(k).case_number  = treatments{it,2};
        RM(k).layout       = layouts{il,1};
        RM(k).apply_main   = layouts{il,2};
        RM(k).apply_margin = layouts{il,3};
        RM(k).margin_pft   = layouts{il,4};
        RM(k).name         = sprintf('%s__%s', RM(k).treatment, RM(k).layout);
        RM(k).geom         = geom;
    end
end

%% A control run has no rock, so "where the rock goes" is meaningless;
%% only the margin vegetation distinguishes control layouts.
for k = 1:numel(RM)
    if strcmp(RM(k).treatment, 'control')
        RM(k).apply_main   = false;
        RM(k).apply_margin = false;
    end
end

end
