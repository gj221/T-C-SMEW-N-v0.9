function RM = Run_Matrix_12()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run_Matrix_12
%%
%% Three treatments crossed with four scenario groups.
%%
%%                     | whole_field  margin_crop  margin_prairie  vertical_1d
%%   ------------------+------------------------------------------------------
%%   control           |      1            4             7             10
%%   basalt            |      2            5             8             11
%%   calcite           |      3            6             9             12
%%
%% The first three groups use the retained 12-column plot and its 5 m
%% margin. The vertical_1d group uses one 3600 m2 column. It has no
%% subdomains or lateral transport.
%%
%% All entries use the same control T&C setup and 2500 mm active
%% biogeochemical depth. Only geometry, amendment placement, feedstock
%% properties and the prairie PFT switch vary. The bounded negative-NO3
%% immobilization logic is shared by all entries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CASE_CONTROL = 5;
CASE_BASALT  = 3;
CASE_CALCITE = 7;

treatments = { ...
    'control', CASE_CONTROL; ...
    'basalt',  CASE_BASALT;  ...
    'calcite', CASE_CALCITE};

%% layout, apply_main, apply_margin, margin_pft, use_subdomains
layouts = { ...
    'whole_field',    true,  true,  '',                   true; ...
    'margin_crop',    false, true,  '',                   true; ...
    'margin_prairie', false, true,  'PFT_Prairie_Margin',true; ...
    'vertical_1d',    true,  false, '',                   false};

geom.field_side_m    = 60;
geom.margin_width_m  = 5;
geom.main_grid_nx    = 2;
geom.main_grid_ny    = 2;
%% 300 mm, matching the workbook value and the published spin-up that
%% produced Data_Kantola_BG.mat. This must not be changed in isolation:
%% ZBIOG*rsd is the gsoil/m2 conversion, so raising the depth dilutes
%% every per-area pool by the same factor when expressed per gram of
%% soil. At 2500 mm the DOC concentration entering the microbial
%% Michaelis-Menten term fell 8.3x, taking M from about 0.25 to 0.03,
%% below the mrb/(Vdb+mrb) growth threshold of 0.20. Maintenance does not
%% scale with concentration, so bacteria and fungi decayed to extinction
%% with a 100-day e-folding time while DOC accumulated unused. Running
%% deeper requires regenerating the initial conditions at that depth.
COMMON_BIOGEO_DEPTH_MM = 300;

RM = struct('id',{},'name',{},'treatment',{},'layout',{}, ...
    'case_number',{},'apply_main',{},'apply_margin',{}, ...
    'margin_pft',{},'use_subdomains',{},'biogeo_depth_mm',{},'geom',{});

k = 0;
for il = 1:size(layouts,1)
    for it = 1:size(treatments,1)
        k = k + 1;
        RM(k).id = k;
        RM(k).treatment = treatments{it,1};
        RM(k).case_number = treatments{it,2};
        RM(k).layout = layouts{il,1};
        RM(k).apply_main = layouts{il,2};
        RM(k).apply_margin = layouts{il,3};
        RM(k).margin_pft = layouts{il,4};
        RM(k).use_subdomains = layouts{il,5};
        RM(k).biogeo_depth_mm = COMMON_BIOGEO_DEPTH_MM;
        RM(k).name = sprintf('%s__%s',RM(k).treatment,RM(k).layout);
        RM(k).geom = geom;
    end
end

for k = 1:numel(RM)
    if strcmp(RM(k).treatment,'control')
        RM(k).apply_main = false;
        RM(k).apply_margin = false;
    end
end

end
