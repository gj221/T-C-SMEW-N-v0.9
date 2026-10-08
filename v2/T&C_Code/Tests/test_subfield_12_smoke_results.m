function test_subfield_12_smoke_results(n_days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saved-result checks for the twelve-scenario smoke run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(n_days); n_days = 1; end

fprintf('\n=== test_subfield_12_smoke_results (%d day) ===\n',n_days);
n_pass = 0;
v2 = fileparts(fileparts(fileparts(mfilename('fullpath'))));
%% Run_Matrix_12 lives in v2/run, not in T&C_Code; without this the test
%% only passes when the caller happens to have v2/run on the path already.
addpath(fullfile(v2,'run'),fullfile(v2,'T&C_Code'));
result_dir = fullfile(v2,'SIM','Subfield_12');
RM = Run_Matrix_12();

metadata_ok = true;
geometry_ok = true;
placement_ok = true;
state_ok = true;
lateral_ok = true;
water_boundary_ok = true;
min_NO3 = inf;
max_heat_residual_relative = 0;
max_erw_residual = 0;

for k = 1:numel(RM)
    file = fullfile(result_dir,sprintf('TC_2D_run%02d_%s__%dd.mat', ...
        k,RM(k).name,n_days));
    if ~isfile(file)
        error('test_subfield_12_smoke_results:MissingResult', ...
            'Missing smoke result: %s',file);
    end
    wanted = {'R','N_col','MOD_PARAM_input','topo','ColCfg', ...
        'P','O','Oice','V','Tdp','PH_d','N2Oflx','LEAK_NO3', ...
        'Rh','Rd','Lk','dth','Domain','lateral_diag', ...
        'Final_state','Rh_d','Rd_d','Lk_d','Hourly_tail','PH_bulk_d'};
    S = load(file,wanted{ismember(wanted,who('-file',file))});
    %% A multi-square result of package y keeps no hourly records: the last
    %% hours are in its final state, and every hour's runoff and drainage
    %% in the daily stores and the hours after the last daily step. Its
    %% daily pH is PH_bulk_d, which PH_d repeats row for row.
    daily_store = ~isfield(S,'O');
    if daily_store
        for f = {'O','Oice','V','Tdp'}
            S.(f{1}) = S.Final_state.hourly.(f{1});
        end
        S.PH_d = S.PH_bulk_d;
    end

    %% Assert the depth the matrix configured, not a literal. Pinning 2500
    %% here encoded the depth that diluted DOC per gram of soil and starved
    %% the microbial pools, so the test asserted the defect and would have
    %% failed the corrected runs. What matters is that every scenario uses
    %% the matrix depth, so treatments stay comparable.
    metadata_ok = metadata_ok && S.R.id == k && ...
        strcmp(S.R.name,RM(k).name) && ...
        strcmp(char(S.MOD_PARAM_input.Location),'Kantola') && ...
        S.MOD_PARAM_input.Zbio == RM(k).biogeo_depth_mm;

    if RM(k).use_subdomains
        geometry_ok = geometry_ok && S.N_col == 12 && ...
            abs(sum(S.topo.cell_area_m2)-3600) < 1e-10 && ...
            abs(sum(S.topo.cell_area_m2(strcmp(S.topo.cell_region,'main')))-2500) < 1e-10;
    else
        geometry_ok = geometry_ok && S.N_col == 1 && ...
            abs(S.topo.cell_area_m2-3600) < 1e-10;
    end

    is_main = strcmp(S.topo.cell_region,'main');
    apply_rock = [S.ColCfg.apply_rock].';
    expected_apply = RM(k).apply_main*is_main + RM(k).apply_margin*(~is_main);
    placement_ok = placement_ok && isequal(apply_rock,logical(expected_apply));
    pft = string({S.ColCfg.PFT_OVERRIDE}).';
    if strcmp(RM(k).layout,'margin_prairie')
        placement_ok = placement_ok && all(pft(is_main) == "") && ...
            all(pft(~is_main) == "PFT_Prairie_Margin");
    else
        placement_ok = placement_ok && all(pft == "");
    end

    state_fields = {'P','O','Oice','V','Tdp','PH_d','N2Oflx','LEAK_NO3'};
    for f = 1:numel(state_fields)
        state_ok = state_ok && all(isfinite(S.(state_fields{f})),'all');
    end
    min_NO3 = min(min_NO3,min(S.P(:,32,:),[],'all'));
    state_ok = state_ok && all(S.O >= 0 & S.O <= 1,'all') && ...
        all(S.Oice >= 0 & S.O+S.Oice <= 1+1e-12,'all') && ...
        all(S.PH_d > 0 & S.PH_d < 14,'all') && min_NO3 >= -1e-12;

    area = S.topo.cell_area_m2(:);
    if daily_store
        surface_expected = squeeze((sum(S.Rh_d+S.Rd_d,1,'omitnan') + ...
            S.Hourly_tail.Rh + S.Hourly_tail.Rd)/1000).*area;
        vertical_expected = squeeze((sum(S.Lk_d,1,'omitnan') + ...
            S.Hourly_tail.Lk)/1000).*area;
    else
        surface_expected = squeeze(sum((S.Rh+S.Rd)/1000,1)).*area;
        vertical_expected = squeeze(sum((S.Lk/1000)*S.dth,1)).*area;
    end
    water_boundary_ok = water_boundary_ok && ...
        max(abs(S.Domain.exit_water_surface_m3(:)-surface_expected(:))) < 1e-9 && ...
        max(abs(S.Domain.exit_water_vert_m3(:)-vertical_expected(:))) < 1e-9 && ...
        max(abs(S.Domain.exit_water_total_m3(:)- ...
        surface_expected(:)-vertical_expected(:))) < 1e-9;

    max_heat_residual_relative = max(max_heat_residual_relative, ...
        max(S.lateral_diag.heat_residual_relative,[],'all'));
    max_erw_residual = max(max_erw_residual, ...
        max(S.lateral_diag.erw_mass_abs_residual_umol,[],'all'));
    lateral_ok = lateral_ok && ...
        all(isfinite(S.lateral_diag.erw_mass_abs_residual_umol),'all') && ...
        all(S.lateral_diag.erw_mass_abs_residual_umol >= 0,'all') && ...
        all(isfinite(S.lateral_diag.erw_mass_throughput_umol),'all') && ...
        all(S.lateral_diag.erw_mass_throughput_umol >= 0,'all');
    if ~RM(k).use_subdomains
        lateral_ok = lateral_ok && ...
            ~any(S.lateral_diag.water_buf_to_main,'all') && ...
            ~any(S.lateral_diag.mass_buf_to_main,'all') && ...
            ~any(S.lateral_diag.heat_buf_to_main_J,'all') && ...
            ~any(S.lateral_diag.erw_mass_throughput_umol,'all');
    end
end

lateral_ok = lateral_ok && max_heat_residual_relative < 1e-10 && ...
    max_erw_residual < 1e-8;
check(metadata_ok,'all results retain the Kantola setup and the matrix active depth');
check(geometry_ok,'all subdomain and one-column geometries are correct');
check(placement_ok,'feedstock and prairie placement match the twelve-run matrix');
check(state_ok,'saved water, heat, chemistry and flux states are finite and physical');
check(water_boundary_ok,'surface runoff and bottom leakage close the water exit ledger');
check(lateral_ok,'closed lateral transport conserves heat and ERW mass');

fprintf('  minimum NO3 pool: %.6g gN/m2\n',min_NO3);
fprintf('  maximum scaled heat residual: %.3e\n',max_heat_residual_relative);
fprintf('  maximum ERW residual: %.3e umol\n',max_erw_residual);
fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_subfield_12_smoke_results:Failed','%s',label);
        n_pass = n_pass + 1;
        fprintf('  PASS  %s\n',label);
    end

end
