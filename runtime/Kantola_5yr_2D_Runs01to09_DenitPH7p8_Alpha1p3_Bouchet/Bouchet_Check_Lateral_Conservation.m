function report = Bouchet_Check_Lateral_Conservation(result_file,error_id_prefix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Confirm the lateral scheme conserved mass and energy over a whole run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test_lateral_routing proves conservation for a single step in isolation.
%% A campaign run takes 43824 steps across twelve columns, and nothing was
%% checking that the residuals stayed small once accumulated. The driver
%% already records them; this reads them back and decides.
%%
%% The subsurface lateral perimeter is closed, so transport only
%% redistributes: summed over every cell and layer the aqueous mass change
%% from lateral exchange must be zero, and so must the energy change. The
%% driver stores those sums per day as
%%
%%   lateral.erw_mass_residual_umol   [NNd x 4]  Na, Al, An, DIC
%%   lateral.heat_residual_J          [NNd x 1]
%%   lateral.heat_residual_relative   [NNd x 1]  already normalised
%%
%% Residuals are scaled by the mass actually transported, because an
%% absolute umol figure means nothing without knowing how much moved.
%%
%% This runs after the result is on disk, so a failure flags the run
%% without destroying it: the MAT stays for inspection.
%%
%% INPUTS
%%   result_file       saved campaign MAT
%%   error_id_prefix   identifier prefix for the calling campaign
%% OUTPUT
%%   report            struct of the residual measures, for the log

if nargin < 2 || isempty(error_id_prefix)
    error_id_prefix = 'Bouchet';
end
error_id_prefix = char(error_id_prefix);

%% Tolerances. A conservative scheme should land near machine precision;
%% these are set far looser so that only a real leak trips them, not the
%% round-off of 43824 accumulated steps.
warn_relative = 1e-9;
fail_relative = 1e-6;

loaded = load(result_file,'Domain');
if ~isfield(loaded,'Domain') || ~isfield(loaded.Domain,'lateral')
    error([error_id_prefix ':NoLateralDiagnostics'], ...
        'The saved result carries no lateral diagnostics: %s',result_file);
end
L = loaded.Domain.lateral;

report = struct();
failures = {};

%% ---------------------------------------------------------------------
%% Aqueous mass: Na, Al, An, DIC
%% ---------------------------------------------------------------------
species = {'Na','Al','An','DIC'};
if isfield(L,'erw_mass_residual_umol') && ~isempty(L.erw_mass_residual_umol)
    residual = L.erw_mass_residual_umol;
    if isfield(L,'erw_mass_buf_to_main_umol') && ...
            ~isempty(L.erw_mass_buf_to_main_umol)
        moved = L.erw_mass_buf_to_main_umol;
    else
        moved = zeros(size(residual));
    end
    for s = 1:min(numel(species),size(residual,2))
        r = sum(residual(:,s),'omitnan');
        scale = max(1,sum(abs(moved(:,min(s,size(moved,2)))),'omitnan'));
        rel = abs(r)/scale;
        report.(['mass_' species{s} '_residual_umol']) = r;
        report.(['mass_' species{s} '_relative'])      = rel;
        fprintf(['  lateral %-4s residual %+.6g umol of %.6g moved ' ...
            '(relative %.3e)\n'],species{s},r,scale,rel);
        if ~isfinite(rel) || rel > fail_relative
            failures{end+1} = sprintf('%s relative residual %.3e', ...
                species{s},rel); %#ok<AGROW>
        elseif rel > warn_relative
            warning([error_id_prefix ':LateralMassDrift'], ...
                'Lateral %s residual is %.3e of transported mass.', ...
                species{s},rel);
        end
    end
else
    failures{end+1} = 'aqueous mass residuals absent';
end

%% ---------------------------------------------------------------------
%% Energy
%% ---------------------------------------------------------------------
if isfield(L,'heat_residual_relative') && ~isempty(L.heat_residual_relative)
    heat_rel = max(abs(L.heat_residual_relative),[],'omitnan');
    report.heat_residual_relative_max = heat_rel;
    if isfield(L,'heat_residual_J')
        report.heat_residual_J_total = sum(L.heat_residual_J,'omitnan');
        fprintf('  lateral heat residual %+.6g J, max relative %.3e\n', ...
            report.heat_residual_J_total,heat_rel);
    else
        fprintf('  lateral heat max relative residual %.3e\n',heat_rel);
    end
    if ~isfinite(heat_rel) || heat_rel > fail_relative
        failures{end+1} = sprintf('heat relative residual %.3e',heat_rel);
    elseif heat_rel > warn_relative
        warning([error_id_prefix ':LateralHeatDrift'], ...
            'Lateral energy residual reached %.3e.',heat_rel);
    end
else
    failures{end+1} = 'energy residuals absent';
end

%% ---------------------------------------------------------------------
%% Adaptive sub-stepping. Saturation means the solver was straining and
%% the fluxes may have been limited rather than resolved.
%% ---------------------------------------------------------------------
for f = {'sub_steps_max','heat_sub_steps_max'}
    if isfield(L,f{1}) && ~isempty(L.(f{1}))
        m = max(L.(f{1}),[],'omitnan');
        report.(f{1}) = m;
        fprintf('  %-19s peak %g\n',f{1},m);
    end
end

%% ---------------------------------------------------------------------
%% Boundary water, reported for the run log rather than tested: the plot
%% loses water vertically by design, so there is no closure to assert.
%% ---------------------------------------------------------------------
D = loaded.Domain;
if isfield(D,'exit_water_total_m3')
    report.exit_water_total_m3 = sum(D.exit_water_total_m3,'omitnan');
    fprintf('  water leaving the plot %.4f m3 (vertical + surface)\n', ...
        report.exit_water_total_m3);
end
if isfield(L,'water_buf_to_main')
    report.water_buf_to_main_m3 = sum(L.water_buf_to_main,'omitnan');
    fprintf('  net lateral water margin->field %.4f m3\n', ...
        report.water_buf_to_main_m3);
end

if ~isempty(failures)
    error([error_id_prefix ':LateralConservation'], ...
        ['The lateral scheme did not conserve over the run: %s. The ' ...
         'result MAT is on disk for inspection.'],strjoin(failures,'; '));
end
fprintf('  lateral conservation OK\n');

end
