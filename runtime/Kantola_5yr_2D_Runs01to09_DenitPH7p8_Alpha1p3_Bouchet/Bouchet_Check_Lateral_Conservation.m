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
%%   lateral.erw_mass_abs_residual_umol  [NNd x 4]  Na, Al, An, DIC
%%   lateral.erw_mass_throughput_umol    [NNd x 4]
%%   lateral.erw_IC_atmosphere_abs_residual_umol  [NNd x 1]
%%   lateral.erw_IC_atmosphere_throughput_umol    [NNd x 1]
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
if isfield(L,'erw_mass_abs_residual_umol') && ...
        ~isempty(L.erw_mass_abs_residual_umol) && ...
        isfield(L,'erw_mass_throughput_umol') && ...
        ~isempty(L.erw_mass_throughput_umol)
    residual = L.erw_mass_abs_residual_umol;
    moved = L.erw_mass_throughput_umol;
    if ~isequal(size(residual),size(moved)) || size(residual,2) ~= numel(species)
        failures{end+1} = 'aqueous mass diagnostics have invalid dimensions';
    else
        for s = 1:numel(species)
            residual_s = residual(:,s);
            moved_s = moved(:,s);
            physical = all(isfinite(residual_s)) && all(residual_s >= 0) && ...
                all(isfinite(moved_s)) && all(moved_s >= 0);
            r = sum(residual_s);
            throughput = sum(moved_s);
            rel = r/max(1,throughput);
            report.(['mass_' species{s} '_residual_umol']) = r;
            report.(['mass_' species{s} '_throughput_umol']) = throughput;
            report.(['mass_' species{s} '_relative'])      = rel;
            fprintf(['  lateral %-4s absolute residual %.6g umol of %.6g moved ' ...
                '(relative %.3e)\n'],species{s},r,throughput,rel);
            if ~physical || ~isfinite(rel) || rel > fail_relative
                failures{end+1} = sprintf('%s relative residual %.3e', ...
                    species{s},rel); %#ok<AGROW>
            elseif rel > warn_relative
                warning([error_id_prefix ':LateralMassDrift'], ...
                    'Lateral %s residual is %.3e of transported mass.', ...
                    species{s},rel);
            end
        end
    end
else
    failures{end+1} = 'uncancelled aqueous mass diagnostics absent';
end

%% ---------------------------------------------------------------------
%% Soil-air piston exchange
%% ---------------------------------------------------------------------
if isfield(L,'erw_IC_atmosphere_abs_residual_umol') && ...
        ~isempty(L.erw_IC_atmosphere_abs_residual_umol) && ...
        isfield(L,'erw_IC_atmosphere_throughput_umol') && ...
        ~isempty(L.erw_IC_atmosphere_throughput_umol)
    residual = L.erw_IC_atmosphere_abs_residual_umol;
    exchanged = L.erw_IC_atmosphere_throughput_umol;
    if ~isvector(residual) || ~isvector(exchanged) || numel(residual) ~= numel(exchanged)
        failures{end+1} = 'ERW atmosphere diagnostics have invalid dimensions';
    else
        residual = residual(:);
        exchanged = exchanged(:);
        physical = all(isfinite(residual)) && all(residual >= 0) && ...
            all(isfinite(exchanged)) && all(exchanged >= 0);
        r = sum(residual);
        throughput = sum(exchanged);
        rel = r/max(1,throughput);
        report.IC_atmosphere_residual_umol = r;
        report.IC_atmosphere_throughput_umol = throughput;
        report.IC_atmosphere_relative = rel;
        fprintf(['  ERW atmosphere absolute residual %.6g umol of %.6g exchanged ' ...
            '(relative %.3e)\n'],r,throughput,rel);
        if ~physical || ~isfinite(rel) || rel > fail_relative
            failures{end+1} = sprintf('atmosphere IC relative residual %.3e',rel); %#ok<AGROW>
        elseif rel > warn_relative
            warning([error_id_prefix ':LateralMassDrift'], ...
                'ERW atmosphere IC residual is %.3e of exchanged mass.',rel);
        end
    end
else
    failures{end+1} = 'ERW atmosphere carbon diagnostics absent';
end

%% ---------------------------------------------------------------------
%% Energy
%% ---------------------------------------------------------------------
if isfield(L,'heat_residual_relative') && ~isempty(L.heat_residual_relative)
    heat_series = L.heat_residual_relative;
    heat_physical = all(isfinite(heat_series)) && all(heat_series >= 0);
    heat_rel = max(abs(heat_series));
    report.heat_residual_relative_max = heat_rel;
    if isfield(L,'heat_residual_J')
        heat_J_physical = all(isfinite(L.heat_residual_J));
        report.heat_residual_J_total = sum(L.heat_residual_J);
        fprintf('  lateral heat residual %+.6g J, max relative %.3e\n', ...
            report.heat_residual_J_total,heat_rel);
        if ~heat_J_physical
            failures{end+1} = 'heat residual contains a nonfinite value';
        end
    else
        fprintf('  lateral heat max relative residual %.3e\n',heat_rel);
    end
    if ~heat_physical || ~isfinite(heat_rel) || heat_rel > fail_relative
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
