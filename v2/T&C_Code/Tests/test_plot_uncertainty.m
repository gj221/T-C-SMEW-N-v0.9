function test_plot_uncertainty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The per-year uncertainty figure draws the envelope members as boxes,
%% the base set as a filled marker, the sensitivity-only members as open
%% markers with legend entries, and the observations with their standard
%% errors, in the treatment colours.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v2 = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(v2,'plotting'));

out_dir = tempname;
mkdir(out_dir);
cleanup = onCleanup(@() tidy(out_dir));

%% Synthetic summary: 36 members (33 in the envelope, three sensitivity-only
%% members far outside it), the five validation years, three runs.
n_members = 36; years = 2016:2020; n_years = numel(years);
run_ids = [10 11 12]; treatments = {'control','basalt','calcite'};
n_runs = numel(run_ids);
in_envelope = true(1,n_members); in_envelope(34:36) = false;
spread = [0; linspace(-1,1,32).'; 4; -3; 3];         %% member 1 is the base
members = repmat(struct('id',0,'tag','','parameter','','level','','envelope',true, ...
    'KNmaxd',[],'AppeloAlpha',[]),1,n_members);
for k = 1:n_members
    members(k).id = k; members(k).tag = sprintf('m%02d',k);
    members(k).envelope = in_envelope(k);
end
members(34).tag = 'knmaxd_literature'; members(34).parameter = 'KNmaxd';
members(34).level = 'literature'; members(34).KNmaxd = 0.217;
members(35).tag = 'alpha_min'; members(35).parameter = 'AppeloAlpha';
members(35).level = 'min'; members(35).AppeloAlpha = 1.825;
members(36).tag = 'alpha_max'; members(36).parameter = 'AppeloAlpha';
members(36).level = 'max'; members(36).AppeloAlpha = 2.35;
ph_base = [6.2 6.25 6.3 6.15 6.0];
ph_shift = [0 0.35 0.55];                             %% control, basalt, calcite
yearly = struct('years',years);
yearly.pH_season_mean = nan(n_members,n_years,n_runs);
yearly.N2O_cumulative_gN_m2 = nan(n_members,n_years,n_runs);
yearly.NO3_leaching_cumulative_gN_m2 = nan(n_members,n_years,n_runs);
yearly.Ca_release_mol_m2 = nan(n_members,n_years,n_runs);
yearly.Mg_release_mol_m2 = nan(n_members,n_years,n_runs);
for r = 1:n_runs
    growth = (r > 1)*cumsum([0 1 1 1 1]);           %% no feedstock in the control
    yearly.pH_season_mean(:,:,r) = ph_base+ph_shift(r)*(0:n_years-1)/(n_years-1)+0.12*spread;
    yearly.N2O_cumulative_gN_m2(:,:,r) = (0.30-0.02*r)*[0.6 1 1.1 0.9 1]+0.08*spread;
    yearly.NO3_leaching_cumulative_gN_m2(:,:,r) = (2.0-0.1*r)*[0.5 1 1.2 0.8 1]+0.4*spread;
    yearly.Ca_release_mol_m2(:,:,r) = 2.5*growth.*(1+0.3*spread);
    yearly.Mg_release_mol_m2(:,:,r) = 1.2*growth.*(1+0.3*spread);
end
yearly.pH_season_signal = yearly.pH_season_mean(:,:,2)-yearly.pH_season_mean(:,:,1);
yearly.CDR_potential_tCO2_ha = 2*(yearly.Ca_release_mol_m2+yearly.Mg_release_mol_m2)*44e-6*1e4;
obs = struct('years',years);
obs.pH.control_0_10 = struct('value',[6.203 6.278 6.313 6.131 5.976],'sem',[0.063 0.076 0.051 0.073 0.054]);
obs.pH.control_10_30 = struct('value',[6.020 6.285 6.315 6.128 6.244],'sem',[0.113 0.083 0.081 0.071 0.078]);
obs.pH.basalt_0_10 = struct('value',[6.127 6.600 6.663 6.919 6.918],'sem',[0.074 0.060 0.042 NaN 0.031]);
obs.pH.basalt_10_30 = struct('value',[6.106 6.403 6.565 6.477 6.509],'sem',[0.071 0.068 0.056 0.048 0.048]);
obs.CDR_potential_tCO2_ha = struct('value',[0 3.82 6.75 11.49 10.51],'sem',[0 1.83 2.51 3.26 3.72]);
obs.Ca_release_mmol_m2 = struct('value',[0 2278 6290 9937 6910],'sem',[0 388 366 353 693]);
obs.Mg_release_mmol_m2 = struct('value',[0 2094 1412 3174 5106],'sem',[0 802 1871 2503 2192]);
summary_file = fullfile(out_dir,'uncertainty_summary.mat');
save(summary_file,'run_ids','treatments','yearly','obs','in_envelope','members');
stale = rmfield(load(summary_file),'in_envelope');
stale_file = fullfile(out_dir,'stale_summary.mat');
save(stale_file,'-struct','stale');
assert_error(@() Plot_Uncertainty(stale_file,out_dir,'off'), ...
    'Plot_Uncertainty:BadSummary');

assert_error(@() Plot_Uncertainty(fullfile(out_dir,'absent.mat'),out_dir,'off'), ...
    'Plot_Uncertainty:MissingSummary');

figure_dir = fullfile(out_dir,'figures');
[fig,plot_file] = Plot_Uncertainty(summary_file,figure_dir,'off');
assert(numel(fig) == 3 && all(isgraphics(fig,'figure')) && ...
    all(strcmp({fig.Visible},'off')),'Three invisible figures are expected.');
assert(numel(plot_file) == 3 && all(cellfun(@isfile,plot_file)), ...
    'Each figure must be written as PNG.');
for k = 1:3
    [folder,name] = fileparts(plot_file{k});
    assert(strcmp(folder,figure_dir) && isfile(fullfile(folder,[name '.fig'])), ...
        'Each PNG must have a FIG companion in the figure directory.');
    assert(isequal(fig(k).Color,[1 1 1]),'Figures must have a white background.');
    assert(isscalar(findall(fig(k),'Type','legend')),'One legend per figure.');
end
treatment_colors = [0 0 0; 0 0.25 0.60; 0.35 0.70 0.90];

%% pH figure: treatment boxes with observations by depth, then the signal.
axes_ph = flipud(findall(fig(1),'Type','axes'));
assert(numel(axes_ph) == 2,'The pH figure has two panels.');
boxes = findall(axes_ph(1),'Type','boxchart');
assert(numel(boxes) == 3,'One box series per treatment in the pH panel.');
for treatment = 1:3
    bx = boxes(arrayfun(@(h) isequal(h.BoxFaceColor,treatment_colors(treatment,:)),boxes));
    assert(isscalar(bx) && strcmp(bx.MarkerStyle,'none') && ...
        numel(bx.YData) == nnz(in_envelope)*n_years, ...
        'Each treatment box series carries every envelope member and year in its colour.');
    assert(max(bx.YData) <= max(yearly.pH_season_mean(in_envelope,:,treatment),[],'all')+1e-12, ...
        'Sensitivity-only members must stay out of the boxes.');
    assert(isequal(bx.WhiskerLineColor,treatment_colors(treatment,:)), ...
        'Whiskers take the treatment colour.');
end
assert(numel(findall(axes_ph(1),'Type','line','Tag','BaseMember')) == 3, ...
    'The base member is marked once per treatment.');
sensitivity_lines = findall(axes_ph(1),'Type','line','Tag','SensitivityMember');
assert(numel(sensitivity_lines) == 3*nnz(~in_envelope) && ...
    all(arrayfun(@(h) isequal(h.MarkerFaceColor,[1 1 1]),sensitivity_lines)) && ...
    numel(unique({sensitivity_lines.Marker})) == nnz(~in_envelope), ...
    'Each sensitivity-only member is an open marker of its own shape, once per treatment.');
knmaxd_marker = sensitivity_lines(arrayfun(@(h) isequal(h.MarkerEdgeColor,treatment_colors(2,:)) && ...
    all(abs(h.YData-yearly.pH_season_mean(34,:,2)) < 1e-12),sensitivity_lines));
assert(isscalar(knmaxd_marker),'The KNmaxd literature member is drawn at its own values.');
lgd = findall(fig(1),'Type','legend');
assert(any(contains(lgd.String,'KNmaxd = 0.217 d^{-1} (literature)')) && ...
    any(contains(lgd.String,'\alpha = 1.825')) && any(contains(lgd.String,'\alpha = 2.35')), ...
    'The legend names the sensitivity-only members.');
ph_obs = findall(axes_ph(1),'Type','errorbar');
assert(numel(ph_obs) == 4,'Control and basalt pH at two depths carry error bars.');
for h = ph_obs.'
    assert(isequal(h.Color,[0 0 0]),'Observations are black.');
end
open_marker = findall(axes_ph(1),'Type','errorbar','Tag','Observed_basalt_pH_10_30');
filled_marker = findall(axes_ph(1),'Type','errorbar','Tag','Observed_basalt_pH_0_10');
assert(isequal(open_marker.MarkerFaceColor,[1 1 1]) && ...
    isequal(filled_marker.MarkerFaceColor,[0 0 0]), ...
    '0-10 cm observations are filled, 10-30 cm open.');
assert(isequaln(filled_marker.YNegativeDelta(:),obs.pH.basalt_0_10.sem(:)) && ...
    isequal(filled_marker.YData(:),obs.pH.basalt_0_10.value(:)), ...
    'pH error bars are the observed SEM.');
signal_boxes = findall(axes_ph(2),'Type','boxchart');
assert(isscalar(signal_boxes) && isequal(signal_boxes.BoxFaceColor,treatment_colors(2,:)) && ...
    numel(signal_boxes.YData) == nnz(in_envelope)*n_years, ...
    'The basalt-minus-control signal is one box series in the basalt colour over the envelope.');
assert(numel(findall(axes_ph(2),'Type','line','Tag','SensitivityMember')) == nnz(~in_envelope), ...
    'The signal panel marks each sensitivity-only member once.');
signal_obs = findall(axes_ph(2),'Type','errorbar');
assert(numel(signal_obs) == 2,'The observed signal is drawn for both depth intervals.');
expected_signal = obs.pH.basalt_0_10.value-obs.pH.control_0_10.value;
signal_0_10 = findall(axes_ph(2),'Type','errorbar','Tag','Observed_signal_pH_0_10');
assert(max(abs(signal_0_10.YData(:)-expected_signal(:))) < 1e-12 && ...
    abs(signal_0_10.YNegativeDelta(1)-hypot(0.074,0.063)) < 1e-12, ...
    'The observed signal carries the propagated SEM.');

%% Nitrogen figure: three treatments per panel, no observations.
axes_n = findall(fig(2),'Type','axes');
assert(numel(axes_n) == 2,'The nitrogen figure has two panels.');
for ax = axes_n.'
    assert(numel(findall(ax,'Type','boxchart')) == 3 && ...
        isempty(findall(ax,'Type','errorbar')) && ...
        numel(findall(ax,'Type','line','Tag','SensitivityMember')) == 3*nnz(~in_envelope), ...
        'Nitrogen panels carry three treatment box series, the sensitivity markers and no observations.');
end

%% Weathering figure: amended treatments only, Dataset S06 in black.
axes_w = findall(fig(3),'Type','axes');
assert(numel(axes_w) == 3,'The weathering figure has three panels.');
for ax = axes_w.'
    boxes = findall(ax,'Type','boxchart');
    assert(numel(boxes) == 2 && ...
        any(arrayfun(@(h) isequal(h.BoxFaceColor,treatment_colors(2,:)),boxes)) && ...
        any(arrayfun(@(h) isequal(h.BoxFaceColor,treatment_colors(3,:)),boxes)), ...
        'Release and CDR panels show basalt and calcite only.');
    assert(isscalar(findall(ax,'Type','errorbar')), ...
        'Each release and CDR panel carries one observed series.');
    assert(numel(findall(ax,'Type','line','Tag','SensitivityMember')) == 2*nnz(~in_envelope), ...
        'Release and CDR panels mark the sensitivity-only members for both amended treatments.');
end
ca_obs = findall(fig(3),'Type','errorbar','Tag','Observed_Ca_release_mol_m2');
assert(max(abs(ca_obs.YData(:)-1e-3*obs.Ca_release_mmol_m2.value(:))) < 1e-12 && ...
    max(abs(ca_obs.YPositiveDelta(:)-1e-3*obs.Ca_release_mmol_m2.sem(:))) < 1e-12, ...
    'Observed release is converted from mmol to mol m-2.');
cdr_obs = findall(fig(3),'Type','errorbar','Tag','Observed_CDR_potential_tCO2_ha');
assert(isequal(cdr_obs.YData(:),obs.CDR_potential_tCO2_ha.value(:)), ...
    'Observed CDR potential is drawn as reported.');
for ax = [axes_ph.' axes_n.' axes_w.']
    assert(isequal(ax.XTick,years),'Year ticks label every panel.');
    assert(~isempty(ax.YLabel.String),'Every panel names its units.');
end

close(fig);
fprintf('PASS: the uncertainty figure draws members, base set and observations with SEM.\n');
end

function tidy(out_dir)
close(findall(groot,'Type','figure'));
if isfolder(out_dir); rmdir(out_dir,'s'); end
end

function assert_error(action,identifier)
try
    action();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected %s but received %s.',identifier,ME.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end
