function [dCEC_dpH, CEC_ref_pH, diag] = ...
    ERW_Resolve_CEC_pH_Sensitivity(bulk_pH_params, scheme_str)
%%% Resolve the optional empirical TC_GT variable-charge sensitivity.
%%%
%%% SCEPTER assigns a fixed mass-specific CEC to each exchanger solid.  Its
%%% aggregate cell capacity can still change if those solid masses change.
%%% This lumped TC Appelo implementation has no evolving solid inventory, so
%%% it keeps the prescribed soil CEC fixed and represents the pH response
%%% through H-site-dependent apparent selectivity.  A separate positive
%%% dCEC_dpH is therefore rejected in Appelo mode.
%%%
%%% The sensitivity is declared mass-specific, in cmolc kg-1 per pH unit, so
%%% that it stays a property of the soil rather than of the lumped
%%% biogeochemical layer thickness.  The biogeochemistry is vertically
%%% lumped into one layer of depth Zbio, so there is a single exchanger and
%%% a single bulk sensitivity; ERW_function converts to the areal umol_c
%%% basis with exactly the factor that built CEC_tot, which keeps the
%%% parameter invariant if Zbio is ever changed.
%%%
%%% bulk_pH_params.dCEC_dpH may be given as
%%%   - a scalar             an empirical TC_GT sensitivity; Appelo accepts
%%%                          only zero
%%%   - a struct with fields the value for the active exchange scheme is
%%%     .appelo, .TC_GT      selected; a scheme absent from the struct gets
%%%                          zero, i.e. fixed CEC for that scheme
%%% The appelo entry is accepted only when it is zero, so older parameter
%%% structures remain readable.  Absent, empty or zero reproduces the fixed
%%% CEC behaviour bit for bit.  Positive values remain an empirical TC_GT
%%% extension and require an explicit reference pH.

dCEC_dpH = 0;
CEC_ref_pH = NaN;
diag = struct('source',"legacy_constant_CEC",'scheme',string(scheme_str), ...
    'dCEC_dpH_cmolc_kg',0,'CEC_ref_pH',CEC_ref_pH,'pH_dependent',false);

if ~ischar(scheme_str) || ~any(strcmpi(scheme_str,{'TC_GT','appelo'}))
    error('ERW:invalidParam', ...
        'Exchange scheme must be TC_GT or appelo to resolve dCEC_dpH.');
end
if ~isstruct(bulk_pH_params) || ~isfield(bulk_pH_params,'dCEC_dpH') || ...
        isempty(bulk_pH_params.dCEC_dpH)
    return
end

spec = bulk_pH_params.dCEC_dpH;
if isstruct(spec)
    %%% Scheme-keyed form: pick the entry for the scheme actually in use.
    known = {'appelo','TC_GT'};
    extra = setdiff(fieldnames(spec), known);
    if ~isempty(extra)
        error('ERW:invalidParam', ...
            'dCEC_dpH struct accepts only the fields appelo and TC_GT; found %s.', ...
            strjoin(extra(:).', ', '));
    end
    match = known{strcmpi(scheme_str, known)};
    if isfield(spec, match) && ~isempty(spec.(match))
        dCEC_dpH = spec.(match);
        diag.source = "scheme_keyed_dCEC_dpH";
    else
        %%% Scheme carries no entry: constant CEC, as before.
        diag.source = "scheme_absent_from_dCEC_dpH";
        return
    end
else
    dCEC_dpH = spec;
    diag.source = "scalar_dCEC_dpH";
end

if ~(isscalar(dCEC_dpH) && isfinite(dCEC_dpH) && dCEC_dpH >= 0)
    error('ERW:invalidParam', ...
        'dCEC_dpH must be a nonnegative finite scalar in cmolc kg-1 per pH unit.');
end

if dCEC_dpH == 0
    diag.dCEC_dpH_cmolc_kg = dCEC_dpH;
    return
end

if strcmpi(scheme_str,'appelo')
    error('ERW:scepterFixedCEC', ...
        ['TC Appelo mode keeps the lumped soil CEC fixed; pH dependence is ' ...
         'represented by the H-site-dependent exchange selectivity.']);
end

if isfield(bulk_pH_params,'CEC_ref_pH') && ~isempty(bulk_pH_params.CEC_ref_pH)
    CEC_ref_pH = bulk_pH_params.CEC_ref_pH;
    if ~(isscalar(CEC_ref_pH) && isfinite(CEC_ref_pH) && ...
            CEC_ref_pH >= 3 && CEC_ref_pH <= 10)
        error('ERW:invalidParam', ...
            'CEC_ref_pH must be a finite scalar between 3 and 10.');
    end
else
    error('ERW:missingParam', ...
        'Positive dCEC_dpH requires CEC_ref_pH on the pore-water pH basis.');
end

diag.dCEC_dpH_cmolc_kg = dCEC_dpH;
diag.CEC_ref_pH = CEC_ref_pH;
diag.pH_dependent = dCEC_dpH > 0;
end
