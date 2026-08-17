function BiogeoPar = Sensitivity_Overrides(BiogeoPar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply sweep overrides to the biogeochemistry parameters.
%%
%% This file is NOT part of the model. Biogeochemistry_Parameter calls it
%% only when this folder has been put on the path, which Sensitivity_Run
%% does for its children and nothing else does. Take v2/Sensitivity off the
%% path, or delete it, and the model runs its own parameters unchanged.
%%
%% Values arrive as TC_BGP_<name> environment variables. Environment is
%% per-process, so parallel sets cannot contaminate one another, which a
%% shared file or a global could not guarantee.
%%
%% A malformed value raises rather than falling back: a sweep that silently
%% runs the baseline and reports it as a treatment is the worst outcome.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for name = Sensitivity_Parameters()
    raw = getenv(['TC_BGP_' name{1}]);
    if isempty(raw); continue; end
    val = str2double(raw);
    if ~isscalar(val) || ~isfinite(val)
        error('Sensitivity_Overrides:BadValue', ...
            'TC_BGP_%s = "%s" is not a finite scalar.',name{1},raw);
    end
    if ~isfield(BiogeoPar,name{1})
        error('Sensitivity_Overrides:UnknownParameter', ...
            ['%s is not a field of BiogeoPar. Overriding it would create ' ...
             'a parameter the model never reads.'],name{1});
    end
    BiogeoPar.(name{1}) = val;
end

end
