function [MOD_base,ERW_treatment,ERW_control] = ...
    Select_Subfield_Parameter_Rows(MOD_all,ERW_all,treatment_case)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hold the T&C setup fixed and select only the ERW feedstock row.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CONTROL_CASE = 5;
MOD_base = MOD_all(MOD_all.Number == CONTROL_CASE,:);
ERW_control = ERW_all(ERW_all.Number == CONTROL_CASE,:);
ERW_treatment = ERW_all(ERW_all.Number == treatment_case,:);
if height(MOD_base) ~= 1 || height(ERW_control) ~= 1 || ...
        height(ERW_treatment) ~= 1
    error('Select_Subfield_Parameter_Rows:MissingCase', ...
        'The control and selected treatment must each have one workbook row.');
end

feedstock_fields = {'Number','ExperimentalDescription','rock_dist_data', ...
    'rock_filepath','rock_type','MgCO3_in','CaCO3_in','M_rock_in', ...
    'mineral','rock_f_in','rho_rock','d_in','psd_in','SSA_in', ...
    'diss_f_original','diss_f','M_rock_event'};
names = ERW_all.Properties.VariableNames;
for k = 1:numel(names)
    name = names{k};
    if ismember(name,feedstock_fields); continue; end
    if ~isequaln(ERW_treatment.(name),ERW_control.(name))
        error('Select_Subfield_Parameter_Rows:NonFeedstockDifference', ...
            'ERW field %s differs from control but is not a feedstock field.',name);
    end
end

end
