function [kinetic_mass_fraction, direct_carbonate_mass_fraction, rock_f_kinetic] = ...
    ERW_Feedstock_Pathway(mineral, rock_f_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Route an applied feedstock between the primary and secondary pools.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The model keeps two distinct calcite inventories and ERW_Total_Calcite_Solid
%% sums them:
%%
%%   primary    calcite grains in M_min, dissolving under their own kinetic
%%              rate law and particle size distribution
%%   secondary  pedogenic carbonate in the CaCO3 pool, precipitating and
%%              dissolving with the Omega_CaCO3 saturation state
%%
%% kinetic_mass_fraction sends applied mass to the first, via M_rock_add.
%% direct_carbonate_mass_fraction sends it to the second, via
%% M_rock_add_CaCO3 and CaCO3(1) = CaCO3_in + M_rock_add_CaCO3(1).
%%
%% Applied feedstock is always primary: it arrives as graded grains whose
%% dissolution is surface-area limited, which is what the kinetic pool
%% represents. This holds for the calcite in a basalt feedstock as much as
%% for a pure calcite amendment, so both return kinetic = 1. Sending
%% feedstock to the secondary pool instead would place it under a
%% saturation law it does not obey, and would conflate rock the operator
%% applied with carbonate the soil produced. That distinction matters for
%% the carbon accounting, because forming secondary carbonate releases CO2
%% while weathering primary carbonate consumes it.
%%
%% The carbonate return is therefore zero for every feedstock currently
%% configured, not unused: it is how a directly seeded carbonate would be
%% represented, such as finely ground agricultural lime that dissolves on
%% contact rather than over seasons. test_pH_conversion and
%% test_erw_application_event_alignment pin both wirings in place.

if iscell(mineral) && isscalar(mineral); mineral = mineral{1}; end
if iscell(rock_f_in) && isscalar(rock_f_in); rock_f_in = rock_f_in{1}; end
if ((ischar(mineral) && isrow(mineral)) || ...
        (isstring(mineral) && isscalar(mineral))) && ...
        startsWith(strtrim(string(mineral)),"[")
    names = strtrim(string(jsondecode(char(mineral))));
else
    names = strtrim(string(cellstr(mineral)));
end
if ischar(rock_f_in) || isstring(rock_f_in)
    fractions = double(jsondecode(char(rock_f_in)));
else
    fractions = rock_f_in;
end
rock_f_kinetic = fractions;
if isempty(fractions) || (isnumeric(fractions) && all(isnan(fractions(:))))
    kinetic_mass_fraction = 0;
    direct_carbonate_mass_fraction = 0;
    return
end

fractions = fractions(:).';
if numel(names) ~= numel(fractions) || any(~isfinite(fractions)) || ...
        any(fractions < 0) || sum(fractions) > 1 + 1e-10
    error('ERW:invalidFeedstockMineralogy', ...
        'Mineral names and finite nonnegative mass fractions must match and sum to at most one.');
end

unsupported = ismember(lower(names),["aragonite","dolomite","magnesite"]);
if any(unsupported)
    error('ERW:unsupportedCarbonate', ...
        'Only calcite has a carbonate mass and stoichiometry pathway in this model.');
end

%% Every configured feedstock is graded rock, so all of its mass — calcite
%% included — dissolves kinetically from the primary pool.
kinetic_mass_fraction = 1;
direct_carbonate_mass_fraction = 0;
rock_f_kinetic = fractions;

end
