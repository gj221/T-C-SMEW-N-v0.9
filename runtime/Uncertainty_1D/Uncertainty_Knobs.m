function K = Uncertainty_Knobs(member,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The parameter values a set runs with: what the member gives, the
%% campaign's set P (Campaign_Parameters) for the rest. A set that moves
%% alpha slides the exchange constants by the campaign's rule, so the
%% seed's measured 80 % base saturation holds at every alpha, and gives
%% them to twelve significant figures. K.options is
%% the Run_One_Subfield option list; for the original it is exactly what
%% the Base and SSP tasks pass. K.SaturationExponent is the value the
%% driver records, 1 when the set leaves it unset.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K.KNmaxd = P.KNmaxd;
if ~isempty(member.KNmaxd); K.KNmaxd = member.KNmaxd; end
K.AppeloAlpha = P.AppeloAlpha;
K.AppeloLogKH = P.AppeloLogKH;
if ~isempty(member.AppeloAlpha)
    K.AppeloAlpha = member.AppeloAlpha;
    K.AppeloLogKH = P.reference_logKH + P.charge*(K.AppeloAlpha-P.reference_alpha)*P.seed_f_H;
    %% Twelve significant figures is the precision a result name carries. A
    %% constant that does not print exactly at that precision adds its
    %% hexadecimal form to the name (number_tag in Run_One_Subfield), and at
    %% alpha 6 two of the five would: the longest 35-year 2-D name would
    %% then pass the 255 characters a file name may have, and fail only
    %% when the finished run is saved.
    K.AppeloLogKH = arrayfun(@(v) str2double(sprintf('%.12g',v)),K.AppeloLogKH);
end
K.Kbot = P.Kbot;
K.BiogeoOverrides = member.BiogeoOverrides;
K.DissolutionFactor = member.DissolutionFactor;
K.SaturationExponent = 1;
if ~isempty(member.SaturationExponent); K.SaturationExponent = member.SaturationExponent; end

K.options = {'AppeloAlpha',K.AppeloAlpha,'AppeloLogKH',K.AppeloLogKH, ...
    'KNmaxd',K.KNmaxd,'Kbot',K.Kbot};
if ~isempty(member.BiogeoOverrides)
    K.options = [K.options {'BiogeoOverrides',member.BiogeoOverrides}];
end
if ~isempty(member.DissolutionFactor)
    K.options = [K.options {'DissolutionFactor',member.DissolutionFactor}];
end
if ~isempty(member.SaturationExponent)
    K.options = [K.options {'SaturationExponent',member.SaturationExponent}];
end
end
