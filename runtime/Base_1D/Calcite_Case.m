function [calcite_case,suffix] = Calcite_Case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TC_CALCITE_CASE names an ERW workbook row for the calcite treatments
%% only, such as 8, the agronomic 5 t/ha rate; unset is the standard
%% campaign. The suffix keeps the case campaign's tag, and so its results,
%% logs and stamp, apart from the standard one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw = strtrim(getenv('TC_CALCITE_CASE'));
calcite_case = [];
suffix = '';
if isempty(raw); return; end
calcite_case = str2double(raw);
if ~(isfinite(calcite_case) && calcite_case == fix(calcite_case) && calcite_case >= 1)
    error('Bouchet_Kantola_CEC13:BadCalciteCase', ...
        'TC_CALCITE_CASE must be unset or a positive integer, not ''%s''.',raw);
end
suffix = sprintf('_case%d',calcite_case);
end
