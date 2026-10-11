function [CEC,diag] = ERW_Resolve_CEC(CEC)
if ~isscalar(CEC) || ~isfinite(CEC) || CEC <= 0
    error('ERW:InvalidCEC', ...
        'A finite positive prescribed CEC is required in cmolc kg-1.');
end
diag = struct('source',"explicit_prescribed_CEC",'estimated',false, ...
    'measured',false,'CEC_cmolc_kg',CEC);
end
