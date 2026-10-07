function logK = ERW_AlOH3_LogK(use_amorphous)
if ~islogical(use_amorphous) || ~isscalar(use_amorphous)
    error('ERW_AlOH3_LogK:InvalidPhase', ...
        'The amorphous-phase selector must be a logical scalar.');
end
if use_amorphous
    logK = 10.8;
else
    logK = 8.04;
end
end
