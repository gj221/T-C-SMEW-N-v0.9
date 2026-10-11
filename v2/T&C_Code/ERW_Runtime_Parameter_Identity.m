function identity = ERW_Runtime_Parameter_Identity(T)
%%% Deterministic identity for one effective ERW parameter-table row.

if ~istable(T) || height(T) ~= 1
    error('ERW:checkpointRuntimeTable', ...
        'The effective ERW runtime parameter must be a one-row table.');
end

names = T.Properties.VariableNames;
parts = cell(1,numel(names));
for k = 1:numel(names)
    value = T.(names{k});
    if iscell(value)
        if ~isscalar(value)
            error('ERW:checkpointRuntimeTable', ...
                'Cell-valued runtime parameters must contain one value.');
        end
        value = value{1};
    end
    parts{k} = [encode_text(names{k}) '=' encode_value(value)];
end
identity = strjoin(parts,';');

end

function text = encode_value(value)
if isnumeric(value)
    if ~isreal(value)
        error('ERW:checkpointRuntimeTable', ...
            'Complex runtime parameters are not supported.');
    end
    if isa(value,'double')
        values = arrayfun(@(x) sprintf('%.17g',x),value(:).', ...
            'UniformOutput',false);
    elseif isa(value,'single')
        values = arrayfun(@(x) sprintf('%.9g',x),value(:).', ...
            'UniformOutput',false);
    elseif isinteger(value)
        values = arrayfun(@(x) sprintf('%d',x),value(:).', ...
            'UniformOutput',false);
    else
        error('ERW:checkpointRuntimeTable', ...
            'Unsupported numeric runtime-parameter type %s.',class(value));
    end
    text = sprintf('%s%s:%s',class(value),mat2str(size(value)), ...
        strjoin(values,','));
elseif islogical(value)
    values = arrayfun(@(x) sprintf('%d',x),value(:).', ...
        'UniformOutput',false);
    text = sprintf('logical%s:%s',mat2str(size(value)),strjoin(values,','));
elseif isstring(value)
    if ~isscalar(value)
        error('ERW:checkpointRuntimeTable', ...
            'String runtime parameters must be scalar.');
    end
    if ismissing(value)
        text = 'string[1 1]:missing';
    else
        text = ['string' encode_text(char(value))];
    end
elseif ischar(value)
    text = encode_text(value);
else
    error('ERW:checkpointRuntimeTable', ...
        'Unsupported runtime-parameter type %s.',class(value));
end
end

function text = encode_text(value)
units = uint16(value);
body = sprintf('%04X',units(:).');
text = sprintf('char%s:%s',mat2str(size(value)),body);
end
