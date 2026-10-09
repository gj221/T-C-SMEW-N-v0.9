function tag = Bounded_File_Tag(tag,knobs_from,knobs_to,prefix,extension)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Keep a result or checkpoint file name within the 255 characters a file
%% name may have on Linux and Windows.
%%
%% Run_One_Subfield names a file prefix + tag + extension, the tag
%% stringing together the run, scenario, calcite case and grid, then one
%% token per knob passed (alpha, the five exchange constants, KNmaxd, the
%% biogeochemistry overrides, the dissolution factor, the saturation
%% exponent, Kbot), then the span (days, start date). tag(knobs_from+1:
%% knobs_to) are the knob tokens. A name that fits is returned as it is,
%% so every name that could be written before stays the same. A name that
%% would not fit keeps its run, scenario, case, grid and span tokens and
%% has its knob tokens replaced by '__knobs' and the first twelve hex
%% digits of their SHA-256 digest; run_provenance carries every knob.
%%
%% The uncertainty ensemble on a field grid needs this: the upper set's
%% knobs and the grid tag together give 258-272 characters, and the move
%% that puts a finished result or a checkpoint in place would fail.
%%
%%   tag = Bounded_File_Tag(tag,knobs_from,knobs_to,prefix,extension)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_name = 255;
if ~(ischar(tag) && ischar(prefix) && ischar(extension))
    error('Bounded_File_Tag:BadInput','tag, prefix and extension must be character vectors.');
end
if ~(isscalar(knobs_from) && isscalar(knobs_to) && knobs_from == fix(knobs_from) && ...
        knobs_to == fix(knobs_to) && knobs_from >= 0 && knobs_from <= knobs_to && ...
        knobs_to <= numel(tag))
    error('Bounded_File_Tag:BadInput', ...
        'The knob tokens must be a range of the tag (%d to %d of %d characters).', ...
        knobs_from,knobs_to,numel(tag));
end
if numel(prefix)+numel(tag)+numel(extension) <= max_name
    return
end
knobs = tag(knobs_from+1:knobs_to);
if ~isempty(knobs)
    digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
    digest.update(typecast(uint8(knobs),'int8'));
    raw = typecast(digest.digest(),'uint8');
    hex = lower(reshape(dec2hex(raw,2).',1,[]));
    tag = [tag(1:knobs_from) '__knobs' hex(1:12) tag(knobs_to+1:end)];
end
if numel(prefix)+numel(tag)+numel(extension) > max_name
    error('Bounded_File_Tag:NameTooLong', ...
        ['%s%s%s has %d characters even with its knobs as a digest; ' ...
        'a file name may have %d.'],prefix,tag,extension, ...
        numel(prefix)+numel(tag)+numel(extension),max_name);
end
end
