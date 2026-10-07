function test_crop_reserve_state_boundary
% Crop reserve write-back must not overwrite the preceding retained state.

code_dir = fileparts(fileparts(mfilename('fullpath')));
sources = {fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_10mins.m')), ...
    fileread(fullfile(code_dir,'MAIN_FRAME_BG_ERW_2D.m'))};

for k = 1:numel(sources)
    assert(isempty(regexp(sources{k}, ...
        'Nreserve_H\(j-1:j[^\)]*\).*CropHeightType','once')), ...
        'CropHeightType must write current reserves only, not j-1:j.');
end

fprintf('PASS: crop reserve write-back preserves the preceding state.\n');
end
