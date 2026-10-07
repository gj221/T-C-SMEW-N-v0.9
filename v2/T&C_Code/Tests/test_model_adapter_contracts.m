function test_model_adapter_contracts
% The 2-D adapters must preserve parameter fields without reinterpretation.

test_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(test_dir));

required = {'Osat','Ohy','L','Pe','O33','alpVG','nVG','lVG', ...
    'Ks_mac','Omac','alpVGM','nVGM','lVGM','s_SVG','bVG'};
soil = struct('Ks_Zs',[1 2],'Phy',-10);
for k = 1:numel(required)
    soil.(required{k}) = k;
end
mapped = Soil_Param_To_Routing(soil,2,1.5);
assert(mapped.SPAR == 2 && isequal(mapped.Ks,soil.Ks_Zs) && ...
    mapped.Phy1 == soil.Phy && mapped.aR == 1.5, ...
    'Routing adapter changed the supplied soil parameters.');

bad_soil = rmfield(soil,'Osat');
assert_throws(@() Soil_Param_To_Routing(bad_soil,2,1), ...
    'Soil_Param_To_Routing:MissingField');
assert_throws(@() Load_Cell_Parameters('missing_parameter_file.m','',struct()), ...
    'Load_Cell_Parameters:MissingParamFile');

legacy_frame = fileread(fullfile(fileparts(test_dir),'MAIN_FRAME_BG.m'));
assert(contains(legacy_frame,'Nuptake_H(j,:),~,~,Puptake_H(j,:)') && ...
       contains(legacy_frame,'Nuptake_L(j,:),~,~,Puptake_L(j,:)'), ...
    'The standalone GT frame shifts BIOGEO_UNIT outputs after Nuptake.');

fprintf('PASS: 2-D parameter adapters preserve their explicit contracts.\n');
end

function assert_throws(fun,identifier)
threw = false;
try
    fun();
catch ME
    threw = strcmp(ME.identifier,identifier);
end
assert(threw,'Expected error %s.',identifier);
end
