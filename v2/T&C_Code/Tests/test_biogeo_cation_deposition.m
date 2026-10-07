function test_biogeo_cation_deposition
% Each element must retain its own measured deposition flux.

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(code_dir);
S = load(fullfile(code_dir,'Inputs','All_deposition_data.mat'));
location = load(fullfile(code_dir,'Inputs','Data_Kantola.mat'),'Lat','Lon');
fields = {'DepCa','DepMg','DepK'};
indices = [1 4 3];
for historical = [0 1]
    B = Biogeochemistry_IO(300,location.Lat,location.Lon,0,historical);
    for k = 1:numel(indices)
        j = indices(k);
        annual = griddata(S.LonCC{j},S.LatCC{j},S.CC{j}, ...
            location.Lon,location.Lat,'nearest');
        actual = B.(fields{k});
        assert(isscalar(actual) && isfinite(actual) && actual >= 0);
        assert(abs(actual*365*10-annual) < 1e-12, ...
            '%s must recover its own annual kg/ha deposition.',fields{k});
    end
end
% Kantola runs without atmospheric base-cation deposition in every period,
% so each element is assigned exactly once, to zero, after B_IO is built.
param = fileread(fullfile(code_dir,'Parameters','MOD_PARAM_Kantola_Maize_Soy.m'));
built = strfind(param,'[B_IO]=Biogeochemistry_IO(');
assert(isscalar(built),'Kantola must build B_IO exactly once.');
for k = 1:numel(fields)
    assigned = regexp(param,['B_IO\.' fields{k} '\s*=[^;\n]*;'],'match');
    assert(isequal(assigned,{sprintf('B_IO.%s = 0;',fields{k})}) && ...
        strfind(param,assigned{1}) > built, ...
        'Kantola must zero %s once, after building B_IO.',fields{k});
end
fprintf('PASS: Ca, Mg and K keep their own fluxes; Kantola zeroes all three.\n');
end
