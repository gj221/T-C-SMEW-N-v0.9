function test_plant_tissue_dilution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Legacy T&C plant nutrients. A reserve below zero is tissue diluted below
%% its nominal content (rXc = 1 + Xreserve/TXns). Growth may draw the
%% reserve down until the non-structural tissue sits at 0.65 of nominal,
%% so availability is Xreserve + 0.35*TXns, clipped at zero. Each nutrient
%% in turn cuts NPP to what its availability supports and the rejected
%% carbon leaves as extra respiration. Translocation stops only when some
%% availability is exhausted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(code_dir);
warning_state = warning('off','all');
restore = onCleanup(@() warning(warning_state)); %#ok<NASGU>

St = Veg_Stoichiometric_Parameter(22,'Maize');
ratios = [St.Nl St.Nr St.Ns St.Nf; St.Phol St.Phor St.Phos St.Phof; ...
    St.Kpotl St.Kpotr St.Kpots St.Kpotf; St.Call St.Calr St.Cals St.Calf; ...
    St.Magl St.Magr St.Mags St.Magf; St.Sill St.Silr St.Sils St.Silf];

%% Availability and relative concentration, Nutrients_Available.
B = [10 0 120 80 5 0 5 0];
tissue = B([1 3 4 5])*(1./ratios');
for fraction = [0.1 -0.2 -0.5]
    reserve = fraction*tissue;
    [r,avail] = nutrients(B,St,reserve);
    assert(fraction > 0 || max(abs(r-(1+fraction))) < 1e-12, ...
        'A reserve deficit must read as tissue diluted below nominal.');
    assert(max(abs(avail-max(reserve+0.35*tissue,0))) < 1e-12, ...
        'Availability must be the reserve plus 0.35 of the non-structural tissue.');
end
for bad = {[-1.2*tissue(1) zeros(1,5)], [NaN zeros(1,5)]}
    rejected = false;
    try
        nutrients(B,St,bad{1});
    catch ME
        rejected = strcmp(ME.identifier,'Nutrients_Available:InvalidReserve');
    end
    assert(rejected, ...
        'A nonfinite reserve or a deficit beyond the tissue must be rejected.');
end
[r_empty,avail_empty] = nutrients(zeros(1,8),St,zeros(1,6));
assert(isequal(r_empty,ones(1,6)) && isequal(avail_empty,zeros(1,6)), ...
    'A plant with no biomass and no reserve has nothing to dilute or supply.');

%% Growth limitation, VEGETATION_DYNAMIC. With every loss term off, a day's
%% N draw is the nominal N of the day's tissue growth.
ample = 1e3*ones(1,6);
grow = struct('B',[50 20 30 40 5 0 0 0],'An',12,'PHE_S',3,'r',0.03,'Rdark',1);
dB_ref = vegetation(grow,St,ample);
npp = sum(dB_ref(1:5));
draw = @(dB) dB([1 3 4 5])'*(1./ratios') + dB(2)./ratios(:,3)';
demand = draw(dB_ref);
assert(npp > 0 && dB_ref(8) == 0 && all(demand > 0), ...
    'The reference day must grow without nutrient limitation.');

limited = ample;
limited(1) = 0.4*demand(1);
dB = vegetation(grow,St,limited);
used = draw(dB);
assert(abs(used(1)-limited(1)) < 1e-12 && ...
    max(abs(dB(1:5)-0.4*dB_ref(1:5))) < 1e-12, ...
    'A limiting nutrient must cut growth to the tissue its availability supports.');
assert(abs(sum(dB(1:5))+dB(8)-npp) < 1e-12, ...
    'Carbon rejected by nutrient limitation must leave as extra respiration.');

limited(2) = 0.3*demand(2);
limited(1) = 0.5*demand(1);
dB = vegetation(grow,St,limited);
assert(max(abs(dB(1:5)-0.3*dB_ref(1:5))) < 1e-12 && ...
    abs(sum(dB(1:5))+dB(8)-npp) < 1e-12, ...
    'The scarcest nutrient must set growth (Liebig).');

%% Translocation from the carbohydrate reserve at leaf-out is not scaled
%% with growth; it stops only when some availability is exhausted. NPP is
%% then zero and each pool pays its own maintenance, as on any day of
%% negative NPP, so a plant without maintenance holds still.
leafout = grow;
leafout.PHE_S = 2;
dB_ample = vegetation(leafout,St,ample);
translocation = -dB_ample(4);
assert(translocation > 0, 'The leaf-out reference must translocate.');
scarce = ample;
scarce(1) = 1e-9;
dB = vegetation(leafout,St,scarce);
assert(abs(dB(4)+translocation) < 1e-12 && dB(8) > 0.99*npp, ...
    'A scarce but positive availability must stop growth, not translocation.');
scarce(1) = 0;
dB = vegetation(leafout,St,scarce);
assert(abs(sum(dB(1:5))) < 1e-12 && abs(dB(8)-npp) < 1e-12, ...
    'An exhausted availability must stop net growth.');
leafout.r = 0;
leafout.Rdark = 0;
dB = vegetation(leafout,St,scarce);
assert(all(abs(dB(1:5)) < 1e-12) && dB(8) > 0, ...
    'An exhausted availability must stop translocation.');

%% A day of negative NPP pays leaf dark respiration from the reserve. That
%% transfer replaces the leaf-out translocation, roots included, so carbon
%% closes whatever the phenological state.
dark = grow;
dark.An = -0.5;
dB_growth_phase = vegetation(dark,St,ample);
dark.PHE_S = 2;
dB_leafout = vegetation(dark,St,ample);
assert(sum(dB_growth_phase(1:5)) < 0 && ...
    max(abs(dB_leafout-dB_growth_phase)) < 1e-12, ...
    'The dark-respiration transfer must not add root carbon.');

fprintf('PASS: plant nutrients follow the legacy tissue-dilution law.\n');
end

function [r,avail] = nutrients(B,St,reserve)
out = cell(1,16);
[out{:}] = Nutrients_Available(B,1,St,reserve(1),reserve(2),reserve(3), ...
    reserve(4),reserve(5),reserve(6),1);
r = [out{3:8}];
avail = [out{11:16}];
end

function dB = vegetation(state,St,avail)
dB = VEGETATION_DYNAMIC(0,state.B',20,18,state.An,state.Rdark,1,1,1,1,18,1,Inf, ...
    0.02,0,St,state.r,1,0.25,5,2,0,0,-5,0,0,90,state.PHE_S,0,0,10,1, ...
    0,0,1,0,0.03,0.1,0,0,0,avail(1),avail(2),avail(3),avail(4), ...
    avail(5),avail(6),0.3,1,0,0,0);
end
