function test_leached_salt_closure
%%% The 0-30 cm box exports what its drainage water carries. Since the ERW
%%% module took over the cation partition, T&C's Ca and Mg pools P(60) and
%%% P(68) are the dissolved pools, so their leaching and passive-uptake
%%% coefficient is 1, as it already is for NO3 and K (Porporato et al.,
%%% 2003, define a as the dissolved fraction of a lumped pool). At 0.4 the
%%% leachate carried 2.5 times less Ca+Mg charge than its anions, the
%%% released alkalinity stayed in the box and the exchanger saturated by
%%% 2019 (budget of package m, 1 October 2026).
%%%
%%% With the cations leaving at their concentration, the legacy steady
%%% closure that returns the leached residual anion as a bare anion would
%%% become an acid input of about 0.65 eq m-2 yr-1 in the control, so the
%%% leached part comes back as the salt it left as through the I_Ca and
%%% I_Mg boundary slots.

test_dir = fileparts(mfilename('fullpath'));
code_dir = fullfile(test_dir,'..');
addpath(code_dir,fullfile(code_dir,'Inputs'),fullfile(code_dir,'Parameters'));
BP = Biogeochemistry_Parameter(1);

%% Leaching: a dissolved pool leaves at concentration times drainage.
B = zeros(1,80);
B(60) = 1.2; B(68) = 0.6; B(52) = 0.3;      % g m-2 in solution
Lk = 20; V = 80;                             % mm d-1 drained of mm stored
[~,~,~,LEAK_K,~,~,~,LEAK_Ca,LEAK_Mg] = Biogeo_Leakage(B,Lk,V,BP);
assert(abs(LEAK_Ca-B(60)*Lk/V) <= 1e-12 && abs(LEAK_Mg-B(68)*Lk/V) <= 1e-12, ...
    'Dissolved Ca and Mg must leach at concentration times drainage (a = 1).');
assert(abs(LEAK_K/B(52)-LEAK_Ca/B(60)) <= 1e-12, ...
    'Every dissolved pool must see the same drainage fraction.');
assert(BP.aCa == 1 && BP.aMg == 1 && BP.aK == 1 && BP.aNO3 == 1, ...
    'The dissolved-fraction coefficients of Ca, Mg, K and NO3 must all be 1.');

%% Passive uptake: the transpiration stream carries the same concentration.
%% A vanishing root mass leaves the active term at zero without tripping the
%% Broot == 0 shortcut.
T = 10; VT = 50; Ccrown = 1;
[~,~,~,~,Ca_up,Mg_up] = Biogeo_uptake(B,1e-9,15,T,VT,Ccrown,0.5,BP);
assert(abs(Ca_up-B(60)*T/VT) <= 1e-9 && abs(Mg_up-B(68)*T/VT) <= 1e-9, ...
    'Passive Ca and Mg uptake must be concentration times transpiration (a = 1).');

%% Salt re-supply: the charge of the leached residual anion that Na did not
%% carry returns on Ca and Mg in their pore-water proportions.
f_Lw = 0.03; An_tot = 1.42e6; Na_tot = 0.111e6; Ca = 9000; Mg = 4000;
[I_Ca,I_Mg] = ERW_Leached_Salt_Resupply(f_Lw,An_tot,Na_tot,Ca,Mg);
charge = f_Lw*(An_tot-Na_tot);
assert(abs(2*I_Ca+2*I_Mg-charge) <= 1e-9*charge, ...
    'The returned salt must carry exactly the leached anion charge Na did not.');
assert(abs(I_Ca/I_Mg-Ca/Mg) <= 1e-12*Ca/Mg, ...
    'The salt cations must return in their pore-water proportions.');
[I_Ca,I_Mg] = ERW_Leached_Salt_Resupply(f_Lw,Na_tot,An_tot,Ca,Mg);
assert(I_Ca == 0 && I_Mg == 0, ...
    'When Na covers the leached anion charge nothing is owed to Ca or Mg.');
[I_Ca,I_Mg] = ERW_Leached_Salt_Resupply(f_Lw,An_tot,Na_tot,0,0);
assert(abs(I_Ca-I_Mg) <= 1e-12*I_Ca && abs(4*I_Ca-charge) <= 1e-9*charge, ...
    'With no Ca or Mg in solution the charge splits equally.');
[I_Ca,I_Mg] = ERW_Leached_Salt_Resupply(0,An_tot,Na_tot,Ca,Mg);
assert(I_Ca == 0 && I_Mg == 0,'No drainage, no salt to return.');
failed = false;
try
    ERW_Leached_Salt_Resupply(-f_Lw,An_tot,Na_tot,Ca,Mg);
catch ME
    failed = strcmp(ME.identifier,'ERW:invalidParam');
end
assert(failed,'Negative or nonfinite inputs must be refused.');

%% Both drivers return the leached salt the same way, so the 1-D frame and
%% the 2-D driver stay bit for bit equal.
for f = {'ERW_Daily_Step_2D.m','MAIN_FRAME_BG_ERW_10mins.m'}
    assert(contains(fileread(fullfile(code_dir,f{1})),'ERW_Leached_Salt_Resupply('), ...
        '%s must return the leached salt through ERW_Leached_Salt_Resupply.',f{1});
end

fprintf('PASS: dissolved Ca and Mg leach at concentration and the leached salt returns neutral.\n');
end
