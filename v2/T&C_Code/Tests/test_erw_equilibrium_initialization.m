function test_erw_equilibrium_initialization
% Initial dissolved Al must invert the selected equilibrium exactly.

test_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(test_dir));

H = 1;
K = [10 8 0.2 0.631];
scale = 1000;

Ca = 1000;
f_Ca = 0.5;
f_Al = 0.1;
K_Ca_Al = 2;
Al3_exchange = f_Al*scale*sqrt(K_Ca_Al*Ca^3/f_Ca^3);
den = H^4+H^3*K(1)+H^2*K(1)*K(2)+ ...
    H*K(1)*K(2)*K(3)+prod(K);
expected_exchange = Al3_exchange*den/H^4;
actual_exchange = Initialize_ERW_Equilibrium(H,Ca,f_Ca,f_Al, ...
    K_Ca_Al,K(1),K(2),K(3),K(4),scale);
assert(abs(actual_exchange/expected_exchange-1) < 1e-12, ...
    'Ca-Al exchange inversion does not recover total dissolved Al.');

conv_mol = 1e6;
logK = 10.8;
Al3_saturation = 10^logK*(H/conv_mol)^3*conv_mol*scale;
expected_saturation = Al3_saturation*den/H^4;
actual_saturation = Initialize_ERW_Equilibrium(H,K(1),K(2),K(3),K(4), ...
    conv_mol,logK,scale);
assert(abs(actual_saturation/expected_saturation-1) < 1e-12, ...
    'Al(OH)3 saturation inversion does not recover total dissolved Al.');
assert(isnan(Initialize_ERW_Equilibrium(0,K(1),K(2),K(3),K(4), ...
    conv_mol,logK,scale)), ...
    'An invalid hydrogen concentration must not return an equilibrium state.');

fprintf('PASS: ERW dissolved-Al initialization inverts both equilibria.\n');
end
