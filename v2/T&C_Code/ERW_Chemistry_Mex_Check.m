function ERW_Chemistry_Mex_Check()
%%% Refuse TC_ERW_KERNEL=mex unless the compiled kernels reproduce the
%%% MATLAB kernels bit for bit.
%%%
%%% Runs once per MATLAB process (or worker), the first time a step asks for
%%% 'mex': each kernel is evaluated compiled and in MATLAB on reference
%%% steps and the outputs must be isequaln.  A missing MEX, one built from
%%% other kernel sources, or one whose math library rounds differently on
%%% this platform stops the run here instead of changing its results.  On a
%%% thread-based worker MATLAB itself refuses the MEX call with its own
%%% message; such pools run TC_ERW_KERNEL=fast.

persistent passed
if ~isempty(passed) && passed
    return
end
if exist('ERW_Appelo_Seed_Kernel_mex','file') ~= 3 || ...
        exist('ERW_Charge_Root_Kernel_mex','file') ~= 3
    error('ERW:mexKernelMissing', ...
        ['TC_ERW_KERNEL=mex needs ERW_Appelo_Seed_Kernel_mex and ' ...
         'ERW_Charge_Root_Kernel_mex on the path; build them with ' ...
         'Build_ERW_Chemistry_Mex on this platform.']);
end
for k = 1:3
    args = reference_seed(k);
    [state_m,status_m] = ERW_Appelo_Seed_Kernel(args{:});
    [state_x,status_x] = ERW_Appelo_Seed_Kernel_mex(args{:});
    if ~(isequaln(state_m,state_x) && isequal(status_m,status_x))
        error('ERW:mexKernelMismatch', ...
            ['ERW_Appelo_Seed_Kernel_mex differs from ERW_Appelo_Seed_Kernel ' ...
             'on reference step %d; rebuild it with Build_ERW_Chemistry_Mex ' ...
             'or run TC_ERW_KERNEL=fast.'],k);
    end
    args = reference_scan(k);
    out_m = cell(1,7);
    out_x = cell(1,7);
    [out_m{:}] = ERW_Charge_Root_Kernel(args{:});
    [out_x{:}] = ERW_Charge_Root_Kernel_mex(args{:});
    if ~isequaln(out_m,out_x)
        error('ERW:mexKernelMismatch', ...
            ['ERW_Charge_Root_Kernel_mex differs from ERW_Charge_Root_Kernel ' ...
             'on reference step %d; rebuild it with Build_ERW_Chemistry_Mex ' ...
             'or run TC_ERW_KERNEL=fast.'],k);
    end
end
passed = true;
end

function args = reference_seed(k)
%%% Kantola-like steps: Appelo constants at alpha 1.3, 3.4 and 6, pore pH
%%% near 5.5, 7 and 8.3, finite Al(OH)3 inventory under the pin.
logKH = [4.71627429367 3.61627429367 8.10254858734 8.41854858734 12.918822881];
alpha = [1.3 3.4 6];
alpha = alpha(k);
K = 10.^(-(logKH + [1 1 2 2 3]*(alpha-1.3)*0.20));
pH = [5.5 7 8.3];
H = 10^(-pH(k))*1e6;
conv_mol = 1e6;
pH_scan = [1e-8:0.25:13.75,14-1e-8];
H_scan = 10.^(-pH_scan)*conv_mol;
soil_water_l = 107.2;
s = 0.699;
bulk_volume_l = soil_water_l/s;
CEC_tot = 5.05e7;
carbonate = [0.371 3.8e-5 0.0072 0.0303 conv_mol];
cation_tot = [1.51e7 1.21e7 2.28e6 3.49e4].*[1 1 1 1]*(0.6+0.4*k);
Al_tot = [2.03e6 1.0e3 0];
finite_al = [Al_tot(k) 1000 10.8 10 7.94 0.1995 0.631];
appelo = [K(3) K(4) K(5) K(2) K(1) alpha];
gas = [9.2 18.2 1/144];
Alk_tot = sum(cation_tot.*[2 2 1 1]) - 6.1e6 + 400*H^-1*soil_water_l;
IC_tot = 1.2e3*k;
args = {H_scan,Alk_tot,IC_tot,cation_tot,soil_water_l,bulk_volume_l,s, ...
    carbonate,CEC_tot,appelo,gas,10.8,1000,finite_al,true};
end

function args = reference_scan(k)
%%% Pore-water charge scans at the same three steps.
conv_mol = 1e6;
H_scan = 10.^(-linspace(1e-8,14-1e-8,281))*conv_mol;
k1 = 0.371; k2 = 3.8e-5; k_w = 0.0072;
CO2_w0 = [40 120 600];
Alk_tot = [2.1e6 5.5e6 1.2e7];
R_alk0 = [2.0e6 5.0e6 1.1e7];
soil_water_l = 107.2;
Al_w0 = [6.2e3 3.1 1e-6];
args = {H_scan,k1,k2,k_w,CO2_w0(k),Alk_tot(k),R_alk0(k),soil_water_l, ...
    Al_w0(k),1000,10,7.94,0.1995,0.631};
end
