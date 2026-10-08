function test_erw_chemistry_kernels(n_seed_cases,n_chain_steps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TC_ERW_KERNEL=fast (and =mex, where built) reproduces the original ERW
%% chemistry bit for bit.
%%
%% ERW_Chemistry_Kernel_Mode selects the original code ('matlab', the
%% default), the restructured MATLAB kernels ('fast') or their MATLAB Coder
%% builds ('mex').  Every comparison is bit for bit (bit_identical: same
%% class, size and IEEE bits, struct fields in the same order, function
%% handles by text and captured values):
%%   1. ERW_Appelo_Seed against ERW_Appelo_Mass_Conserving_State on steps
%%      built to put the charge root anywhere in pore pH 3.5 to 9.5, over
%%      temperature, ionic strength, alpha 1.3/3.4/6, CEC, soil water,
%%      untreated and treated cation pools and the Al(OH)3 inventory
%%      (none, inventory-limited, saturated), plus steps whose root lies
%%      outside the scan, where both must raise the same error;
%%   2. ERW_Charge_Root_Kernel against the original scan and bisection;
%%   3. all seven ERW_function outputs on the step-22053 fixture perturbed
%%      across pore pH, the Al pin on and off, temperature, saturation and
%%      alpha, and along chains of steps that feed each state to the next.
%% The optional arguments set the number of seed cases (default 600) and
%% chain steps (default 12); larger values give the full report.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(n_seed_cases), n_seed_cases = 600; end
if nargin < 2 || isempty(n_chain_steps), n_chain_steps = 12; end
fprintf('\n=== test_erw_chemistry_kernels ===\n');
here = fileparts(mfilename('fullpath'));
code_dir = fileparts(here);
addpath(code_dir,fullfile(code_dir,'Inputs'),fullfile(code_dir,'Parameters'));

old_mode = getenv('TC_ERW_KERNEL');
restore_mode = onCleanup(@() setenv('TC_ERW_KERNEL',old_mode));
warn_state = warning('off','ERW:NoConverge');
restore_warn = onCleanup(@() warning(warn_state));

modes = {'fast'};
if exist('ERW_Appelo_Seed_Kernel_mex','file') == 3 && ...
        exist('ERW_Charge_Root_Kernel_mex','file') == 3
    modes{end+1} = 'mex';
else
    fprintf('  SKIP  mex: ERW_Appelo_Seed_Kernel_mex / ERW_Charge_Root_Kernel_mex not built here\n');
end
setenv('TC_ERW_KERNEL','');
assert(strcmp(ERW_Chemistry_Kernel_Mode(struct()),'matlab'), ...
    'test_erw_chemistry_kernels:Default', ...
    'With TC_ERW_KERNEL unset the original code must run.');
assert(strcmp(ERW_Chemistry_Kernel_Mode(struct('erw_kernel','fast')),'fast'), ...
    'test_erw_chemistry_kernels:Override', ...
    'bulk_pH_params.erw_kernel must select the kernel.');

%% 1. The reduced Appelo solve.
rng(20261008,'twister');
n_err = 0;
pH_root = zeros(1,n_seed_cases);
for c = 1:n_seed_cases
    args = seed_case(c,n_seed_cases);
    [ref,ref_err] = call_capture(@ERW_Appelo_Mass_Conserving_State,args{1:14});
    for m = 1:numel(modes)
        [got,got_err] = call_capture(@ERW_Appelo_Seed,modes{m},args{:});
        assert(bit_identical(got,ref) && isequal(got_err,ref_err), ...
            'test_erw_chemistry_kernels:Seed', ...
            'Seed case %d (%s) differs from ERW_Appelo_Mass_Conserving_State.', ...
            c,modes{m});
    end
    if isempty(ref_err)
        pH_root(c) = -log10(ref.H/1e6);
    else
        n_err = n_err + 1;
        pH_root(c) = NaN;
    end
end
valid = isfinite(pH_root);
fprintf('  PASS  seed: %d cases identical (%s); roots at pore pH %.2f to %.2f, %d errors raised identically\n', ...
    n_seed_cases,strjoin(modes,', '),min(pH_root(valid)),max(pH_root(valid)),n_err);

%% 2. The pore-water charge scan and bisection.
n_scan_cases = max(100,round(n_seed_cases/3));
for c = 1:n_scan_cases
    args = scan_case(c);
    ref = cell(1,8);
    [ref{:}] = call_capture_n(8,@original_charge_scan,args{:});
    for m = 1:numel(modes)
        got = cell(1,8);
        if strcmp(modes{m},'mex')
            [got{:}] = call_capture_n(8,@kernel_charge_scan_mex,args{:});
        else
            [got{:}] = call_capture_n(8,@kernel_charge_scan,args{:});
        end
        assert(bit_identical(got,ref), 'test_erw_chemistry_kernels:Scan', ...
            'Scan case %d (%s) differs from the original scan.',c,modes{m});
    end
end
fprintf('  PASS  charge scan: %d cases identical (%s)\n',n_scan_cases,strjoin(modes,', '));

%% 3. ERW_function on the fixture step.
[args0,names] = fixture_args(code_dir);
idx = @(name) find(strcmp(names,name),1);
variants = fixture_variants(args0,names);
for v = 1:numel(variants)
    a = variants{v};
    setenv('TC_ERW_KERNEL','matlab');
    [ref,ref_err] = call_capture_n(7,@ERW_function,a{:});
    for m = 1:numel(modes)
        setenv('TC_ERW_KERNEL',modes{m});
        [got,got_err] = call_capture_n(7,@ERW_function,a{:});
        assert(isequal(got_err,ref_err), 'test_erw_chemistry_kernels:FunctionError', ...
            'ERW_function variant %d (%s) raised a different error.',v,modes{m});
        [same,where] = bit_identical(got,ref);
        assert(same, 'test_erw_chemistry_kernels:Function', ...
            'ERW_function variant %d (%s) differs at %s.',v,modes{m},where);
    end
end
setenv('TC_ERW_KERNEL','');
fprintf('  PASS  ERW_function: %d fixture variants, all 7 outputs identical (%s)\n', ...
    numel(variants),strjoin(modes,', '));

chains = {args0, set_arg(args0,idx('An_tot'),args0{idx('An_tot')}*1.3), ...
    set_arg(set_arg(args0,idx('An_tot'),args0{idx('An_tot')}*2.5), ...
        idx('temp_soil'),6)};
for c = 1:numel(chains)
    ref_traj = run_chain(chains{c},names,'matlab',n_chain_steps);
    for m = 1:numel(modes)
        got_traj = run_chain(chains{c},names,modes{m},n_chain_steps);
        assert(bit_identical(got_traj,ref_traj), 'test_erw_chemistry_kernels:Chain', ...
            'Chain %d (%s) leaves the original trajectory.',c,modes{m});
    end
end
setenv('TC_ERW_KERNEL','');
fprintf('  PASS  %d chains of %d steps identical (%s)\n',numel(chains), ...
    n_chain_steps,strjoin(modes,', '));
end

%% ------------------------------------------------------------------------
function args = seed_case(c,n_cases)
%%% One reduced-solve call as ERW_function builds it, with the charge root
%%% placed at a chosen pore pH (test_appelo_mass_conserving_state's
%%% construction) and a jitter so it sits at no constructed value.
conv_mol = 1e6;
T_ref = 298.15;
T_K = 273.15 + 2 + 28*rand;
pk1 = -(-356.309-0.0609*T_K+21834.37/T_K+126.8339*log10(T_K)-1684915/(T_K)^2);
pk2 = -(-107.887-0.032528*(T_K)+5151.79/(T_K)+38.92561*log10(T_K)-563713.9/(T_K)^2);
pk_w = -(-283.971+13323/(T_K)-0.0507*(T_K)+102.24447*log10(T_K)-1119669/(T_K)^2);
k1 = 10^(-pk1)*conv_mol;
k2 = 10^(-pk2)*conv_mol;
k_w = 10^(-pk_w)*conv_mol^2;
k_H = 0.83*exp(2400*(1/T_K-1/T_ref));
I = 10^(-3+2*rand);
[k_cond,K_AlOH3] = ERW_Conditional_Constants( ...
    [k1,k2,k_w,10,7.94328,0.199526,0.630957],10.8,I);
k1 = k_cond(1); k2 = k_cond(2); k_w = k_cond(3);
K1 = k_cond(4); K2 = k_cond(5); K3 = k_cond(6); K4 = k_cond(7);
alphas = [1.3 3.4 6];
logKH = [4.71627429367 3.61627429367 8.10254858734 8.41854858734 12.918822881];
alpha = alphas(mod(c,3)+1);
[K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,alpha] = Resolve_Appelo_Constants( ...
    struct('logKH',logKH + [1 1 2 2 3]*(alpha-1.3)*0.20,'alpha',alpha));
appelo = [K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,alpha];
n_por = 0.45 + 0.1*rand;
Zrw = 0.3;
s = 0.12 + 0.88*rand;
soil_water_l = n_por*Zrw*s*1000;
bulk_volume_l = n_por*Zrw*1000;
CEC_tot = 5.05e7*10^(-0.3+0.6*rand);
treated = rand < 0.5;
cation_tot = [1.51e7 1.21e7 2.28e6 3.49e4].*10.^(-0.5+rand(1,4));
if treated
    cation_tot(1:2) = cation_tot(1:2).*(1.5+2*rand(1,2));
end
Al_choices = [0 1e3 3e4 2.03e6 2e7];
Al_tot = Al_choices(randi(numel(Al_choices)));
Al_CEC_scale = 1000;
finite_al = struct('Al_tot',Al_tot,'Al_CEC_scale',Al_CEC_scale, ...
    'K_AlOH3',K_AlOH3,'K1',K1,'K2',K2,'K3',K3,'K4',K4);
al3_pin_fn = @(Hp) 10^(K_AlOH3-12) * max(real(Hp), 1e-12)^3 * Al_CEC_scale;
al3_fn = @(Hp) al3_pin_fn(Hp)/Al_CEC_scale;
D = 1.3824*(1-min(s,0.999))^(10/3)*n_por^(4/3);
gas = [D/0.15*1000, 18.18*(0.8+0.4*rand), 1/144];
carbonate = [k1,k2,k_w,k_H,conv_mol];
%% Pore pH spread over 3.5-9.5; a few steps far outside the scan.
if c > n_cases - max(3,round(0.02*n_cases))
    pH_ref = [-0.6 14.8 15.5];
    pH_ref = pH_ref(mod(c,3)+1);
else
    pH_ref = 3.5 + 6*(c-1)/max(1,n_cases-1);
end
H_ref = 10^(-pH_ref)*conv_mol;
CO2_ref = 10^(1+2*rand);
[f_ref,~,Al3_ref] = ERW_Appelo_State_At_H(H_ref,cation_tot,al3_fn(H_ref), ...
    CEC_tot,soil_water_l,alpha,K_Na_H,K_K_H,K_Ca_H2,K_Mg_H2,K_Al_H3,finite_al);
if any(~isfinite(f_ref))
    f_ref = zeros(1,6);
    Al3_ref = 0;
end
R_alk_ref = sum(f_ref(1:4))*CEC_tot;
Alk_ref = k1*CO2_ref/H_ref + 2*k1*k2*CO2_ref/H_ref^2 - H_ref + k_w/H_ref;
Al_net = Al3_ref*(3+2*K1/H_ref+K1*K2/H_ref^2-K1*K2*K3*K4/H_ref^4);
Alk_tot = (R_alk_ref + (Alk_ref-Al_net)*soil_water_l)*(1+1e-4*randn);
IC_final = CO2_ref*bulk_volume_l*(s*(1+k1/H_ref+k1*k2/H_ref^2)+(1-s)/k_H);
Fs_ref = gas(1)*(CO2_ref/k_H-gas(2));
IC_tot = max(IC_final+Fs_ref*gas(3),0);
H_anchor = H_ref*10^(0.5*randn);
cec_fn = @(Hx) CEC_tot;
args = {Alk_tot,IC_tot,cation_tot,soil_water_l,bulk_volume_l,s, ...
    carbonate,CEC_tot,appelo,gas,al3_fn,H_anchor,finite_al,cec_fn, ...
    K_AlOH3,Al_CEC_scale};
end

function args = scan_case(c)
%%% ERW_function's pore-water charge scan inputs over the same range.
conv_mol = 1e6;
T_K = 273.15 + 2 + 28*rand;
pk1 = -(-356.309-0.0609*T_K+21834.37/T_K+126.8339*log10(T_K)-1684915/(T_K)^2);
pk2 = -(-107.887-0.032528*(T_K)+5151.79/(T_K)+38.92561*log10(T_K)-563713.9/(T_K)^2);
pk_w = -(-283.971+13323/(T_K)-0.0507*(T_K)+102.24447*log10(T_K)-1119669/(T_K)^2);
[k_cond] = ERW_Conditional_Constants([10^(-pk1)*conv_mol,10^(-pk2)*conv_mol, ...
    10^(-pk_w)*conv_mol^2,10,7.94328,0.199526,0.630957],10.8,10^(-3+2*rand));
k1 = k_cond(1); k2 = k_cond(2); k_w = k_cond(3);
K = k_cond(4:7);
soil_water_l = 20 + 140*rand;
CO2_w0 = 10^(0.5+2.5*rand);
pH0 = 3.5 + 6*rand;
H0 = 10^(-pH0)*conv_mol;
Al_w0 = 10^(-1.2)*H0^3*1000*10^(-2+4*rand);
R_alk0 = 5e7*rand;
Alk_pore = k1*CO2_w0/H0 + 2*k1*k2*CO2_w0/H0^2 - H0 + k_w/H0;
Alk_tot = R_alk0 + Alk_pore*soil_water_l*(1+0.2*randn);
if mod(c,25) == 0
    Alk_tot = -1e12;                       %% no root inside the scan
end
H_scan = 10.^(-linspace(1e-8,14-1e-8,281))*conv_mol;
args = {H_scan,k1,k2,k_w,CO2_w0,Alk_tot,R_alk0,soil_water_l,Al_w0,1000, ...
    K(1),K(2),K(3),K(4)};
end

function [R_scan,lo_H,hi_H,lo_r,iter,mid_H,R_mid,H_robust] = original_charge_scan( ...
    H_scan,k1,k2,k_w,CO2_w0,Alk_tot,R_alk0,soil_water_l,Al_w0,Al_CEC_scale,K1,K2,K3,K4)
%%% ERW_function's scan and bisection before the kernels, verbatim.
al_charge_at = @(H,Al_w) (Al_w/Al_CEC_scale)* ...
    (3*H^4+2*K1*H^3+K1*K2*H^2-K1*K2*K3*K4)/ ...
    (H^4+K1*H^3+K1*K2*H^2+K1*K2*K3*H+K1*K2*K3*K4);
charge_res_fn = @(XH) (k1*CO2_w0/XH + 2*k1*k2*CO2_w0/(XH^2) - XH + k_w/XH) - ...
    ((Alk_tot-R_alk0)/soil_water_l+al_charge_at(XH,Al_w0));
R_scan = nan(size(H_scan));
for ii = 1:numel(H_scan)
    R_scan(ii) = charge_res_fn(H_scan(ii));
end
valid_scan = isfinite(R_scan);
sign_change = find(valid_scan(1:end-1) & valid_scan(2:end) & ...
    (R_scan(1:end-1) == 0 | R_scan(2:end) == 0 | sign(R_scan(1:end-1)) ~= sign(R_scan(2:end))));
sign_change = sign_change([true,diff(sign_change) > 1]);
if numel(sign_change) == 1
    lo_H = H_scan(sign_change);
    hi_H = H_scan(sign_change+1);
    lo_r = charge_res_fn(lo_H);
    for iter = 1:60
        mid_H = 0.5*(lo_H + hi_H);
        R_mid = charge_res_fn(mid_H);
        if ~isfinite(R_mid)
            hi_H = mid_H;
        elseif abs(R_mid) < max(1e-8, 1e-10) || abs(hi_H-lo_H)/max(mid_H,eps) < 1e-12
            break
        elseif sign(R_mid) == sign(lo_r)
            lo_H = mid_H;
            lo_r = R_mid;
        else
            hi_H = mid_H;
        end
    end
    H_robust = 0.5*(lo_H + hi_H);
else
    error('ERW:nonUniquePorePH','found %d',numel(sign_change));
end
end

function varargout = kernel_charge_scan(varargin)
varargout = cell(1,8);
[varargout{:}] = kernel_charge_scan_with(@ERW_Charge_Root_Kernel,varargin{:});
end

function varargout = kernel_charge_scan_mex(varargin)
varargout = cell(1,8);
[varargout{:}] = kernel_charge_scan_with(@ERW_Charge_Root_Kernel_mex,varargin{:});
end

function [R_scan,lo_H,hi_H,lo_r,iter,mid_H,R_mid,H_robust] = ...
    kernel_charge_scan_with(kernel,H_scan,varargin)
%%% The kernel and ERW_function's own handling of its scan.
[R_scan,lo_H,hi_H,lo_r,iter,mid_H,R_mid] = kernel(H_scan,varargin{:});
valid_scan = isfinite(R_scan);
sign_change = find(valid_scan(1:end-1) & valid_scan(2:end) & ...
    (R_scan(1:end-1) == 0 | R_scan(2:end) == 0 | sign(R_scan(1:end-1)) ~= sign(R_scan(2:end))));
sign_change = sign_change([true,diff(sign_change) > 1]);
if numel(sign_change) == 1
    H_robust = 0.5*(lo_H + hi_H);
else
    error('ERW:nonUniquePorePH','found %d',numel(sign_change));
end
end

function [args,names] = fixture_args(code_dir)
%%% The ERW_function call refused at step 22053 (test_erw_appelo_exchange_floor).
erw_source = fileread(fullfile(code_dir,'ERW_function.m'));
F = load(fullfile(code_dir,'Tests','Fixtures','erw_appelo_al_floor_step22053.mat'),'fixture');
F = F.fixture;
header_start = regexp(erw_source,'ERW_function\(','end','once') + 1;
header_end = regexp(erw_source,'\)\r?\n\r?\n%%%% INPUT','start','once') - 1;
names = strtrim(strsplit(regexprep( ...
    erw_source(header_start:header_end),'\.\.\.',''),','));
args = cell(1,numel(names));
for k = 1:numel(names)
    if ~strcmp(names{k},'~')
        args{k} = F.inputs.(names{k});
    end
end
end

function variants = fixture_variants(a0,names)
%%% Pore pH from about 4 to 11, the Al pin on and off, cold and warm soil,
%%% dry and wet soil, alpha 1.3/3.4/6 and a treated pool.
idx = @(name) find(strcmp(names,name),1);
iAn = idx('An_tot'); iT = idx('temp_soil'); is = idx('s'); isp = idx('s_pre');
iCa = idx('Ca_tot'); iMg = idx('Mg_tot'); ibp = idx('bulk_pH_params');
variants = {};
for m = [0.98 1 1.15 1.5 2 3 4.5 6 8]
    variants{end+1} = set_arg(a0,iAn,a0{iAn}*m); %#ok<AGROW>
end
for m = [1.15 3 6]
    a = set_arg(a0,iAn,a0{iAn}*m);
    bp = a{ibp}; bp.Al_solubility_pin = false;
    variants{end+1} = set_arg(a,ibp,bp); %#ok<AGROW>
end
for T = [2 12 30]
    variants{end+1} = set_arg(set_arg(a0,iAn,a0{iAn}*2),iT,T); %#ok<AGROW>
end
for sv = [0.35 0.55 0.95]
    a = set_arg(set_arg(a0,iAn,a0{iAn}*1.5),is,sv);
    variants{end+1} = set_arg(a,isp,sv*1.01); %#ok<AGROW>
end
for alpha = [1.3 3.4 6]
    a = set_arg(a0,iAn,a0{iAn}*2);
    bp = a{ibp}; bp.appelo.alpha = alpha;
    variants{end+1} = set_arg(a,ibp,bp); %#ok<AGROW>
end
for m = [1.6 2.4]
    a = set_arg(a0,iCa,a0{iCa}*m);
    a = set_arg(a,iMg,a0{iMg}*m);
    variants{end+1} = set_arg(a,iAn,a0{iAn}*m*1.4); %#ok<AGROW>
end
end

function traj = run_chain(a,names,mode,n_steps)
%%% Successive steps, each fed the state the previous one returned, as
%%% ERW_Daily_Step_2D feeds row dt_m-1 to row dt_m.
setenv('TC_ERW_KERNEL',mode);
idx = @(name) find(strcmp(names,name),1);
carry = {'Na','An','CO2_air','Ca_tot','Mg_tot','K_tot', ...
    'Na_tot','An_tot','Al_tot','Si_tot','Al','AlOH4','R_alk','H','IC_tot', ...
    'Fs','DIC','CaCO3','MgCO3','f_CEC','f_Ca','f_Mg','f_K','f_Na','f_Al', ...
    'f_H','PO4_tot','NO3_tot','NH4_tot'};
carry = carry(ismember(carry,names));
traj = cell(n_steps,1);
for step = 1:n_steps
    a{idx('time_step')} = a{idx('time_step')} + 1;
    a{idx('s_pre')} = a{idx('s')};
    try
        [~,~,B,W] = ERW_function(a{:});
    catch ME
        traj{step} = struct('identifier',ME.identifier,'message',ME.message);
        traj = traj(1:step);
        break
    end
    traj{step} = B;
    for k = 1:numel(carry)
        a{idx(carry{k})} = B.(carry{k});
    end
    a{idx('W_CaCO3')} = W.W_CaCO3;
    a{idx('W_MgCO3')} = W.W_MgCO3;
end
setenv('TC_ERW_KERNEL','');
end

function a = set_arg(a,k,value)
a{k} = value;
end

function [out,err] = call_capture(fn,varargin)
out = [];
err = [];
try
    out = fn(varargin{:});
catch ME
    err = struct('identifier',ME.identifier,'message',ME.message);
end
end

function varargout = call_capture_n(n,fn,varargin)
%%% n outputs of fn, or the error it raised as the last output.
varargout = cell(1,n);
outs = cell(1,n);
try
    [outs{:}] = fn(varargin{:});
    err = [];
catch ME
    err = struct('identifier',ME.identifier,'message',ME.message);
end
if nargout == 2
    varargout{1} = outs;
    varargout{2} = err;
else
    varargout = outs;
    if ~isempty(err)
        varargout{1} = err;
    end
end
end

function [same,where] = bit_identical(a,b,where,depth)
%%% a and b hold the same value bit for bit: class, size, the IEEE bits of
%%% every float (signed zeros and NaN payloads included), struct fields in
%%% the same order, cells element by element, function handles by their
%%% text and, for anonymous ones, by the values they captured (two calls
%%% create distinct handles, which isequaln never equates).
if nargin < 3, where = 'value'; end
if nargin < 4, depth = 0; end
same = strcmp(class(a),class(b)) && isequal(size(a),size(b));
if ~same
    where = [where ' (class or size)'];
    return
end
if isfloat(a)
    same = isreal(a) == isreal(b) && ...
        isequal(typecast(real(a(:)),float_bits(a)),typecast(real(b(:)),float_bits(b))) && ...
        isequal(typecast(imag(a(:)),float_bits(a)),typecast(imag(b(:)),float_bits(b)));
elseif isstruct(a)
    fa = fieldnames(a);
    same = isequal(fa,fieldnames(b));
    if ~same
        where = [where ' (field names)'];
        return
    end
    for e = 1:numel(a)
        for f = 1:numel(fa)
            [same,where_f] = bit_identical(a(e).(fa{f}),b(e).(fa{f}), ...
                sprintf('%s(%d).%s',where,e,fa{f}),depth);
            if ~same, where = where_f; return, end
        end
    end
elseif iscell(a)
    for e = 1:numel(a)
        [same,where_e] = bit_identical(a{e},b{e},sprintf('%s{%d}',where,e),depth);
        if ~same, where = where_e; return, end
    end
elseif isa(a,'function_handle')
    same = strcmp(func2str(a),func2str(b));
    if same && depth < 6
        fa = functions(a);
        fb = functions(b);
        if isfield(fa,'workspace') && isfield(fb,'workspace')
            [same,where_w] = bit_identical(fa.workspace,fb.workspace, ...
                [where ' workspace'],depth+1);
            if ~same, where = where_w; return, end
        end
    end
else
    same = isequaln(a,b);
end
if ~same
    where = [where ' (value)'];
end
end

function c = float_bits(x)
if isa(x,'single')
    c = 'uint32';
else
    c = 'uint64';
end
end
