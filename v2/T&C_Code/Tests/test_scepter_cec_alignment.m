function test_scepter_cec_alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCEPTER-source exchange equation and the fixed lumped Kantola CEC.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== test_scepter_cec_alignment ===\n');
n_pass = 0;
here = fileparts(mfilename('fullpath'));
code_dir = fileparts(here);
v2_dir = fileparts(code_dir);
repo_root = fileparts(v2_dir);
campaign_dir = fullfile(repo_root,'runtime', ...
    'Calibration_Sweep_1D');
addpath(code_dir,campaign_dir);

%% The lumped TC Appelo treatment keeps soil CEC fixed and varies selectivity
%% through alpha and f_H.  SCEPTER itself fixes each solid's mass-specific CEC.
assert_error(@() ERW_Resolve_CEC_pH_Sensitivity( ...
    struct('dCEC_dpH',2.9,'CEC_ref_pH',6.1),'appelo'), ...
    'ERW:scepterFixedCEC');
check(ERW_Resolve_CEC_pH_Sensitivity(struct(),'appelo') == 0, ...
    'TC Appelo mode retains fixed lumped soil CEC');

%% SCEPTER 1.0.2 source-equation reference, converted from mol/L to umol/L.
H = 1;
Na = 500; K = 200; Ca = 1000; Mg = 500;
Al_w = 1e-6;
[f_Ca,f_Mg,f_K,f_Na,f_Al,f_H] = Exchange_Fractions_GT( ...
    H,Ca,Mg,K,Na,Al_w,ones(1,5),0,0,0,0,1,'appelo');
actual = [f_Na,f_K,f_Ca,f_Mg,f_Al,f_H];
expected = [0.0025001760,0.0125901406,0.5345658002, ...
    0.1291133575,2.1232530e-9,0.3212305234802847];
check(max(abs(actual-expected)) < 5e-10, ...
    'Appelo fractions reproduce the SCEPTER reference case');
check(abs(sum(actual)-1) < 1e-11, ...
    'SCEPTER reference exchange fractions close');

%% Direct Energy Farm control-soil CEC is the common pretreatment capacity.
workbook = fullfile(code_dir,'Inputs','ERW_MOD_PARAM_parameters.xlsx');
T = readtable(workbook,'Sheet','ERW_MOD');
cec = T.CEC_tot_in(1:5);
check(all(abs(cec-13) < eps(13)), ...
    'all Kantola treatments use the measured 13 cmolc kg-1 capacity');

%% Alpha changes hold selectivity at the seed's f_H = 0.20, so every alpha
%% starts from the measured 80% base saturation.
cases = Sweep_Cases;
reference_alpha = 1.3;
reference_logKH = [4.71627429367 3.61627429367 8.10254858734 ...
    8.41854858734 12.918822881];
charge = [1 1 2 2 3];
for alpha = unique([cases.AppeloAlpha])
    these = cases([cases.AppeloAlpha] == alpha);
    expected_logKH = reference_logKH + ...
        charge*(alpha-reference_alpha)*0.20;
    check(all(arrayfun(@(x) max(abs(x.AppeloLogKH-expected_logKH)) < 1e-12,these)), ...
        sprintf('alpha %.3g holds the seed selectivity',alpha));
end

fprintf('  ---------------------------------------\n');
fprintf('  %d passed, 0 failed\n',n_pass);

    function check(ok,label)
        assert(ok,'test_scepter_cec_alignment:Failed','%s',label);
        n_pass = n_pass+1;
        fprintf('  PASS  %s\n',label);
    end

    function assert_error(f,expected_id)
        try
            f();
        catch ME
            assert(strcmp(ME.identifier,expected_id), ...
                'Expected %s, got %s: %s',expected_id,ME.identifier,ME.message);
            fprintf('  PASS  rejects incompatible input with %s\n',expected_id);
            n_pass = n_pass+1;
            return
        end
        error('Expected error %s was not raised.',expected_id);
    end
end
