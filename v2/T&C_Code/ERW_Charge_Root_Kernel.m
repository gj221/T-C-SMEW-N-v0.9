function [R_scan,lo_H,hi_H,lo_r,iter,mid_H,R_mid] = ERW_Charge_Root_Kernel( ...
    H_scan,k1,k2,k_w,CO2_w0,Alk_tot,R_alk0,soil_water_l,Al_w0, ...
    Al_CEC_scale,K1,K2,K3,K4) %#codegen
%%% The pore-water charge scan and bisection of ERW_function, with the same
%%% arithmetic in the same order and without its anonymous-function calls.
%%%
%%% ERW_function evaluates charge_res_fn, which calls al_charge_at, at the
%%% 281 scan points and then bisects the one bracketed root.  This kernel
%%% evaluates the two handles' expressions verbatim on the same operands,
%%% the scan elementwise over all points at once (MATLAB's elementwise
%%% +, -, *, / and .^ round exactly as its scalar forms), so R_scan and
%%% every bisection variable are bit-identical to the original's
%%% (test_erw_chemistry_kernels).  H_scan is the caller's
%%% 10.^(-linspace(1e-8,14-1e-8,281))*conv_mol; the caller still forms
%%% valid_scan and sign_change from R_scan with its own code, so a scan with
%%% no or several roots stops exactly where it stopped before.
%%%
%%% The bisection outputs are the loop's final lo_H, hi_H, lo_r, iter,
%%% mid_H and R_mid when the scan holds exactly one run of sign changes,
%%% and NaN otherwise.

n_scan = numel(H_scan);
XH = reshape(H_scan,1,n_scan);
H = XH;
Al_w = Al_w0;
al_charge = (Al_w/Al_CEC_scale)* ...
    (3*H.^4+2*K1*H.^3+K1*K2*H.^2-K1*K2*K3*K4)./ ...
    (H.^4+K1*H.^3+K1*K2*H.^2+K1*K2*K3*H+K1*K2*K3*K4);
R_scan = (k1*CO2_w0./XH + 2*k1*k2*CO2_w0./(XH.^2) - XH + k_w./XH) - ...
    ((Alk_tot-R_alk0)/soil_water_l+al_charge);
lo_H = NaN; hi_H = NaN; lo_r = NaN; iter = NaN; mid_H = NaN; R_mid = NaN;

%% sign_change = find(...); sign_change([true,diff(sign_change) > 1])
n_runs = 0;
first = 0;
previous = -1;
for k = 1:n_scan-1
    r0 = R_scan(k);
    r1 = R_scan(k+1);
    if isfinite(r0) && isfinite(r1) && ...
            (r0 == 0 || r1 == 0 || sign(r0) ~= sign(r1))
        if n_runs == 0 || k - previous > 1
            n_runs = n_runs + 1;
            if n_runs == 1
                first = k;
            end
        end
        previous = k;
    end
end
if n_runs ~= 1
    return
end

lo_H = H_scan(first);
hi_H = H_scan(first+1);
lo_r = R_scan(first);
for iter = 1:60
    mid_H = 0.5*(lo_H + hi_H);
    R_mid = charge_res(mid_H,k1,k2,k_w,CO2_w0,Alk_tot,R_alk0, ...
        soil_water_l,Al_w0,Al_CEC_scale,K1,K2,K3,K4);
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
end

function R = charge_res(XH,k1,k2,k_w,CO2_w0,Alk_tot,R_alk0,soil_water_l, ...
    Al_w0,Al_CEC_scale,K1,K2,K3,K4)
%%% charge_res_fn(XH) with al_charge_at(XH,Al_w0) written out.
H = XH;
Al_w = Al_w0;
al_charge = (Al_w/Al_CEC_scale)* ...
    (3*H^4+2*K1*H^3+K1*K2*H^2-K1*K2*K3*K4)/ ...
    (H^4+K1*H^3+K1*K2*H^2+K1*K2*K3*H+K1*K2*K3*K4);
R = (k1*CO2_w0/XH + 2*k1*k2*CO2_w0/(XH^2) - XH + k_w/XH) - ...
    ((Alk_tot-R_alk0)/soil_water_l+al_charge);
end
