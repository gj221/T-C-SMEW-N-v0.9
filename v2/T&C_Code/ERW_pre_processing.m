function [excel_var_names, ERW_preprocess_output] = ...
    ERW_pre_processing(vars,ERW_preprocess_data)

% excel_input.fX_conc_convention:
% "conc_to_f_CEC" (concentrations -> CEC fractions)
% "total_to_f_CEC_and_conc" (total element pools -> f_CEC + conc)
% "f_CEC_to_conc" (CEC fractions -> concentrations)
% "f_CEC_inventory" (fractions define conserved initial inventories)

% excel_input.rock_filepath: datafile for rock distribution

for ind_vars =1:length(vars)
    if isvarname(vars{ind_vars})
        eval([vars{ind_vars},'=ERW_preprocess_data.',vars{ind_vars},';'])
    end
end
if ~exist('Al_CEC_scale','var')
    Al_CEC_scale = conv_Al;
end

%% fX_conc_convention

if excel_input.fX_conc_convention == "f_CEC_inventory"
    f_all = [f_Ca(1), f_Mg(1), f_K(1), f_Na(1), f_Al(1), f_H(1)];
    if any(~isfinite(f_all)) || any(f_all < 0) || ...
            f_Ca(1) <= 0 || abs(sum(f_all) - 1) > 1e-6
        error('ERW:InvalidCECFractions', ...
            'f_CEC_inventory requires nonnegative fractions summing to one, with positive f_Ca.');
    end

elseif excel_input.fX_conc_convention == "conc_to_f_CEC"

    % f_CEC from porewater concentrations via Gaines-Thomas exchange
    if ~isfinite(Na(1))
        Na(1) = 0;
    end
    if ~isfinite(Al_w(1))
        Al_w(1) = Initialize_ERW_Equilibrium(H(1), K1, K2, K3, K4, conv_mol, K_AlOH3, Al_CEC_scale);
    end
    Al = (H(1)^4/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4)) * Al_w(1);

    K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H];
    [f_Ca(1), f_Mg(1), f_K(1), f_Na(1), f_Al(1), f_H(1)] = ...
        Exchange_Fractions_GT(H(1), Ca(1), Mg(1), K(1), Na(1), Al_w(1), ...
        K_CEC, K1, K2, K3, K4, Al_CEC_scale);


elseif excel_input.fX_conc_convention == "total_to_f_CEC_and_conc"

    f_base_CEC = 0.31;
    f_acid_CEC = 1 - f_base_CEC;
    f_acid = f_acid_CEC;

    Ca_tot(1) = 214*1e-6*rho_bulk*Zrw/MM_Ca;
    Mg_tot(1) = 34.7*1e-6*rho_bulk*Zrw/MM_Mg;
    K_tot(1) = 40*1e-6*rho_bulk*Zrw/MM_K;
    Na_tot(1) = 9*1e-6*rho_bulk*Zrw/MM_Na;

    f_Ca0 = 0.8*(1 - f_acid);
    f_Mg0 = 0.1*(1 - f_acid);
    f_K0 = 0.05*(1 - f_acid);
    f_Na0 = 0.05*(1 - f_acid);
    f_H0 = f_acid*0.3;
    f_Al0 = f_acid*0.7;

    Si(1) = 0;

    equations = @(p)[p(1)*n*Zrw*s(1)*1000+(p(12)/3)*CEC_tot*Al_CEC_scale-p(3);
        p(2)-(H(1)^4/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4))*p(1);
        p(4)*n*Zrw*s(1)*1000+p(8)/2*CEC_tot-Mg_tot(1);
        p(5)*n*Zrw*s(1)*1000+p(11)/2*CEC_tot-Ca_tot(1);
        p(6)*n*Zrw*s(1)*1000+p(9)*CEC_tot-Na_tot(1);
        p(7)*n*Zrw*s(1)*1000+p(10)*CEC_tot-K_tot(1);
        p(12) - (p(2)/Al_CEC_scale)*(p(11)^3/(K_Ca_Al*p(5)^3))^(1/2);
        p(13) - H(1)*(p(11)/(K_Ca_H*p(5)))^(1/2);
        p(13) + p(12) - f_acid;
        p(8) - p(4)*(p(11)/(K_Ca_Mg*p(5)));
        p(9) - p(6)*(p(11)/(K_Ca_Na*p(5)))^(1/2);
        p(10) - p(7)*(p(11)/(K_Ca_K*p(5)))^(1/2);
        1-(p(11)+p(12)+p(8)+p(9)+p(10)+p(13))];

    Mg0 = (Mg_tot(1)-f_Mg0/2*CEC_tot(1))/(n*Zrw*s(1)*1000);
    Ca0 = (Ca_tot(1)-f_Ca0/2*CEC_tot(1))/(n*Zrw*s(1)*1000);
    K0 = (K_tot(1)-f_K0*CEC_tot(1))/(n*Zrw*s(1)*1000);
    Na0 = (Na_tot(1)-f_Na0*CEC_tot(1))/(n*Zrw*s(1)*1000);
    Al0 = (f_Al0*Al_CEC_scale)*((K_Ca_Al*Ca0^3)/f_Ca0^3)^(1/2);
    Al_w0 = Al0/(H(1)^4/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4));
    Al_tot0 = (f_Al0/3)*CEC_tot(1)*Al_CEC_scale+Al_w0*n*Zrw*s(1)*1000;

    x0 = [Al_w0, Al0, Al_tot0, Mg0, Ca0, Na0, K0, f_Mg0, f_Na0, f_K0, f_Ca0, f_Al0, f_H0];

    options = optimoptions(fsolve_options, 'FunctionTolerance', 1e-14);
    sol = fsolve(equations, x0, options);

    Al_w(1) = sol(1);
    Al(1) = sol(2);
    Al_tot(1) = sol(3);
    Mg(1) = sol(4);
    Ca(1) = sol(5);
    Na(1) = sol(6);
    K(1) = sol(7);
    f_Mg(1) = sol(8);
    f_Na(1) = sol(9);
    f_K(1) = sol(10);
    f_Ca(1) = sol(11);
    f_Al(1) = sol(12);
    f_H(1) = sol(13);

    if exist('opt_TC_biochem','var') && opt_TC_biochem == 1
        K(1) = K_in;
        Ca(1) = Ca_in;
        Mg(1) = Mg_in;
        Si(1) = Si_in;
    end


elseif excel_input.fX_conc_convention == "f_CEC_to_conc"
    f_all = [f_Ca(1), f_Mg(1), f_K(1), f_Na(1), f_Al(1), f_H(1)];
    if any(~isfinite(f_all)) || any(f_all < 0) || ...
            f_Ca(1) <= 0 || f_H(1) <= 0 || abs(sum(f_all) - 1) > 1e-6
        error('ERW:InvalidCECFractions', ...
            'f_CEC_to_conc requires nonnegative fractions summing to one, with positive f_Ca and f_H.');
    end
    Ca(1) = (f_Ca(1)/K_Ca_H)*(H(1)/f_H(1))^2;
    Mg(1) = (Ca(1)/f_Ca(1))*K_Ca_Mg*f_Mg(1);
    K(1) = ((Ca(1)/f_Ca(1))*K_Ca_K)^(1/2)*f_K(1);
    Si(1) = 0;
    Na(1) = ((Ca(1)/f_Ca(1))*K_Ca_Na)^(1/2)*f_Na(1);
    Al(1) = ((Ca(1)/f_Ca(1))^3*K_Ca_Al)^(1/2)*f_Al(1)*Al_CEC_scale;
    Al_w(1) = Al(1)/(H(1)^4/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4));


elseif excel_input.fX_conc_convention == "Kantola"

    kantola_scheme = [];
    if exist('ERW_TC_const_input','var') && isstruct(ERW_TC_const_input) && ...
            isscalar(ERW_TC_const_input) && ...
            isfield(ERW_TC_const_input,'bulk_pH_params') && ...
            isstruct(ERW_TC_const_input.bulk_pH_params) && ...
            isscalar(ERW_TC_const_input.bulk_pH_params) && ...
            isfield(ERW_TC_const_input.bulk_pH_params,'scheme')
        kantola_scheme = ERW_TC_const_input.bulk_pH_params.scheme;
    end
    kantola_scheme_is_text = ...
        (ischar(kantola_scheme) && isrow(kantola_scheme) && ~isempty(kantola_scheme)) || ...
        (isstring(kantola_scheme) && isscalar(kantola_scheme) && ...
         ~ismissing(kantola_scheme));
    if ~kantola_scheme_is_text || ...
            ~strcmpi(string(kantola_scheme), "appelo")
        error('ERW:KantolaRequiresAppelo', ...
            ['The Kantola exchange seed requires bulk_pH_params.scheme = ' ...
             '''appelo''; TC_GT is retained only for legacy-compatible ' ...
             'non-Kantola paths.']);
    end

    % case study: Kantola.
    %
    % The legacy T&C-SMEW seed prescribes an exchange composition
    % (f_CEC_sat = 0.8, f_H = 0.15) and back-derives the pore-water
    % concentrations from the Gaines-Thomas inverses, then fits
    % K_Ca_H/K_Ca_Mg/K_Ca_K to reproduce that same state.  Under GT that
    % is self-consistent: the fitted constants absorb whatever the
    % assumed composition implies, so the absolute Na and Al values never
    % have to be physical.
    %
    % This model runs Appelo/SCEPTER selectivity with FIXED literature
    % constants (Resolve_Appelo_Constants), which removes that freedom.
    % The assumed composition and the derived concentrations become
    % physical claims, and two of them are false: handing sodium every
    % site the prescribed Ca/Mg/K do not claim gives f_Na = 0.76, an ESP of
    % 76% (a sodic soil, not a Mollisol), and the GT inverses evaluated at
    % f_Ca ~ 0.003 return 0.5 M Na and 89 mM Al, the latter 2e5 times the
    % Al(OH)3 solubility at pH 6.
    %
    % The prescribed 80% base-saturation target is retained for
    % this soil, so it is kept.  What changes is where the unmeasured
    % remainder goes.  The spin-up's exchangeable pools (P(:,61/69/53))
    % come back essentially equal to the solution pools and account for
    % only ~4% of the CEC, so the remainder is real and has to be
    % assigned; a Mollisol carries it as Ca and Mg, not sodium.
    % Proportions below are an explicitly assumed exchangeable-base split for a
    % near-neutral Mollisol and are exposed as constants so a site with
    % its own exchange data can override them.
    %% Pretreatment base saturation is approximately 80% in Fig. 1F.
    f_CEC_sat = 0.80;
    %% The 20% exchangeable acidity is carried entirely by H.
    %%
    %% It was previously split f_H = 0.15, f_Al = 0.05, but that 5% Al is not
    %% supportable here. This configuration runs with Al_solubility_pin = 1,
    %% so dissolved Al is held at amorphous Al(OH)3 saturation, and at the
    %% site pH of about 6 that supports f_Al near 4e-4 -- some 125x short of
    %% 0.05. The Appelo activity term 10^(alpha*3*f_H) cannot close that gap
    %% either: the whole range from the SCEPTER default alpha of 3.4 down to
    %% the 1.3 used here spans only 8.8x. The seed was therefore prescribing
    %% an exchange composition the chemistry strips on the first solve,
    %% moving the state off the measured 80% anchor before the run starts.
    %%
    %% Carrying the acidity as H is also the better soil-chemical claim at
    %% this pH: exchangeable Al only dominates below roughly pH 5.5, and a
    %% near-neutral Mollisol at pH 6.2 holds its exchangeable acidity as H
    %% and hydroxy-Al rather than as exchangeable Al3+. Al remains free to
    %% populate the complex from solution if the chemistry supports it; it
    %% is simply no longer asserted at initialisation.
    f_H(1) = 1 - f_CEC_sat;
    f_Al(1) = 0;

    % Pore-water Na is prescribed or deposition-derived, so it is set first and the
    % exchange fraction follows from the Appelo law rather than the other
    % way round.  It must be strictly positive: the log-space exchange
    % residual log(f_Na) - (log K_Na_H + log[Na] + log(f_H/[H]) + ...)
    % has no root at [Na] = 0, where both sides clamp onto different
    % floors and leave a fixed ~13 log-unit residual that pins the scaled
    % error at 1 and stalls every other equation in the block.  Use the
    % deposition-supported steady state when the site supplies a rain-Na
    % flux, otherwise an explicit fallback assumption of 200 umol/L
    % (~4.6 mg/L), which is mid-range for a leached non-sodic soil.
    if ~(isfinite(Na(1)) && Na(1) > 0)
        Na(1) = 200;
    end
    % The seed must be built on the same exchanger the run applies. Calling
    % Resolve_Appelo_Constants([]) takes the SCEPTER default alpha of 3.4,
    % but a run supplies its own through bulk_pH_params.appelo (1.3 for the
    % Kantola campaign). The activity term is 10^(alpha*Z*f_H), so at
    % f_H = 0.15 the two differ by 10^((3.4-1.3)*0.15) = 2.1 for sodium and
    % by the cube of that, 8.8, for aluminium. A seed built at 3.4 therefore
    % prescribes an exchange composition the run cannot hold, and the first
    % solve strips it: f_Al falls from the prescribed 0.05 to ~4e-4, leaving
    % exchangeable acidity carried entirely by H, which titrates far more
    % readily than Al and overstates the pH response to feedstock.
    appelo_params_kant = [];
    if exist('ERW_TC_const_input','var') && isstruct(ERW_TC_const_input) && ...
            isfield(ERW_TC_const_input,'bulk_pH_params') && ...
            isstruct(ERW_TC_const_input.bulk_pH_params) && ...
            isfield(ERW_TC_const_input.bulk_pH_params,'appelo')
        appelo_params_kant = ERW_TC_const_input.bulk_pH_params.appelo;
    else
        warning('ERW:KantolaSeedAlphaDefault', ...
            ['Kantola seed could not read the run Appelo parameters and is ' ...
             'using the SCEPTER defaults. The seed exchanger will not match ' ...
             'the run if the run overrides alpha or logKH.']);
    end
    [K_Na_H_kant, ~, ~, ~, ~, alpha_kant] = ...
        Resolve_Appelo_Constants(appelo_params_kant);
    f_Na(1) = K_Na_H_kant*Na(1)*(f_H(1)/H(1))*10^(alpha_kant*f_H(1));
    if ~(isfinite(f_Na(1)) && f_Na(1) >= 0 && f_Na(1) <= f_CEC_sat)
        error('ERW:InfeasibleKantolaNaExchange', ...
            ['Pore-water Na implies exchange fraction %.6g outside the ' ...
             'available base-saturation interval [0, %.6g].'], ...
            f_Na(1),f_CEC_sat);
    end

    % Ca/Mg/K carry the remaining prescribed base saturation in the
    % assumed exchangeable split for a near-neutral Mollisol; exposed as
    % constants so a site with its own exchange data can override them.
    kantola_base_split = [0.79, 0.17, 0.04];   % Assumed Ca Mg K partition
    kantola_base_split = kantola_base_split/sum(kantola_base_split);
    f_base_kant = f_CEC_sat - f_Na(1);
    f_Ca(1) = f_base_kant*kantola_base_split(1);
    f_Mg(1) = f_base_kant*kantola_base_split(2);
    f_K(1)  = f_base_kant*kantola_base_split(3);

    % Dissolved Al is set by Al(OH)3 solubility, which is what controls
    % it in any soil near pH 6 (Lindsay 1979 p.51); SCEPTER holds the
    % same phase as its equilibrium solid 'amal'.  Never by the GT
    % inverse, which is unbounded as f_Ca falls.
    Al_w(1) = Initialize_ERW_Equilibrium(H(1), K1, K2, K3, K4, ...
        conv_mol, K_AlOH3, Al_CEC_scale);
    Al_den_kant = H(1)^4 + H(1)^3*K1 + H(1)^2*K1*K2 + ...
        H(1)*K1*K2*K3 + K1*K2*K3*K4;
    Al(1) = (H(1)^4/Al_den_kant)*Al_w(1);

    % estimates of K_CEC
    K_Ca_H = (f_Ca(1)/Ca(1))*(H(1)/f_H(1))^2;
    K_Ca_Mg = (f_Ca(1)/Ca(1))*(Mg(1)/f_Mg(1));
    K_Ca_K = (f_Ca(1)/Ca(1))*(K(1)/f_K(1))^2;


else
    error('ERW:UnknownCECConvention', 'Unknown fX_conc_convention: %s', ...
        char(string(excel_input.fX_conc_convention)));
end



K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H];



%% rock_dist_data
rock_filepath_txt = strtrim(string(excel_input.rock_filepath));
rock_f_check = [];
if exist('rock_f_in','var')
    rock_f_check = rock_f_in;
end
rock_dist_active = false;
if isnumeric(rock_f_check)
    rock_dist_active = any(isfinite(rock_f_check(:))) && any(rock_f_check(:) > 0);
end
if rock_dist_active && ~ismissing(rock_filepath_txt) && strlength(rock_filepath_txt) > 0 && ...
        lower(rock_filepath_txt) ~= "nan"
    % particle size distribution - diameter and class weight
    % data from json file
    rock_data = fileread(char(rock_filepath_txt));
    rock_data = jsondecode(rock_data);

    rock_type = excel_input.rock_type;

    feedstock_dataset = rock_data.datasetColl(find(strcmp({rock_data.datasetColl.name}, rock_type), 1));

    x1 = zeros(1, length(feedstock_dataset.data));
    y1 = zeros(1, length(feedstock_dataset.data));
    for a = 1:length(feedstock_dataset.data)
        x1(a) = feedstock_dataset.data(a).value(1);
        y1(a) = feedstock_dataset.data(a).value(2);
    end

    rock_dist_mode = lower(strtrim(string(excel_input.rock_dist_data)));
    if rock_dist_mode == "cdf_um" || rock_dist_mode == "kantola"
        d_in = x1 * 1e-6; % [m]
        CDF_in = y1;
        d_in(1) = 1e-6;
        psd_in = ERW_Particle_Mass_Fractions(d_in,CDF_in);
    elseif rock_dist_mode == "mass_fraction_um"
        d_in = x1 * 1e-6;
        psd_in = y1 / sum(y1);
    else
        error('ERW:UnknownRockDistribution', ...
            'rock_dist_data must be cdf_um or mass_fraction_um when rock_filepath is set.');
    end
end

%% save values
ERW_preprocess_output = struct;

excel_var_names = {'f_Ca';'f_Mg';'f_K';'f_Na';'f_Al';'f_H';
    'K_Ca_Mg';'K_Ca_K';'K_Ca_Na';'K_Ca_Al';'K_Ca_H';
    'Ca';'Mg';'K';'Na';'An';'Al';'Si';'H';'Al_w';
    'd_in';'psd_in';'pH';'Al_CEC_scale'};

for var_index = 1:size(excel_var_names,1)
    vn = cell2mat(excel_var_names(var_index));
    if exist(vn, 'var')
        ERW_preprocess_output.(vn) = eval(vn);
    end
end

end
