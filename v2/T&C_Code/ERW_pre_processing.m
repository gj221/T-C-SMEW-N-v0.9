function [excel_var_names, ERW_preprocess_output] = ...
    ERW_pre_processing(vars,ERW_preprocess_data)

% excel_input.fX_conc_convention:
% "Amann" (for case study "Amann")
% "total_to_f_CEC_and_conc" (for case study "te_pas")
% "f_CEC_to_conc" (for case study "Dietzen")
% "conc_to_f_CEC" (for case study "Kelland")
% "Kelland"


% excel_input.rock_dist_data:
% "Amann"
% "te_Pas"
% "Kelland"

% excel_input.rock_filepath: datafile for rock distribution
% "Amann_psd.json"

for ind_vars =1:length(vars)
    eval([vars{ind_vars},'=ERW_preprocess_data.',vars{ind_vars},';'])
end

%% fX_conc_convention

if excel_input.fX_conc_convention == "Amann"

    % case study: Amann
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % INPUT CEC fractions
    % f_Ca = 0.87-(0.005);
    % f_Mg = 0.05;
    % f_K = 0.05;
    % f_Na = 0.03;
    % f_Al = 0.002;
    % f_H = 0.003;

    % estimates of concentrations
    Ca(1) = (Mg(1)/f_Mg(1))/K_Ca_Mg*f_Ca(1);
    K(1) = ((Ca(1)/f_Ca(1))*K_Ca_K)^(1/2)*f_K(1);

    if excel_input.opt_TC_biochem == 1
        Ca(1) = Ca_in; % [micromol/l]
        K(1) = K_in;
    end

    Na(1) = ((Ca(1)/f_Ca(1))*K_Ca_Na)^(1/2)*f_Na(1);
    Al(1) = ((Ca(1)/f_Ca(1))^3*K_Ca_Al)^(1/2)*f_Al(1)*conv_Al;
    Al_w(1) = Al(1)/(H(1)^4/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4));

    % estimates of K_CEC
    K_Ca_H = (f_Ca(1)/Ca(1))*(H(1)/f_H(1))^2;

elseif excel_input.fX_conc_convention == "total_to_f_CEC_and_conc"

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_base_CEC = 0.31;
    f_acid_CEC = 1 - f_base_CEC;
    f_acid = f_acid_CEC;

    Ca_tot(1) = 214*1e-6*rho_bulk*Zrw/MM_Ca; % mg/kg_soil converted to mol-conv/m2
    Mg_tot(1) = 34.7*1e-6*rho_bulk*Zrw/MM_Mg;
    K_tot(1) = 40*1e-6*rho_bulk*Zrw/MM_K;
    Na_tot(1) = 9*1e-6*rho_bulk*Zrw/MM_Na;

    % initial guess
    f_Ca0 = 0.8*(1 - f_acid);
    f_Mg0 = 0.1*(1 - f_acid);
    f_K0 = 0.05*(1 - f_acid);
    f_Na0 = 0.05*(1 - f_acid);
    f_H0 = f_acid*0.3;
    f_Al0 = f_acid*0.7;
    
    Si(1) = 0;

    equations = @(p)[p(1)*n*Zrw*s(1)*1000+(p(12)/3)*CEC_tot*conv_Al-p(3);
        p(2)-(H(1)^4/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4))*p(1);
        p(4)*n*Zrw*s(1)*1000+p(8)/2*CEC_tot-Mg_tot(1);
        p(5)*n*Zrw*s(1)*1000+p(11)/2*CEC_tot-Ca_tot(1);
        p(6)*n*Zrw*s(1)*1000+p(9)*CEC_tot-Na_tot(1);
        p(7)*n*Zrw*s(1)*1000+p(10)*CEC_tot-K_tot(1);
        p(12) - (p(2)/conv_Al)*(p(11)^3/(K_Ca_Al*p(5)^3))^(1/2);
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
    Al0 = (f_Al0*conv_Al)*((K_Ca_Al*Ca0^3)/f_Ca0^3)^(1/2);
    Al_w0 = Al0/(H(1)^4/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4));
    Al_tot0 = (f_Al0/3)*CEC_tot(1)*conv_Al+Al_w0*n*Zrw*s(1)*1000;

    x0 = [Al_w0, Al0, Al_tot0, Mg0, Ca0, Na0, K0, f_Mg0, f_Na0, f_K0, f_Ca0, f_Al0, f_H0];

    % system solution
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

    if excel_input.opt_TC_biochem == 1
        K(1) = K_in;
        Ca(1) = Ca_in; % [micromol/l]
        Mg(1) = Mg_in;
        Si(1) = Si_in;
    end


elseif excel_input.fX_conc_convention == "f_CEC_to_conc"

    % case study: Dietzen
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CEC fractions
    % f_Ca = 0.05;
    % f_Mg = 0.03;
    % f_K = 0.01;
    % f_Na = 0.01;
    % f_Al = 0.4;
    % f_H = 0.5;

    % estimates of concentrations
    Ca(1) = (f_Ca(1)/K_Ca_H)*(H(1)/f_H(1))^2;
    Mg(1) = (Ca(1)/f_Ca(1))*K_Ca_Mg*f_Mg(1);
    K(1) = ((Ca(1)/f_Ca(1))*K_Ca_K)^(1/2)*f_K(1);
    Si(1) = 0;

    if excel_input.opt_TC_biochem == 1
        Ca(1) = Ca_in; % [micromol/l]
        Mg(1) = Mg_in;
        K(1) = K_in;
        Si(1) = Si_in;
    end

    Na(1) = ((Ca(1)/f_Ca(1))*K_Ca_Na)^(1/2)*f_Na(1);
    Al(1) = ((Ca(1)/f_Ca(1))^3*K_Ca_Al)^(1/2)*f_Al(1)*conv_Al;
    Al_w(1) = Al/(H(1)^4/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4));


elseif excel_input.fX_conc_convention == "conc_to_f_CEC"

    % case study: Kelland
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assign values once get updates from soil samples
    % initial soil ion concentration (e.g., control plots) (mol-conv/l)
    % Ca(1) = Ca_in; % [micromol/l]
    % Mg(1) = Mg_in;
    % K(1) = K_in;
    % Si(1) = Si_in;
    %
    % Na(1) = Na_in;
    % Al_w(1) = Al_w_in; % [nanomol/l]
    Al=(H.^4./(H.^4+H.^3.*K1+H.^2.*K1.*K2+H.*K1.*K2.*K3+K1.*K2.*K3.*K4)).*Al_w;%mol-conv/l

    eqf_Ca=@(f_Ca)(1-(f_Ca+(Al(1)./conv_Al).*((f_Ca.^3./(K_Ca_Al.*Ca(1).^3)).^(1/2))+Mg(1).*(f_Ca./(K_Ca_Mg.*Ca(1)))+Na(1).*((f_Ca./(K_Ca_Na.*Ca(1))).^(1/2))+K(1).*((f_Ca./(K_Ca_K.*Ca(1))).^(1/2))+H(1).*((f_Ca./(K_Ca_H.*Ca(1))).^(1/2))));
    [f_Ca(1), fval, exitflag,~] = fsolve(eqf_Ca,0.2,fsolve_options);
    if exitflag<=0
        error('Equation is not solved for CEC saturation')
    end

    f_Al(1) = (Al(1)/conv_Al)*((f_Ca(1)^3/(K_Ca_Al*Ca(1)^3))^(1/2)); % fraction of aluminium contributing to CEC state
    f_Mg(1) = Mg(1)*(f_Ca(1)/(K_Ca_Mg*Ca(1))); % fraction of magnesium contributing to CEC state
    f_Na(1) = Na(1)*((f_Ca(1)/(K_Ca_Na*Ca(1)))^(1/2)); % fraction of sodium contributing to CEC state
    f_K(1) = K(1)*((f_Ca(1)/(K_Ca_K*Ca(1)))^(1/2)); % fraction of potassium contributing to CEC state
    f_H(1) = H(1)*((f_Ca(1)/(K_Ca_H*Ca(1)))^(1/2));% fraction of hydrogen contributing to CEC state


elseif excel_input.fX_conc_convention == "Kantola"

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT CEC fractions
    f_Mg = (ERW_biogeochem_input.excMg*n*Zrw*s(1)*1000)/CEC_tot(1)*2;
    f_Ca = (ERW_biogeochem_input.excCa*n*Zrw*s(1)*1000)/CEC_tot(1)*2;
    f_K = (ERW_biogeochem_input.excK*n*Zrw*s(1)*1000)/CEC_tot(1);
    f_H = 0.15;

    f_CEC_sat = 0.8;
    f_Na = f_CEC_sat-f_Ca-f_Mg-f_K;
    f_Al = (1-f_CEC_sat)-f_H;

    % estimates of concentrations
    Na(1) = ((Ca(1)/f_Ca(1))*K_Ca_Na)^(1/2)*f_Na(1);
    Al(1) = ((Ca(1)/f_Ca(1))^3*K_Ca_Al)^(1/2)*f_Al(1)*conv_Al;
    Al_w(1) = Al(1)/(H(1)^4/(H(1)^4+H(1)^3*K1+H(1)^2*K1*K2+H(1)*K1*K2*K3+K1*K2*K3*K4));

    % estimates of K_CEC
    K_Ca_H = (f_Ca(1)/Ca(1))*(H(1)/f_H(1))^2;
    K_Ca_Mg = (f_Ca(1)/Ca(1))*(Mg(1)/f_Mg(1));
    K_Ca_K = (f_Ca(1)/Ca(1))*(K(1)/f_K(1))^2;

elseif excel_input.fX_conc_convention == "Kelland"

    % unused
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conc in leachate from untreated experiment (see Table S5)
    Ca(1) = 800; % [micromol/l]
    Mg(1) = 70;
    K(1) = 30;
    Na(1) = 190;
    Si(1) = 200;
    Al_w(1) = 2500; % [nanomol/l]

    Ca_tot(1) = Ca_tot_in;
    Mg_tot(1) = Mg_tot_in;
    K_tot(1) = K_tot_in;
    Na_tot(1) = Na_tot_in;

    Al=(H.^4./(H.^4+H.^3.*K1+H.^2.*K1.*K2+H.*K1.*K2.*K3+K1.*K2.*K3.*K4)).*Al_w;%mol-conv/l

    % CEC
    f_Mg(1) = (Mg_tot(1) - Mg(1)*n*Zrw*s(1)*1000)*2/CEC_tot;
    f_K(1) = (K_tot(1) - K(1)*n*Zrw*s(1)*1000)*2/CEC_tot;
    f_Na(1) = (Na_tot(1) - Na(1)*n*Zrw*s(1)*1000)/CEC_tot;

    equations = @(p)[Al_w(1)*n*Zrw*s(1)*1000+(p(3)/3)*CEC_tot*conv_Al-p(1);
        Ca(1)*n*Zrw*s(1)*1000+p(5)/2*CEC_tot+p(2)-Ca_tot(1);
        p(3) - (Al(1)/conv_Al)*(p(5)^3/(K_Ca_Al*Ca(1)^3))^(1/2);
        p(4) - H(1)*(p(5)/(K_Ca_H*Ca(1)))^(1/2);
        1-(p(5)+p(3)+f_Mg(1)+f_Na(1)+f_K(1)+p(4))];

    % initial guess
    f_H0 = 1e-3;
    f_Al0 = 1e-3;
    f_Ca0 = 0.8*(1 - f_H0 + f_Al0);
    Al_tot0 = (f_Al0/3)*CEC_tot*conv_Al+Al_w(1)*n*Zrw*s(1)*1000;
    CaCO30 = Ca_tot(1)-Ca(1)*n*Zrw*s(1)*1000+f_Ca0/2*CEC_tot;

    x0 = [Al_tot0, CaCO30, f_Al0, f_H0, f_Ca0];

    % system solution
    options = optimoptions(fsolve_options, 'FunctionTolerance', 1e-14);
    sol = fsolve(equations,x0, options);
    Al_tot = sol(1);
    CaCO3 = sol(2);
    f_Al(1) = sol(3);
    f_H(1) = sol(4);
    f_Ca(1) = sol(5);

    % estimating soil-dependent K_CEC
    K_Ca_Mg = (f_Ca(1)/Ca(1))*(Mg(1)/f_Mg(1));
    K_Ca_K = (f_Ca(1)/Ca(1))*(K(1)/f_K(1))^2;
    K_Ca_Na = (f_Ca(1)/Ca(1))*(Na(1)/f_Na(1))^2;

end



K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H];



%% rock_dist_data
if excel_input.rock_filepath ~= 'NaN'
    % particle size distribution - diameter and class weight
    % data from json file
    rock_data = fileread(excel_input.rock_filepath);
    rock_data = jsondecode(rock_data);

    % "Amann" - 'coar' or 'fine'
    % "Pas" - 'olivine','basalt','albite' or 'anorthite'
    % "Kelland" - 'Default Dataset'
    rock_type = excel_input.rock_type;

    feedstock_dataset = rock_data.datasetColl(find(strcmp({rock_data.datasetColl.name}, rock_type), 1));

    x1 = zeros(1, length(feedstock_dataset.data));
    y1 = zeros(1, length(feedstock_dataset.data));
    for a = 1:length(feedstock_dataset.data)
        x1(a) = feedstock_dataset.data(a).value(1);
        y1(a) = feedstock_dataset.data(a).value(2);
    end


    if excel_input.rock_dist_data == "Amann"

        d_in = x1 * 1e-6; % [m]
        psd_in = y1/sum(y1); % [-]

    elseif excel_input.rock_dist_data == "te_Pas"

        d_in = x1 * 1e-6; % [m]
        CDF_in = y1;

        psd_in = diff(CDF_in)./diff(d_in);
        psd_in(psd_in < 0) = 0;
        psd_in = [CDF_in(1)/d_in(1) psd_in];
        psd_in = psd_in/sum(psd_in);

    elseif excel_input.rock_dist_data == "Kelland"

        d_in = flip(x1);
        psd_in = flip(y1);

        d_in = d_in(11:end-1);
        psd_in = psd_in(11:end-1);

        d_in = [10 d_in]* 1e-6;
        psd_in = [100-sum(psd_in) psd_in]/100;

    elseif excel_input.rock_dist_data == "Kantola"

        d_in = x1 * 1e-6; % [m]
        CDF_in = y1;

        d_in(1) = 1e-6;
        CDF_in(1) = 1e-6;

        psd_in = diff(CDF_in)./diff(d_in);
        psd_in(psd_in < 0) = 0;
        psd_in = [CDF_in(1)/d_in(1) psd_in];
        psd_in = psd_in/sum(psd_in);


    end
end

%% save values
ERW_preprocess_output = struct;

excel_var_names = {'f_Ca';'f_Mg';'f_K';'f_Na';'f_Al';'f_H';
    'K_Ca_Mg';'K_Ca_K';'K_Ca_Na';'K_Ca_Al';'K_Ca_H';
    'Ca';'Mg';'K';'Na';'An';'Al';'Si';'H';'Al_w';
    'd_in';'psd_in';'pH'};

for var_index = 1:size(excel_var_names,1)
    ERW_preprocess_output.(cell2mat(excel_var_names(var_index))) = eval(cell2mat(excel_var_names(var_index)));
end

end
