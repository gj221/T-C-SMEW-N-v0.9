function [MOD_PARAM_PFT_pars] = PFT_selection_parameters(site_info, n_veg, stature_type,crop_type)


par_names = site_info.Properties.VariableNames(1,10:end);
add_term = {'_H','_L'};

%% Initialisation and assign values
for n_par = 1:length(par_names)
    for n_term = 1:2
        if contains(['Sp_LAI'],par_names{n_par})
            par_names_new = [par_names{n_par},add_term{n_term},'_In']; 
    
        elseif contains(['Knit'],par_names{n_par})
            par_names_new = [par_names{n_par},add_term{n_term}(2:end)]; 
    
        elseif matches(["Veg_Optical_Parameter"],string(par_names{n_par}))
            
        else
            par_names_new = [par_names{n_par},add_term{n_term}]; 
        end
    
        eval([par_names_new,'=zeros(1,',char(string(n_veg)),')',';'])

        for n_stat_col = 1:length(stature_type)
            
            if matches(["Veg_Optical_Parameter"],string(par_names{n_par}))
                PFT_type = site_info(n_stat_col,:).Veg_Optical_Parameter;
                eval(['[PFT_opt_',char(stature_type(n_stat_col)),'(',char(string(n_stat_col)),')]=Veg_Optical_Parameter(',char(string(PFT_type)),');'])

            else
                if contains(add_term{n_term},char(stature_type(n_stat_col)))
                
                    eval([par_names_new,'(',char(string(n_stat_col)),')=site_info(', char(string(n_stat_col)),',:).',par_names{n_par},';'])
                
                end

            end

        end

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%

for len_var = 1:n_veg
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% HIGH VEGETATION 
    [Stoich_H(len_var)]=Veg_Stoichiometric_Parameter(Nl_H(len_var),crop_type);
    [ParEx_H(len_var)]=Exudation_Parameter(0);
    [Mpar_H(len_var)]=Vegetation_Management_Parameter;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LOW VEGETATION
    [Stoich_L(len_var)]=Veg_Stoichiometric_Parameter(Nl_L(len_var),crop_type);   
    [ParEx_L(len_var)]=Exudation_Parameter(0); 
    [Mpar_L(len_var)]=Vegetation_Management_Parameter; 
end

vars = who();
MOD_PARAM_PFT_pars = struct();

for ind_vars =1:length(vars)
    MOD_PARAM_PFT_pars.(vars{ind_vars})=eval(vars{ind_vars});
end


end