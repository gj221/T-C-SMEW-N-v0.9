function sp = Soil_Param_To_Routing(Soil_Param, SPAR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Soil_Param_To_Routing
%%
%% Adapt the Soil_Param struct built by Restating_parameters.m to the
%% field names Conductivity_Suction expects inside Darcy_Lateral_Topology.
%%
%% Two names differ between the two conventions and are the only reason
%% this adapter exists:
%%     Soil_Param.Ks_Zs -> sp.Ks     (saturated conductivity per layer)
%%     Soil_Param.Phy   -> sp.Phy1   (air-dry matric potential, VG only)
%% SPAR is a free-standing workspace variable in the parameter files, so
%% it is passed in rather than read off the struct.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp        = Soil_Param;
sp.SPAR   = SPAR;
sp.Ks     = Soil_Param.Ks_Zs;
sp.Phy1   = Soil_Param.Phy;

required = {'SPAR','Ks','Osat','Ohy','L','Pe','O33','alpVG','nVG','lVG', ...
            'Ks_mac','Omac','alpVGM','nVGM','lVGM','Phy1','s_SVG','bVG'};
missing = required(~isfield(sp, required));
if ~isempty(missing)
    error('Soil_Param_To_Routing:MissingField', ...
        'Soil_Param is missing fields required by Conductivity_Suction: %s', ...
        strjoin(missing, ', '));
end

end
