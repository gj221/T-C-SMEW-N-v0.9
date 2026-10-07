function [SAB1_c,SAB2_c,SAD1_c,SAD2_c,PARB_c,PARD_c] = Topographic_Radiation_Projection( ...
    SAB1,SAB2,SAD1,SAD2,PARB,PARD,Datam_i,DeltaGMT,Lon,Lat,t_bef,t_aft, ...
    slope,aspect_deg,svf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topographic_Radiation_Projection
%%
%% Per-cell shortwave for one hour on sloping ground: the direct beam
%% projected onto each cell's slope and aspect (Essery & Marks 2007; the
%% T&C spatial-radiation form T&C-BG-2D uses), the diffuse scaled by the
%% sky-view factor, plus a terrain-reflected share of the global
%% radiation.  One definition, shared by the 2-D driver (every cell, one
%% hour) and the field spin-up's plot-scale replay of a single cell
%% (Lateral_Replay_Prepare), so a replayed cell receives exactly the
%% radiation it received in the field run.
%%
%%   SAB1..PARD      station values for this hour (scalars)
%%   Datam_i         [Y M D H] of this hour
%%   t_bef, t_aft    solar averaging window (TOPO_SOLAR_WINDOW in the 2-D
%%                   driver: the forcing's own hour-ending convention)
%%   slope           [N x 1] tangent of the cell slope
%%   aspect_deg      [N x 1] cell aspect [deg]
%%   svf             [N x 1] sky-view factor
%%
%% Returns [N x 1] per-cell values.  On flat ground (slope 0, svf 1) the
%% beam factor is 1 and the diffuse gain 0, so every value equals the
%% station value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_col = numel(slope);
SvF_c = svf(:);
[h_S_i,~,zeta_S_i] = SetSunVariables( ...
    Datam_i,DeltaGMT,Lon,Lat,t_bef,t_aft);
beta = atan(slope(:));
asp  = aspect_deg(:)*pi/180;
if isfinite(h_S_i) && sin(h_S_i) > 0.10 && ...
        (SAB1+SAB2+PARB) > 0
    cos_fst = cos(beta)*sin(h_S_i) + ...
        sin(beta).*cos(h_S_i).*cos(zeta_S_i - asp);
    beam_factor = max(0,cos_fst)./sin(h_S_i);   %% ~1 on flat ground
else
    beam_factor = ones(N_col,1);
end
diff_gain = (1-SvF_c)*0.15;   %% terrain-reflected fraction of global
SAB1_c = SAB1*beam_factor;  SAB2_c = SAB2*beam_factor;
PARB_c = PARB*beam_factor;
SAD1_c = SAD1*SvF_c + diff_gain*(SAB1+SAD1);
SAD2_c = SAD2*SvF_c + diff_gain*(SAB2+SAD2);
PARD_c = PARD*SvF_c + diff_gain*(PARB+PARD);
end
