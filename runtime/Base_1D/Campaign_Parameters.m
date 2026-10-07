function P = Campaign_Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter set for the all-12 Kantola campaign.
%%
%% Alpha 3.4 was adopted on 2026-10-02 (package u) in place of the 1.3
%% fitted on 26 September: it is the Appelo (1994) proton-loading slope,
%% the value SCEPTER carries as its default, and the uncertainty ensemble
%% runs 1.3 and 6 around it (runtime/Uncertainty_1D). KNmaxd 0.0048 and
%% free drainage are kept from the joint 1-D sweep of 26 September
%% (Calibration_Sweep_1D, arrays 27559013 to 27559017): KNmaxd 0.0048
%% reproduced the total chamber N2O at alpha 1.3 (0.97 of observed) and is
%% carried to alpha 3.4 unchanged, and free drainage is the legacy boundary
%% under which Data_Kantola_BG.mat was spun up (no new spin-up). In the
%% sweep's case table the set is combination a3p4_kn0p0048_kbfree_fd1
%% (dissolution factor 1), combination id 65, cases 193 to 195 for runs 10
%% to 12. Every wrapper, validator and assertion of this campaign reads the
%% values from here; no other campaign file spells a number.
%%
%% Alpha is formed exactly as Sweep_Cases forms it, as the last point of
%% linspace(1.3,3.4,5), so that runs 10 to 12 carry the options and result
%% filenames of those sweep cases.
%%
%% CEC (13 cmolc kg^-1) and the Appelo exchange scheme are not options of
%% Run_One_Subfield: they come from the pinned parameter workbook and
%% MOD_PARAM_Kantola_Maize_Soy.m already installed on Bouchet. They are
%% listed here as expectations that the provenance assertion enforces.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_grid = linspace(1.3,3.4,5);
P.AppeloAlpha = alpha_grid(5);
reference_alpha = 1.3;
reference_logKH = [4.71627429367 3.61627429367 8.10254858734 ...
    8.41854858734 12.918822881];
charge = [1 1 2 2 3];
%% Selectivity is 10^(-logKH + alpha*charge*f_H), so shifting logKH by
%% charge*(alpha - reference_alpha)*seed_f_H holds it fixed at the seed's
%% f_H (1 - the measured 0.80 base saturation). The pH 6.1 inversion then
%% keeps 80% base saturation at every alpha, and alpha only sets how the
%% exchanger responds once pH moves.
seed_f_H = 0.20;
P.AppeloLogKH = reference_logKH + charge*(P.AppeloAlpha-reference_alpha)*seed_f_H;
%% The rule's anchors, for a run that moves alpha (the uncertainty
%% ensemble's sensitivity members): logKH = reference_logKH +
%% charge*(alpha - reference_alpha)*seed_f_H, exactly as the sweep forms it.
P.reference_alpha = reference_alpha;
P.reference_logKH = reference_logKH;
P.charge = charge;
P.seed_f_H = seed_f_H;
P.KNmaxd = 0.0048;                       % d^-1
P.Kbot = NaN;                            % free drainage (Leakage_Bottom: Lk = Ko), as spun up
P.CEC_cmolc_kg = 13;                     % expectation: workbook value
P.CECpHSensitivity = 0;                  % expectation: constant CEC
P.exchange_scheme = 'appelo';            % expectation: parameter script
P.NDays = 1826;
P.StartDate = datenum(2016,1,1);         %#ok<DATNM>
P.sweep_campaign = ...
    'Calibration_Sweep_1D';
P.sweep_combination_tag = 'a3p4_kn0p0048_kbfree_fd1';
P.sweep_combination_id = 65;
P.sweep_case_ids = [193 194 195];        % runs 10, 11, 12
end
