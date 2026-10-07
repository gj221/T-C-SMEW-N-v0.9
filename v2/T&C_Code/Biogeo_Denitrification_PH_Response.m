function [gross_factor,n2_to_n2o] = Biogeo_Denitrification_PH_Response(PH)
% Soil-pH factors used by the denitrification scheme.

if ~isreal(PH) || any(~isfinite(PH(:))) || any(PH(:) <= 0 | PH(:) >= 14)
    error('Biogeo_Denitrification_PH_Response:InvalidPH', ...
        'Soil pH must be finite and satisfy 0 < pH < 14.');
end

% Gross response form: Blanc-Betes et al. (2021), Eq. 5,
% doi:10.1111/gcbb.12757.
% Breakpoint pH 7.5. Chosen because it is the upper bound of a directly
% tested range rather than an extrapolation: Olaya-Abril et al. (2021),
% doi:10.1016/j.jprot.2021.104169, measure denitrification proteomics over
% pH 6.5-7.5 and report an optimum of 7.2, and the Simek & Cooper (2002)
% review, doi:10.1046/j.1365-2389.2002.00461.x, places the general soil
% denitrification optimum in the pH 7.0-8.0 band. 7.5 sits inside both.
% The previous value of 7.8 came from soil-specific evidence in Simek et al.
% (2002), doi:10.1016/S0038-0717(02)00059-7, but lies OUTSIDE the range any
% of the cited studies actually tested. Anderson et al. (2017),
% doi:10.1016/j.soilbio.2017.06.002, showed base-cation identity can alter
% the response, so this breakpoint must not be treated as universal.
pH_opt = 7.5;
gross_factor = min(1,max(0.001,(PH-4)/(pH_opt-4)));
% Product ratio: Blanc-Betes et al. (2021), Eq. 9. The increasing pH
% response is also supported by Liu et al. (2010),
% https://academic.oup.com/femsec/article/72/3/407/566642, and Henault et al.
% (2019), doi:10.1038/s41598-019-56694-3.
n2_to_n2o = 0.0016*exp(1.006*PH);
end
