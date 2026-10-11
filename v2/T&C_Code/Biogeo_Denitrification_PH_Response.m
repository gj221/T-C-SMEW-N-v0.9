function [gross_factor,n2_to_n2o] = Biogeo_Denitrification_PH_Response(PH)
% Soil-pH factors used by the denitrification scheme.

if ~isreal(PH) || any(~isfinite(PH(:))) || any(PH(:) <= 0 | PH(:) >= 14)
    error('Biogeo_Denitrification_PH_Response:InvalidPH', ...
        'Soil pH must be finite and satisfy 0 < pH < 14.');
end

% Gross response shape adapted from Blanc-Betes et al. (2021), Eq. 5,
% doi:10.1111/gcbb.12757. That paper used a pH 7.0 breakpoint.
% The pH 7.8 substitution is soil-specific evidence from Simek et al. (2002),
% doi:10.1016/S0038-0717(02)00059-7. Bremner and Shaw (1958),
% doi:10.1017/S0021859600032779, found rapid activity at pH 8.0-8.6, while
% Anderson et al. (2018), doi:10.7717/peerj.6090, found delayed activity at
% pH 8.8. The breakpoint is a treatment sensitivity, not a universal optimum.
% Olaya-Abril et al. (2021), doi:10.1038/s41598-021-96559-2, is deliberately
% NOT cited for this breakpoint. It reports an optimum of 7.2 for Paracoccus
% denitrificans and tests only pH 6.5-7.5, so it cannot support an optimum at
% 7.8, which lies outside its range; it is also a pure-culture proteomic study
% in mineral medium rather than a bulk-soil rate measurement. Its evidence
% bears on the product ratio below, not on this gross response.
pH_opt = 7.8;
gross_factor = min(1,max(0.001,(PH-4)/(pH_opt-4)));
% Product ratio: Blanc-Betes et al. (2021), Eq. 9. The increasing pH response
% is also supported by Liu et al. (2010),
% https://academic.oup.com/femsec/article/72/3/407/566642, and Henault et al.
% (2019), doi:10.1038/s41598-019-56694-3. Olaya-Abril et al. (2021),
% doi:10.1038/s41598-021-96559-2, gives the mechanism: at pH 6.5 NosZ
% synthesis and intracellular copper are both suppressed, impairing assembly
% of the Cu-dependent N2O reductase, so N2O escapes unreduced. That is the
% process this term represents. Eq. 9 carries no floor of its own.
% The published 0.1 minimum belongs to the complete N2:N2O product in Eq. 6
% and is applied there, in Biogeo_Denitrification_Product_Ratio.
n2_to_n2o = 0.0016*exp(1.006*PH);
end
