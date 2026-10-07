function cases = Sweep_Cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fixed 1-D Appelo-alpha, KNmaxd and bottom-boundary sweep at dissolution
%% factor 1.
%%
%% Package t (2 October 2026) returns to package o's three grids. Packages
%% p to s swept the dissolution factor (Run_One_Subfield's
%% DissolutionFactor, the workbook diss_f) as a fourth dimension because
%% the basalt run released 12.6 eq m-2 by August 2020 against an exchanger
%% acidity of 10.3 eq m-2. That excess was one rate constant: the amphibole
%% ran on the Palandri-Kharaka tremolite row, about 60 times the measured
%% rate, and t replaces it (AMPHIBOLE_KINETICS_2026-10-02.md). The factor
%% is carried as the single value 1, so tags, provenance and the summary
%% column keep the form of p to s. 112 combinations, 336 five-year
%% trajectories. The reference tuple (alpha 1.3, KNmaxd 0.0048, free
%% drainage) is case 1 to 3, so the one-year control gate is unchanged,
%% and the last case is the extreme pilot. Since package u the all-12
%% campaign runs alpha 3.4 at the same KNmaxd and boundary, cases 193 to
%% 195 (runtime/Base_1D/Campaign_Parameters.m).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = [linspace(1.3,3.4,5) 4.5 6];   % the first five form the earlier grid exactly
KNmaxd = [0.0048 0.0016 0.0032 0.0064];
Kbot = [NaN 0.07 0.03 0.01];
DissolutionFactor = 1;                 % held; the uncertainty ensemble varies it
runs = 10:12;
reference_alpha = 1.3;
reference_logKH = [4.71627429367 3.61627429367 8.10254858734 ...
    8.41854858734 12.918822881];
charge = [1 1 2 2 3];
seed_f_H = 0.20;   % logKH moves hold selectivity at the seed's f_H, keeping 80% base saturation
CEC_cmolc_kg = 13;
n = numel(alpha)*numel(KNmaxd)*numel(Kbot)*numel(DissolutionFactor)*numel(runs);
cases = repmat(struct('id',0,'combination_id',0,'tag','', ...
    'combination_tag','','run_id',0,'AppeloAlpha',0, ...
    'AppeloLogKH',reference_logKH,'CEC_cmolc_kg',CEC_cmolc_kg, ...
    'KNmaxd',0,'Kbot',NaN,'DissolutionFactor',1),1,n);
id = 0;
combination_id = 0;
for ia = 1:numel(alpha)
    for ik = 1:numel(KNmaxd)
        for ib = 1:numel(Kbot)
            for id_f = 1:numel(DissolutionFactor)
                combination_id = combination_id+1;
                combination_tag = sprintf('a%s_kn%s_kb%s_fd%s', ...
                    number_tag(alpha(ia)),number_tag(KNmaxd(ik)), ...
                    kbot_tag(Kbot(ib)),number_tag(DissolutionFactor(id_f)));
                for ir = 1:numel(runs)
                    id = id+1;
                    cases(id).id = id;
                    cases(id).combination_id = combination_id;
                    cases(id).combination_tag = combination_tag;
                    cases(id).tag = sprintf('%s_run%d',combination_tag,runs(ir));
                    cases(id).run_id = runs(ir);
                    cases(id).AppeloAlpha = alpha(ia);
                    cases(id).AppeloLogKH = reference_logKH + ...
                        charge*(alpha(ia)-reference_alpha)*seed_f_H;
                    cases(id).KNmaxd = KNmaxd(ik);
                    cases(id).Kbot = Kbot(ib);
                    cases(id).DissolutionFactor = DissolutionFactor(id_f);
                end
            end
        end
    end
end
end
function tag = number_tag(value)
tag = strrep(sprintf('%.6g',value),'.','p');
tag = strrep(tag,'-','m');
tag = strrep(tag,'+','');
end
function tag = kbot_tag(value)
if isnan(value)
    tag = 'free';
else
    tag = number_tag(value);
end
end
