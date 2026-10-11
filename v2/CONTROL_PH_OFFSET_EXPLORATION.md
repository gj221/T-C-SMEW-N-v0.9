# Control pH offset: where the base cations go and what could close the gap

Working note, 11 Oct 2026. Branch `bouchet-campaign-gating`. No model edits.
Carries on from the Codex alpha sweep of 11 Oct (`runtime/alpha_20261011`,
not in git). The Si pool mismatch is left out on purpose, it is a known issue.

## What was established before this note

The alpha sweep showed that most of the bulk pH error in the Bouchet 2-D runs
sits in the control, not in the basalt response. Over 2016-2020 the control
loses about 2.59 mol charge m^-2 of Ca+Mg from the solution plus exchange
inventory, and that loss is numerically balanced (charge residuals ~1e-8
umol charge L^-1, exchange fractions sum to one). Alpha only sets how a given
inventory maps onto pH (dpH/dalpha ~ -(f_H - 0.20)); it does not change the
inventory. Deposition as currently implemented adds only ~0.072 mol charge
m^-2 over the five years.

## What the site data say

`v2/Experimental observations/Web plot digitizer/Data_Beerling.mat`, control
plots, pre / 2017 / 2018 / 2019 / 2020:

| | pre | 2017 | 2018 | 2019 | 2020 |
|---|---|---|---|---|---|
| pH 0-10 cm | 6.20 | 6.28 | 6.31 | 6.13 | 5.98 |
| pH 10-30 cm | 6.02 | 6.28 | 6.32 | 6.13 | 6.25 |
| base saturation % 0-10 cm | 79.5 | 80.8 | 85.1 | 79.1 | 72.9 |
| base saturation % 10-30 cm | 81.0 | 77.1 | 76.5 | 76.1 | 82.9 |
| exchangeable Ca (digitised, mg kg^-1) | 3396 | 3226 | 3269 | 3014 | 3524 |

The control was not limed after 2012 (Kantola et al. 2023, GCB), so a
no-lime control in the model is right. But the observed control barely
moves: pH within about 0.2, base saturation within about 7 points at 0-10 cm
and flat at 10-30 cm, exchangeable Ca flat within the digitising noise. With
CEC 19.5 cmolc kg^-1 and 0.3 m of soil at ~1300 kg m^-3 the exchanger holds
~7.6 mol charge m^-2, so a 7 point drop in base saturation is at most ~0.5
mol charge m^-2. The model's 2.59 is five times that. So this is not only
a missing replenishment problem, the loss itself is too large, or the soil
buffers it in a way the model does not represent, or both.

Site management from the two papers, for the parameter check below:
N as 28 % UAN before maize planting, 168 kg N ha^-1 (Kantola 2023) or 202
(Beerling 2024), none in soybean years. Grain harvested by machine. Basalt
(and in the limed subplots, lime) chisel ploughed in each November. Soil pH
on 0-10 and 10-30 cm samples. Exchangeable acidity by a phosphate buffer
method, CEC as bases plus exchangeable acidity.

## How the model moves Ca and Mg (code as of f6c86b2)

Pools (g element m^-2): B60/B68 solution, B61/B69 exchangeable, B62/B70
fixed, B63/B71 primary mineral, B56:59 / B64:67 litter and SOM.
`BIOGEOCHEMISTRY_DYNAMIC3.m` lines 1104-1150 and 1287-1301:

```text
dB60 = Incoming + Ip_lea + decomp - uptake + Kmin_rel + Kfix_sol + Kex_sol - LEAK
LEAK_Ca   = aCa * B60 * min(1, Lk/V)          aCa = aMg = 0.4   (Biogeo_Leakage)
Camin_rel = K3ca * B63                        K3ca = 1/140000 d^-1, K3mg = 1/175000
Cafix_sol = K2ca*B62 - K1ca*B60               zeroed under opt_ERW with 'appelo'
Caex_sol  = Kdca*B61 - Kaca*B60               zeroed under opt_ERW
Incoming  = (FertCa + DepCa)/(Zbio*rsd)       DepCa = DepMg = DepK  (Biogeochemistry_IO.m:101)
```

Under the Appelo scheme `Biogeo_Cation_Partition_Rates` sets all four
exchange and fixation rates to zero, so solution and exchange are
repartitioned only by the daily ERW solve. The ERW receives the day's change
in (B60+B61) and (B68+B69) as a boundary flux (`ERW_Daily_Step_2D.m:332`),
so leaching, uptake, litter return, weathering and deposition all reach the
chemistry through that delta, and exchange itself is internal.

Initial pools from `MOD_PARAM_Kantola_Maize_Soy.m`: B62 = B70 = 0.78 g m^-2,
B63 = B71 = 502 g m^-2. B60/B61 and B68/B69 are overwritten by
`ERW_Initialize_TC_Cation_Pools` from the ERW seed.

Fertiliser: `FertN(81:100) = 0.84` g N m^-2 d^-1 in the maize years, so
168 kg N ha^-1, added as NH4 only (`BIOGEOCHEMISTRY_DYNAMIC3.m:721`).
Residue: `Mpar_L(1).fract_left_AB = 0.5` for maize, `fract_left_fr = 0`.

## Size of each term, from the code constants alone

These are back-of-envelope numbers to rank the leads. The budget script
below replaces them with the retained series.

- Primary mineral release: K3ca*B63 = 3.6e-3 g Ca m^-2 d^-1 and
  K3mg*B71 = 2.9e-3 g Mg m^-2 d^-1, together ~0.15 mol charge m^-2 yr^-1,
  ~0.75 over five years. That is already at the top of what loess-derived
  Mollisols weather (order 0.5-1 kmol charge ha^-1 yr^-1), so there is no
  room to raise it honestly.
- Fixed pools: 0.78 g m^-2 each is 0.04 mol charge m^-2 for Ca+Mg
  together. Even with K2 = 1/13500 d^-1 switched back on the release is
  ~6e-5 g m^-2 d^-1. Dead end, whichever scheme is active.
- Deposition: NADP wet Ca in east-central Illinois is of order 1-2 kg ha^-1
  yr^-1 with Mg a fraction of that; dry deposition roughly doubles it. Call
  it 0.15-0.3 kmol charge ha^-1 yr^-1, 0.1-0.15 mol charge m^-2 over five
  years. The placeholder `DepCa = DepMg = DepK` gives 0.072. Worth replacing
  with the IL11 Bondville numbers for correctness, it will not close the gap.
- Exchange kinetics: cation exchange equilibrates in minutes to hours, so a
  daily equilibrium step is the standard treatment (PHREEQC, SCEPTER do the
  same). Kinetics would shift timing within a day and supply no base. Dead
  end for the inventory; the lead that remains is what the exchanger is, not
  how fast it moves.
- Harvest: with `fract_left_AB = 0.5` half the maize stover Ca and Mg
  leaves the field each maize year. Stover of ~7 t ha^-1 at ~0.25 % Ca and
  ~0.2 % Mg holds ~2 kmol charge ha^-1, so the removal is ~0.1 mol charge
  m^-2 per maize year, ~0.3 over the three maize years. Neither paper
  reports stover removal; the site harvests grain. Setting
  `fract_left_AB = 1` (as for soybean and the prairie margin) removes a loss
  with no field support.
- Leaching: `0.4 * B60 * Lk/V` is 40 % of the drainage-weighted solution
  concentration. B60 is the porewater Ca at the exchange equilibrium, so at
  ~1-2 mmol L^-1 and ~300 mm yr^-1 drainage this is ~0.25-0.5 mol charge
  m^-2 yr^-1 for Ca and Mg together. This is the term that can carry the
  2.59 on its own and it has no site constraint behind the 0.4.

The N side sets the acid load that pushes bases into solution. Net acid
from an ammoniacal fertiliser is about one mol H^+ per mol NO3 leached plus
the ash alkalinity of what is harvested (Helyar and Porter 1989, in Soil
Acidity and Plant Growth; Bolan and Hedley 2003, in Handbook of Soil
Acidity). With the model's 68 kg N ha^-1 yr^-1 nitrate loss (noted at the
KNmaxd calibration as the upper end for tile-drained maize-soy; long-term
Embarras and Upper Salt Fork averages are 31 and 23, David et al. 1997 JEQ
26:1038) that is ~0.5 mol charge m^-2 per maize year from nitrate alone.
Three maize years give ~1.5, which with harvest gets near 2.59. Two things
in the model make it larger than the field: nitrate leaching about twice the
regional record, and all fertiliser N entered as NH4 whereas UAN is one
quarter nitrate and one half urea, so the nitrifiable share is 75 %.

## Candidate changes, ranked

1. Charge-balanced drainage. The ERW porewater is electroneutral at each
   solve. The T&C leak then removes 100 % of NO3 (aNO3 = 1), 100 % of DIC,
   Na and An (ERW side, at concentration times Lw), but only 40 % of Ca and
   Mg. The drained water therefore carries excess anion charge and the
   leached Ca and Mg are tied to their own solution pool rather than to the
   anions actually leaving. Leaching the porewater composition as it stands
   (aCa = aMg = 1 against the ERW solution concentration, which is what
   SCEPTER and any reactive transport code do) removes an untuned 0.4 and
   makes cation loss follow the nitrate and bicarbonate flux. Whether this
   raises or lowers the loss depends on how much of the 2.59 is leaching, so
   run the budget first.
2. Nitrate leaching magnitude. If leaching dominates the budget the fix is
   on the N side: the 68 kg N ha^-1 yr^-1 against ~30 observed in the
   region, and the resin lysimeter N losses measured at the site (Kantola
   2023, 50 cm depth) are the constraint to use. The fertiliser form
   (UAN, not NH4) is a 25 % correction on the nitrifiable input and should
   go in with it.
3. Residue return. `fract_left_AB = 1` for maize. Small, mechanistically
   justified, no downside.
4. Buffering. Even a loss of ~0.5-1 mol charge m^-2 shifts this soil's pH
   by only one or two tenths: Illinois lime tables need roughly 9-11 t ECCE
   ha^-1 per pH unit on dark silt loams (Illinois Agronomy Handbook ch. 8),
   i.e. ~20 mol charge m^-2 per pH unit at plough depth, several times the
   permanent CEC. The difference is the pH-dependent charge of organic
   matter and oxides, which the model does not carry: CEC_tot is fixed at
   19.5 and the only sinks for H^+ are the exchanger and Al(OH)3. A
   variable-charge term, CEC(pH) = CEC_perm + CEC_var*g(pH) with g fitted
   to the site's exchangeable acidity (Beerling 2024 Fig. 1F), would give
   the control the flat pH the field shows without touching the base
   inventory. This is the only lead that also explains why the basalt
   response at alpha 3.4 is near the observed one while the control is not,
   because a stronger buffer damps both. It is also the largest change.
5. Deposition from NADP IL11. Correctness only.

Not worth pursuing: fixed pools, exchange kinetics, raising K3ca/K3mg.

## Dimensional check of the terms used above

| Variable | Expected SI | In the code |
|---|---|---|
| B60, B61, B63 | kg m^-2 | g m^-2; divided by Zbio*rsd inside DYNAMIC3 to g gsoil^-1 |
| LEAK_Ca | kg m^-2 s^-1 | g m^-2 d^-1; aCa dimensionless, Lk/V dimensionless (mm d^-1 over mm) |
| K1ca, K2ca, K3ca | s^-1 | d^-1 |
| DepCa | kg m^-2 s^-1 | g m^-2 d^-1, from a kg ha^-1 yr^-1 grid /10/365 |
| FertN | kg m^-2 s^-1 | g N m^-2 d^-1, 0.84 for 20 d = 16.8 g m^-2 = 168 kg ha^-1 |
| dCa to ERW | mol m^-2 s^-1 | umol m^-2 per window, g/40*1e6 then / window days |
| CEC_tot | mol charge m^-2 | umol charge m^-2: cmolc kg^-1 *1e-5 * g m^-3 * m * 1e6 |
| f_H, alpha | dimensionless | dimensionless |

No unit mismatch found in these paths. The one definitional inconsistency,
Si as all of B76 against aSi, is the known issue and is excluded.

## Three ways this analysis can be wrong

1. The 2.59 may be mostly uptake not returned, not leaching. Then candidate
   1 does little and the question moves to litter Ca return and the 0.9
   lambda that sends litter cations straight to solution.
2. The digitised exchangeable Ca may not be mg kg^-1 (the Mg series is the
   same magnitude as Ca, which is odd for this soil). The base saturation
   series does not depend on that and says the same thing.
3. The buffer estimate from lime tables is for 17 cm of soil and a
   different extraction. If the site's exchangeable acidity in Beerling
   Fig. 1F is small, candidate 4 shrinks.

## What to run next

`v2/T&C_Code/Tests/Kantola_Base_Cation_Budget.m` takes a saved 2-D record
and a column index and returns, in mol charge m^-2, the change in solution
plus exchange Ca+Mg and its split into leaching, plant uptake, harvest
reserve export, primary release and the closed return term (litter and
decomposition return plus deposition). Run it on the alpha 3.4 control for
the full 1826 days:

```matlab
S = load('TC_2D_run01_<tag>.mat');
B = Kantola_Base_Cation_Budget(S, 1, 1826)
```

If `B.leaching` is most of 2.59, candidates 1 and 2 come first. If
`B.net_plant_removal` is, candidate 3 and the litter path come first. In
either case candidate 4 is the one that reconciles the control pH with a
non-trivial loss, and should be sized against the site's exchangeable
acidity before anything is fitted.
