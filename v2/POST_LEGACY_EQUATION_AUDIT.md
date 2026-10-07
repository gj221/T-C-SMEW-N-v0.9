# Post-legacy equation and mechanism audit

## Scope and acceptance rule

Legacy reference:

`C:\Users\Gregory\OneDrive\Yale\TC-SMEW-N\Other models\T&C-SMEW`

Audited model:

`C:\Users\Gregory\OneDrive\Yale\TC-SMEW-N\T&C-SMEW-N v0.9`

This audit covers equations, state transitions, and boundary couplings that are
absent from the legacy model or materially changed in the current model. A
mechanism is accepted only when all of the following hold:

1. The equation has a physically interpretable state and boundary.
2. Every term has consistent dimensions after the code's documented unit
   conversions.
3. Internal transfers cancel in the relevant element or energy ledger.
4. Source-pool availability is enforced before the state is updated.
5. A focused regression checks the equation, its limiting cases, or its
   conservation identity.
6. Empirical parameters are either traceable to literature or explicitly
   identified as model priors.

Passing a regression does not establish observational validity. That requires
the 92-day initialization gate and fresh full control, basalt, and calcite
trajectories.

## Highest-risk failure points

1. Pore-water and bulk-extract pH have different roles. Pore-water pH drives
   aqueous speciation and mineral kinetics. Reconstructed bulk-extract pH
   drives the T&C microbial response functions. Substituting either pH state
   for the other changes the model and invalidates the campaign fingerprint.
2. Hourly operator-split lateral transport can under-resolve a rapid reversal of
   water flux. It is conservative and refinement-convergent, but it is only
   first-order accurate and must not be used with an enlarged coupling step
   without repeating the refinement test.
3. Mass-action constants are conditioned at the step's ionic strength with
   Davies coefficients, and a step is solved again when its own ionic
   strength moves more than 1 %. Davies holds to about 0.5 mol L^-1 and A is
   kept at its 25 C value, so concentrated fertilizer microsites are outside
   what the activity model covers.

## Equation audit

### 1. ERW water volume and aqueous inventory

For an active soil depth `Zbio` in mm,

```text
W = Zbio * max(0, theta - theta_hyg)
M_i = C_i * W
```

`W` is L m^-2 because 1 mm water equals 1 L m^-2. `C_i` is umol L^-1
and `M_i` is umol m^-2. The previous layer-weighted expression multiplied
layer thickness twice and changed with mesh refinement. It was replaced by
`Biogeo_Mobile_Water_Volume`.

Status: corrected and tested.

### 2. Hydrogen concentration and carbonate speciation

```text
H = 10^(6-pH)                                      [umol L^-1]
h = H / 1e6                                        [mol L^-1]
CO2aq = DIC * h^2 / (h^2 + K1*h + K1*K2)
HCO3  = DIC * K1*h / (h^2 + K1*h + K1*K2)
CO3   = DIC * K1*K2 / (h^2 + K1*h + K1*K2)
```

The acid-base ratios are dimensionless when `h`, `K1`, and `K2` use molar
units. Species sum exactly to DIC. `K1`, `K2` and `Kw` follow soil
temperature (Plummer and Busenberg, 1982) and are conditioned at the step's
ionic strength by `ERW_Conditional_Constants` (`K1/g1^2`, `K2/g2`,
`Kw/g1^2`), so the solved species are concentrations and the reported pH is
`-log10(g1*h)`.

Status: tested (`test_erw_conditional_constants`,
`test_erw_activity_ph_reporting`).

### 3. Al hydrolysis and finite Al(OH)3 inventory

For free Al3+ and hydrogen in the same concentration basis,

```text
Al_total_aq = Al3 * (H^4 + b1*H^3 + b2*H^2 + b3*H + b4) / H^4
Al3_gibbsite = 10^logK * (H/1e6)^3 * 1e6
```

The current initializer inverts these expressions rather than treating free
Al3+ as total dissolved Al. Secondary `AlOH3_solid` is retained as a finite
state, so precipitation cannot create an unlimited Al source and dissolution
cannot overdraw the solid. Crystalline gibbsite uses `logK = 8.04`; amorphous
Al(OH)3 uses `logK = 10.8` only when selected explicitly.

Status: corrected and tested for both exchange-controlled and solid-controlled
initialization.

### 3a. Aluminium leaching scale

The ERW solver carries dissolved aluminium on the `Al_CEC_scale` basis
(`conv_Al = 1e3`, nmol L^-1 against umol L^-1 for every other ion), and
`Al_tot` on the matching nmol m^-2 basis: the seed divides
`Al_tot - (f_Al/3)*CEC*Al_CEC_scale` by the pore-water volume, the hydrolysis
closure forms `AlOH4 = Al_w - (Al + AlOH + AlOH2 + AlOH3)` by subtraction, and
the balance accumulates `EW * Al_CEC_scale`. Leaching removes the mobile
species (free Al3+ and Al(OH)4-) in the drained water,

```text
Leak_Al = max(0, min( (Al + AlOH4) * W * f_Lk,  Al_tot ))
```

with `W` the pore-water volume in L m^-2 and `f_Lk` the daily drained
fraction, the same construction as `Leak_Na` and `Leak_An` in the legacy
steady-background branch both drivers run (`keyword_add = 1`, section 9).
The `keyword_add = 0` branch adds a third cap, `Lw * 1000 * (Al + AlOH4)`.
The first term
is the mobile inventory on `Al_tot`'s basis; when all dissolved Al is mobile
it equals `Al_tot * f_Lk`.

Legacy status: the legacy TC-SMEW line (`Other models/T&C-SMEW`,
`ERW_function.m` lines 263 and 271) multiplies that first term by `1e-3`,
which puts it on the umol basis against the other two terms and holds
aluminium leaching about a thousand times below the drainage every other
species sees. The factor was inherited verbatim into this model. Removing it
is a departure from the legacy equation and is recorded here as one. It
affects every run through the aluminium inventory, the exchangeable-Al
fraction and the Al(OH)3 solid target; at this site's pH range the dissolved
pool is small, so the effect is expected to be modest and is measured by a
1-D control/basalt/calcite comparison against the activity-coefficient
state (`tc_erw_al_leaching_scale_20260922_v1`).

Status: corrected and tested (`test_erw_al_leaching_scale` pins both source
lines, replays the seed and closure on the solver scale, and checks that the
balance drains `f_Lk` of the mobile inventory per day on the weathering
basis).

### 3b. Carbonate saturation on the ion activity product

`k1` and `k2` are conditioned at the step's ionic strength (section 2), so
the solved `HCO3` and `CO3` are concentrations, like `Ca`. The calcite ion
activity product is then

```text
a_Ca * a_CO3 = (gamma_div*Ca) * (gamma_div*CO3)
Omega = gamma_div^2 * Ca * CO3 / K_CaCO3(T)
```

with `gamma_div` the Davies coefficient at the step's ionic strength.
`K_CaCO3(T)` is the Plummer and Busenberg (1982) calcite solubility,

```text
log10 K_CaCO3 = -171.9065 - 0.077993 T + 2839.319/T + 71.595 log10(T)
```

with T in K (-8.48 at 25 C, -8.41 at 10 C), taken at the step's soil
temperature and at the seed temperature for initialization. The fixed
`10^-8.35` it replaces sits near aragonite and held calcite Omega 7 to 26 %
low. Every step records its own `Omega_CaCO3`.

An earlier revision formed the product with `gamma_div/gamma_mono^2`,
because its speciation used thermodynamic constants with concentrations.
That factor no longer applies once the constants are conditioned. This is
the same construction SCEPTER 1.0.2 uses (per-reaction `fkeq`, pH reported
as an activity).

Legacy status: departure. The legacy chemistry carried no activity
coefficients and a fixed `K_CaCO3 = 10^-8.35`.

Status: tested (`test_erw_omega_activity_basis`, `test_calcite_*`).

### 4. Gaines-Thomas exchange and base saturation

The exchange state is solved on charge-equivalent fractions. The complete
state satisfies

```text
f_H + f_Al + f_Na + f_K + f_Ca + f_Mg = 1
base_saturation = f_Na + f_K + f_Ca + f_Mg
exchange_acidity = f_H + f_Al
```

Selectivity equations use one internally consistent concentration basis. The
Al/Ca inversion now receives `f_Al`, not `f_Ca`. Cation totals used by the
T&C/ERW daily coupling include solution plus exchange inventory, so internal
solution-exchange transfer cancels.

A change of the Appelo alpha moves each full-reaction `logKH` by
`Z*(alpha - 1.3)*0.20`. Selectivity is `10^(-logKH + alpha*Z*f_H)`, so the
shift holds it fixed at the seed's `f_H = 0.20` and the pH 6.1 inversion keeps
the measured 80 % base saturation at any alpha (0.798 at 1.3 and at 2.35).
The earlier shift, `Z*(alpha - 1.3)`, held selectivity fixed at `f_H = 1`
and let the seed fall to 69 % base saturation at alpha 2.35.

Status: corrected and tested for fraction closure, cation closure, pH response,
finite Al, and base saturation at the campaign alpha (`test_kantola_seed_charge`).

### 5. Bulk-pH extraction operator and biological pH

Pore-water pH is the aqueous chemistry state. Bulk pH is reconstructed from
the same finite aqueous and exchange totals using the configured 1:2 DI-water
extraction operator:

```text
pH_bulk = Porewater_to_Bulk_pH(finite totals, exchange state)
pH_offset = 0
```

The Kantola initial condition uses the measured pretreatment bulk pH as the
boundary for `Invert_Bulk_pH`. This solves the initial pore-water hydrogen and
exchange state without applying the same pH observation twice. Later bulk pH
is forward-calculated from the evolving field state. It is retained for
comparison with soil-extract observations and is also the pH supplied to the
T&C decomposition, nitrification and denitrification response functions.
Pore-water pH remains the input to aqueous speciation, exchange and mineral
kinetics.

This is a finite extraction calculation, not an additive calibration. The
extract uses 25 C laboratory constants conditioned at its own ionic strength
and reports an activity pH. It is still an idealized laboratory operator,
because equilibration time is not represented. Undissolved primary and
secondary calcite are left out of the extract, following SCEPTER's
laboratory extraction; the calcite branch of `Porewater_to_Bulk_pH` is
switched off on purpose.

Status: initialization, forward extraction and the separate aqueous and
biological pH call paths are tested.

### 6. Mineral saturation index

For a reaction with stoichiometric exponents `nu_i`,

```text
IAP = product(a_i^nu_i)
Omega = IAP / Ksp
SI = log10(Omega)
```

Each charged ion enters as its Davies activity at the step's ionic strength,
proton included; neutral H4SiO4 takes none. Corrected stoichiometry includes:

```text
albite:       NaAlSi3O8
anorthite:    CaAl2Si2O8
zoisite:      Ca2Al3Si3O12(OH)
epidote:      Ca2Al2FeSi3O12(OH)
hydroxyapatite: Ca5(PO4)3OH
```

Unsupported Fe-bearing or custom phases return `NaN` saturation rather than a
fabricated number. Their dissolution can retain a kinetic affinity of one only
when the caller explicitly treats them as irreversible primary minerals.

Albite and anorthite are written as their Al+++ reactions, so their
constants are the llnl.dat log K (2.7645 and 26.578) rooted per proton. The
constants used before belonged to the legacy kaolinite-product reactions and
put albite Omega 23 times too high. Albite and quartz still read
supersaturated in the current runs, because ERW receives the whole T&C Si
pool as dissolved Si (inherited from legacy).

Status: stoichiometric defects corrected and tested against direct IAP
reconstruction and an llnl.dat equilibrium composition
(`test_erw_mineral_kinetics`).

### 7. Primary-mineral shrinking-particle kinetics

For each particle class,

```text
r = r0 * (M/M0)^(1/3)
A = A0 * (M/M0)^(2/3)
dM/dt = -A * k(T,pH) * affinity
```

The particle-size distribution is formed from CDF mass increments, not CDF
slopes. Primary mineral mass is finite and dissolution is limited by remaining
inventory before update. Native mineral units, molar conversion, and g m^-2
are kept separate.

Status: corrected and tested for particle distribution, molar-mass shape,
finite inventory, and dissolution direction.

### 8. Calcite and carbonate solids

Primary calcite only dissolves. Dissolved carbonate can precipitate into a
separate secondary `CaCO3` pool. Dissolution follows the Plummer-Wigley-Parkhurst
surface-reaction structure. Secondary calcite precipitation uses the Kirk et al.
soil-volume rate coefficient:

```text
alpha = 5e-10 mol dm^-3 soil s^-1
      = 43200 umol m^-3 soil d^-1
Rprecip = alpha * f(Omega), Omega > 1
```

Direct anhydrous MgCO3 precipitation is disabled at ambient field temperature;
its dissolution pathway remains available. Precipitation and dissolution are
limited so `Omega` cannot be crossed by more material than is chemically or
physically available. Calcite solubility follows soil temperature
(section 3b).

Status: corrected and tested for density, DIC stoichiometry, event routing,
primary inventory, secondary separation, flux sign, and Ca closure.

### 9. T&C to ERW daily element coupling

The first retained ERW state is already synchronized to the first T&C state.
It is not charged with the raw seed difference again. For later retained days,

```text
DeltaX_j = X_TC(j) - X_TC(j-1)
rate_X = DeltaX_j / Delta_t_window
sum(rate_X * dt_window) = DeltaX_j
```

Each daily delta is applied once. A repeated final partial flush carries zero
new T&C delta. K, Ca, and Mg totals include solution plus exchange pools;
therefore exchange is an internal redistribution rather than a boundary flux.

Background Na and residual anion follow the legacy steady closure
(`keyword_add = 1`) in both drivers. Each day `I_Na` and `I_An` return that
day's leaching and uptake loss, so both hold their seed inventories. The 1-D
frame used to switch this off for Al-pinned runs while the 2-D driver kept it
on, and the equivalence test copied the 1-D setting into both, so it never ran
the production branch. It now checks that its baseline used the closure
production reads from the workbook.

Status: duplicate-final-delta and partial-window defects corrected. Focused
coupling and 1-D/2-D equivalence tests cover the boundary, on the production
background closure.

### 9a. Weathered phosphorus return

```text
P_release_step   = sum(nu_P .* EW_step) * dtw            [umol P m^-2]
P_release_window = sum(P_release_step) over the integrated rows
P43_after = P43_before + P_release_window*31/1e6         [g P m^-2], aP < 1
P43_after = PO4_tot_end*31/1e6                           [g P m^-2], aP = 1
```

`nu_P` is `min_st(:,7)` in mol P per mol mineral and `EW` is the mineral
dissolution rate in umol mineral m^-2 d^-1, so the product is umol P m^-2 d^-1
and `dtw` converts it to an areal increment. Multiplying by the molar mass of
elemental P and dividing by 1e6 gives g P m^-2, the unit of the T&C mineral-P
state `P(:,43)`. Layer thickness does not enter this conversion. Weathered P
also enters the alkalinity balance once, as `dAlk_PO4_w = -P_release`.

The summed rows are `max(dm_start,2):dm_end`, exactly the row set the ERW
window integrates. The carried anchor row `max(dm_start,2)-1` is therefore
never added twice, and a resumed or repeated final window carries no new
source. `ERW_Return_Weathered_Phosphorus` recomputes the expected inventory
`aP*P43_before/31*1e6 + P_release_window` and raises `ERW:phosphorusBudget`
when the returned `PO4_tot` does not match, so an inconsistent window is an
error rather than a silent gain. When `aP = 1` the dissolved pool is the whole
pool, and the operator writes that level instead of adding the increment again.

Both the twelve-cell runs and the vertical 1-D runs reach this operator through
the single `Run_One_Subfield` to `MAIN_FRAME_BG_ERW_2D` to `ERW_Daily_Step_2D`
path, so one owner performs the transfer in both layouts. The legacy
single-column `MAIN_FRAME_BG_ERW_10mins` frame used by `v2/run/ERW_run.m` does
not carry the return and must not be used for phosphorus comparisons.

Every feedstock row selected by `Run_Matrix_12` has `min_st(:,7) = 0`, so the
current campaign returns exactly zero basalt-derived P. The Blue Ridge
metabasalt row lists albite, ferroactinolite, epidote, chlorite, quartz and
calcite after Lewis et al. (2021), doi:10.1016/j.apgeochem.2021.105023, which
reports that apatite, the P-bearing phase in the other studied basalts, was not
detected in Blue Ridge, and that no K-bearing mineral was detected in it
either. Zero P and zero K release is faithful to that measured mineralogy
rather than a missing term. The residual limitation is that the rock still
carries bulk P and K that XRD did not resolve into a named phase. This model
releases nutrients only from named mineral stoichiometry, so P held in glass or
in undetected accessory phases is not released. An apatite fraction may be
added only from a feedstock measurement that assigns P to a phase. It must not
be fitted to the field outcome.

Status: operator verified for units, single transfer, resume safety, budget
closure and the P-free boundary. The zero-P feedstock is a documented
representation limit of the current campaign, not a defect in the transfer.

### 9b. Dissolved cation export and the background salt closure

The box exports what its drainage water carries. T&C leaches each dissolved
pool with its drained fraction of the active-layer water,

```text
LEAK_X = min(B_X, a_X * B_X * min(1, Lk/VT))        [g m^-2 d^-1]
```

where `a_X` is the fraction of the pool in solution (Porporato et al., 2003).
`P(60)` and `P(68)` hold only dissolved Ca and Mg since the ERW module took
over the exchange partition (section 9), so `a_Ca = a_Mg = 1`, as for NO3
and K. The legacy 0.4 belonged to a lumped pool; under it the leachate
carried 2.5 times less Ca+Mg charge than anion charge, 84 % of the alkalinity
the rock released stayed in the box and the exchanger saturated by 2019
(`LEACHING_CLOSURE_REVIEW_2026-10-01.md`). The same coefficient scales
passive uptake and the 2-D lateral mobility.

The steady closure of the residual anion and Na inventories (`keyword_add =
1`) returns each day what leached or was taken up. The ERW books its own An
and Na fluxes, so the cycle is alkalinity-neutral; but with Ca and Mg
leaving at their concentration the anion charge Na did not carry out leaves
on them, and the bare anion returning would be an acid input. The leached
part therefore returns with its cations through the existing boundary slots
(`ERW_Leached_Salt_Resupply`, both drivers):

```text
2 I_Ca + 2 I_Mg = f_Lw * max(0, An_tot - Na_tot)   [umol charge m^-2 d^-1]
I_Ca : I_Mg = Ca : Mg (pore water; equal when both are zero)
```

booked in `dAlk_boundary`. The uptake part keeps the legacy booking.

Status: corrected and tested (`test_leached_salt_closure`,
`test_lateral_routing`). Every pin and the equivalence baseline change.

### 10. Gross nitrification pH response

```text
f_nit(pH) = clamp(0.56 + atan(pi*0.45*(pH-5))/pi, 0, 1)
f_T_nit = ((70-T)/(70-38))^12 * exp(12*(T-38)/(70-38))
f_water_nit = 4*Se*(1-Se)
R_nit = k_nit * f_T_nit * f_water_nit * f_nit(pH_bulk) * NH4
```

`k_nit = 0.086 d^-1`. The pH factor is dimensionless, bounded and monotonic.
The moisture response is zero at `Se = 0` and `Se = 1`, with a maximum of one
at `Se = 0.5`. The total rate remains first-order in the internal NH4 pool
`B(31)`, matching the Dickinson et al. (2002) and legacy T&C-SMEW reservoir
loss. It uses the same-day reconstructed bulk-extract pH, not pore-water pH.

`f_nit(pH)` is the Parton et al. (1996) Figure 2c curve transcribed exactly,
`0.56 + atan(pi*0.45*(-5+pH))/pi`. That curve is a fit to 13 tropical soils
spanning pH 4.0 to 6.8 (Motavalli et al., 1996) plus Gilmour (1984), and the
authors caution it may not suit acid soils with high nitrification. The
treatment pH of this campaign sits above the fitted range, but the function has
effectively saturated by its upper calibration limit: 0.94 at pH 6.8 against
0.98 at pH 8. The extrapolation therefore carries little of the campaign
signal.

Parton et al. (2001) adds a mineralization term, `FNO3 = Net_min*K1 +
Kmax*NH4*F(t)*F(WFPS)*F(pH)` with `K1 = 0.20` and `Kmax = 0.10 d^-1`. This
model retains the legacy TC-SMEW first-order NH4 form and implements only the
second term, with `k_nit = 0.086 d^-1` against Parton's `Kmax = 0.10 d^-1`.
Because the omitted term is pH-independent, all nitrification here passes
through `f_nit(pH)`, so total nitrification is more pH-sensitive than DAYCENT
or CLM5. Any pathway that depends on pH-driven NH4 drawdown, including the
indirect pH response of volatilization, is correspondingly stronger.

Status: source coefficient, temperature, moisture, pH, NH4 dependence and the
shared NH4 availability cap are tested.

### 11. Nitrification N2O fraction

```text
f_N2O_nit = min(1, 721.86 * exp(-2.387*pH))
N2O_nit = R_nit * f_N2O_nit
```

The fraction is dimensionless and cannot exceed gross nitrification. This is
the Val Martin et al. implementation used for the Energy Farm comparison.
The NH4 pool loses the complete gross nitrification flux. The nitrate pool
receives `R_nit - N2O_nit`, so emitted N2O-N and nitrate production close to
gross nitrified N.

The primary source is the Inatomi et al. (2019) meta-analysis, adopted into
CLM5 by Val Martin et al. (2023) in place of the constant `K2` fraction that
Parton et al. (2001) and earlier CLM5 versions used. Val Martin et al. evaluate
it at cropland pH 5.8 to 6.2, where it returns 3 to 7e-4 and is close to the
6e-4 constant it replaced, and they state that further work is needed under
other soil conditions. At pH 7.8 to 8 the function returns about 4e-6. Unlike
the nitrification pH response, this function is still falling steeply across
the extrapolated interval, so a large part of any modelled N2O reduction
originates outside the range in which the function was evaluated. This is a
property of the published function, not of its transcription here.

Status: tested for bounds, monotonic bulk-pH response, donor limitation and N
partition closure.

### 12. Gross denitrification and N2:N2O partition

```text
DOC_v = (B(12) + B(13)) * rho_soil                      [g C m^-3]
f_DOC = DOC_v / (DOC_v + Kc), Kc = 17.5 g C m^-3
f_den_pH = clamp((pH_bulk - 4)/(7.5 - 4), 0.001, 1)
f_T_den = exp(308.56*(1/68.02 - 1/(T + 46.02)))
f_water_den = Se^2
R_den_potential = KNmaxd*f_T_den*f_water_den*f_den_pH*f_DOC*NO3
R_den = min(NO3_available_rate, R_den_potential)

f_ratio_pH = 0.0016*exp(1.006*pH_bulk)
f_NO3 = (1 - (0.5 + atan(pi*0.01*(NO3_ppm - 190))/pi))*25
f_resp = 13 + 30.78*atan(pi*0.07*(Resp_kg_C_ha_d - 13))/pi
p_DEN = 13^(2.2*WFPS)
f_WFPS = 1.4*13^(-17/p_DEN)
ratio_N2_to_N2O = max(0.1, min(f_NO3,f_resp)*f_WFPS*f_ratio_pH)
N2O_den = R_den / (1 + ratio_N2_to_N2O)
N2_den  = R_den - N2O_den
```

`KNmaxd = 0.0048 d^-1` is the campaign fit and a reference rate at 22 C,
where `f_T_den = 1`; the legacy Dickinson et al. (2002) value is 0.217 d^-1.
`f_water_den` is the legacy T&C `Se^2`. `B(12)`, `B(13)` and `B(32)` are internal g element per
g dry-soil pools. Multiplication by dry bulk density converts the DOC pool to
the volumetric basis required by `Kc`. Layer depth must not appear in this
half-saturation conversion.

The total nitrate dependence is the first-order Dickinson et al. (2002) and
legacy T&C-SMEW reservoir loss. Li et al. (1992) use `Kn` inside a double-Monod
denitrifier growth model that also carries denitrifier biomass, yield and
maintenance. Those states are absent here. `Kn = 82.5 g N m^-3` is retained as
an inactive provenance field; multiplying a nitrate Monod fraction by the
legacy reservoir loss would create the unsupported `k*N^2/(N+Kn)` response.
Setting a total-rate `Vmax = k*Kn` would introduce an uncalibrated capacity and
is not used. The later PnET-DNDC table reports `Kc = 17 g C m^-3` and
`Kn = 83 g N m^-3`. The inherited TC-SMEW values are 17.5 and 82.5 on that
same volumetric basis. This correction fixes the unit basis without silently
retuning either value.

Blanc-Betes et al. (2021) used the same clipped linear gross-pH shape with a
pH 7.0 breakpoint. The active pH 7.5 breakpoint is the upper end of a range
that was actually tested: Olaya-Abril et al. (2021) measured denitrification
over pH 6.5 to 7.5 with an optimum at 7.2, and Simek and Cooper (2002) put the
general soil optimum between pH 7.0 and 8.0. The earlier 7.8 came from
soil-specific evidence in Simek et al. (2002) and lies outside any tested
range. Alpha 2.35 and KNmaxd 0.0048 were first fitted while 7.8 was active;
the refit of 26 September, which chose alpha 1.3 and kept KNmaxd 0.0048, ran
at 7.5. The
breakpoint is not a universal optimum and must not be described as the
unmodified Blanc-Betes equation.

The nitrate, respiration and WFPS terms were checked against Parton et al.
(1996) Figure 5 and its text. The `min` is the published structure, not a
transplant from the rate equation: Weier et al.'s data led the authors to
control the combined nitrate and respiration effect on the ratio by Liebig's
law of the minimum, with each function giving the maximum ratio while the other
substrate is non-limiting. `f_NO3`, `f_resp` and `f_WFPS` reproduce Figures 5b,
5c and 5a with their published coefficients.

Retaining the Parton terms rather than adopting Blanc-Betes Eqs. 7 and 8 is a
positive choice, not a default. Blanc-Betes contribute two pH scalars, Eq. 5 on
the gross rate and Eq. 9 on the product ratio; Eqs. 7 and 8 are the Del Grosso
et al. (2000) ratio machinery they inherited unchanged and did not validate.
They state that the ability of Fr(pH) to reproduce the pH response of the
N2:N2O ratio was tested independently against direct observations from McMillan
et al. (2016) and Mukumbuta et al. (2018), so Eq. 9 is a separable pH scalar
rather than a coefficient that is only meaningful beside Eqs. 7 and 8. Taking
Eqs. 5 and 9 onto the legacy TC-SMEW ratio is therefore the same operation
Blanc-Betes performed on DayCent, applied to a different host model. Adopting
Eq. 7 would instead import a diffusivity-fitted coefficient `k` that scales the
product ratio and is not calibrated here, and adopting it coherently would
require Eqs. 1 to 4, replacing the first-order nitrate reservoir this campaign
established. Their validation is also site-specific evidence for the retained
choice: adding the pH response to denitrification cut Energy Farm N2O RMSE from
27.7 to 10.3 for the Blueridge feedstock used here.

The N2:N2O ratio is a hybrid and is documented as one. `f_ratio_pH` is
Blanc-Betes Eq. 9 unchanged, and Eq. 9 carries no floor of its own. The
nitrate, respiration and WFPS terms are the legacy TC-SMEW DayCent forms from
Parton et al. (1996) and Del Grosso et al. (2000); the Blanc-Betes Eqs. 7 and 8
substitutions are not used, because the legacy form is retained wherever the
replacement is not required. The published 0.1 minimum in Blanc-Betes Eq. 6
applies to the complete product, so it is applied there. Flooring `f_ratio_pH`
alone was not equivalent: it held the pH term up below pH 4.1, and it still
allowed a dry or nitrate-rich cell to drive the product to zero and report
denitrification as pure N2O. The composite must not be described as the
published Blanc-Betes ratio.

The Parton nitrate, heterotrophic-respiration and WFPS functions affect only
the N2:N2O product split. They do not change the gross nitrate sink. N2 and
denitrification-derived N2O sum exactly to gross denitrification. Total saved
N2O also includes the separately tracked nitrification-derived N2O flux.

Status: tested for the first-order nitrate response, volumetric DOC
half-saturation, temperature, moisture, pH, availability, partition closure and
nonnegative gas fluxes.

### 12a. Empirical ammonia volatilization

```text
VOL = kv * NH4, kv = 8.6e-5 d^-1
```

`VOL` is first-order in the internal NH4 pool `B(31)` and is converted to
g N m^-2 d^-1 with the other reported N fluxes. This retains Dickinson et al.
(2002) and legacy T&C-SMEW. There is no explicit TAN pool, NH3/NH4 equilibrium,
surface resistance or direct pH multiplier. pH changes volatilization only
indirectly because pH-responsive nitrification competes for NH4 and changes the
NH4 inventory available on later days. The shared availability cap prevents
nitrification plus volatilization from overdrawing NH4.

Status: numerically conservative legacy loss diagnostic. It is not a complete
mechanistic NH3 volatilization model.

### 13. Biological N fixation carbon payment

```text
BfixN = Cfix_today / Costfix(T)
```

Both numerator and cost use g C m^-2 d^-1 and g C per g N, so the result is
g N m^-2 d^-1. The current model uses the same-day fixation allocation
`Rexmy(3)`, not its 365-day moving average. This prevents fixation after the
carbon payment stops. `Costfix(T)` uses the Fisher et al. (2010) coefficients
(`s = -6.25`, `c = 25.15`, 7.5 to 12.5 g C g^-1 N) in both places it is
evaluated: `Root_Exudation`, where the plant sizes the carbon it pays, and
`Biogeo_Bio_fixation`, where the soil credits the fixed N. Legacy used
`s = -30` (36 to 60 g C g^-1 N) in both. A version that changed only the soil
side credited 4.8 times the N the plant paid for. The legacy model does not
retain a separate symbiont-respiration pool, which remains a structural
limitation.

Status: current boundary corrected and tested, with one cost on both sides
(`test_biological_fixation_carbon_allocation`).

### 14. Plant reserve and crown-area ledgers

The reserve is signed mass, never a sentinel: below zero it records tissue
diluted below its nominal content (legacy T&C). Initialization uses zero
unless an explicit reserve is supplied. Per-crown fluxes convert to plot-area
mass exactly once:

```text
F_plot = Ccrown * F_crown
Delta reserve = inputs - remobilization - export
```

Harvest export is accumulated in plot-area units. The state update writes only
the current retained row; it does not overwrite the previous row. Standing
leaf nutrient closure uses one concentration for prior storage, same-day
senescence input, litter loss, and export.

At harvest the whole remaining nutrient reserve leaves the field through
`AgrHarNut`, as in legacy, and none of it goes to litter. That matters
because the legacy planting reserve (20 g N m^-2 for maize, 15 for soybean) is
far above real seed N.

Growth may dilute the non-structural tissue to 0.65 of its nominal content,
as in legacy (`rCa_min` and `rMg_min`, both 0.65, for Ca and Mg). For each
element X, with `X_i` the organ C:X ratios and `f_i` the allocation:

```text
TXns   = B1/Xl + B3/Xr + B4/Xc + B5/Xf               g X m^-2
rXc    = 1 + Xreserve/TXns        (Xreserve < 0)     -
Xavail = max(Xreserve + 0.35*TXns, 0)                g X m^-2
Xcons  = sum(f_i*NPP/X_i) + transfers                g X m^-2 d^-1
f_red  = (Xavail/dtd)/Xcons       (Xcons > Xavail/dtd)
NPP   <- f_red*NPP,  WResp <- WResp + (1-f_red)*NPP  g C m^-2 d^-1
```

The elements act in turn, so the scarcest one sets growth, and
translocation stops once any availability is zero. A deficit larger than
`TXns` would be negative plant N and is rejected. On a day of negative NPP
the reserve pays leaf dark respiration (`Tr_l`), which replaces the
leaf-out translocation, roots included; legacy kept the leaf-out `Tr_r`
there, root carbon the reserve never paid. At harvest the organs leave at
the concentration the plant has that day, `rXM = 1 + min(Xreserve,0)/Xharv`
with `Xharv = dtd*(RB1/Xl + RB3/Xr + RB4/Xc + RB5/Xf)`, so a deficit leaves
with the grain and the residue and the reserve export is zero to round-off.

Status: legacy sentinel, crown-zero division, row-overwrite, and harvest-unit
defects corrected and tested; harvest reserve export and tissue dilution
restored to legacy (`test_crop_harvest_reserve_ledger`,
`test_crop_harvest_scale_reserve_export`, `test_plant_tissue_dilution`,
`test_kantola_n_magnitude_audit`).

### 15. Lateral Darcy water flow

For edge half-distances `d_a` and `d_b`,

```text
Q = A_face * (h_a-h_b) / (d_a/K_a + d_b/K_b)
Delta theta_c = sum(signed Q * dt) / (A_cell * Zlayer)
```

`Q` is m^3 h^-1. The perimeter is closed. Every face is shared by two cells,
so water gain and loss cancel exactly. Adaptive substeps limit donor and
receiver water without clipping volume. On a flat field, elevation head
cancels and the lateral gradient is matric head only. `K_a` and `K_b` are
horizontal conductivities, `K*aR`, as in T&C's own lateral flow (`aR = 1` at
Kantola); `aR` was previously dropped.

Status: tested for direction, unequal-cell resistance, closure, bounds,
heterogeneous soils, applied-volume identity and anisotropy.

### 16. Frozen-pore impedance

```text
F_ice = theta_ice / (theta_liquid-mobile + theta_ice)
K_frozen = K_liquid * 10^(-5*F_ice)
```

Lateral and vertical flow both use the legacy T&C form (`HYDROLOGIC_UNIT`),
`10^(-5*F_ice)` with the Hansson et al. (2004) ice fraction written above. The
lateral solver used `Omega = 6` from Swenson et al. (2012) before, which let
lateral and vertical frozen conductivity differ by up to `10^(-F_ice)`.
Liquid plus ice cannot exceed mobile pore volume.

Status: tested; lateral and vertical impedance now agree.

### 17. Lateral solute advection

```text
Delta M_i = Q * dt * C_i,donor
```

`Q*dt` is m^3 and `C_i` is mass or amount per m^3, so `Delta M_i` is mass or
amount. Export is limited by donor inventory. The beginning-of-step donor
concentration makes this a first-order operator split. A fixed-total-time test
shows monotonic convergence at 1, 2, 4, and 64 substeps.

Status: conservative and refinement-convergent. No lateral dispersion or
molecular diffusion is represented.

### 18. Lateral heat and freeze-thaw enthalpy

```text
H = cv*T - rho_i*Lf*theta_ice
dtheta_liq/dT = Lf * (1000*Cw) / (g*Tref)
Capp = dH/dT
     = cv + (rho_i*Lf + T*(cv_w-cv_i))*dtheta_liq/dT
E_adv = rho_w*cp_w*Q*dt*T_donor
E_cond = G_edge*(T_a-T_b)*dt
```

`Cw` is the water-retention derivative in mm^-1. `Lf/(g*Tref)` is m K^-1;
the factor 1000 converts it to mm K^-1. Thus `dtheta_liq/dT` is K^-1 and
`Capp` is J m^-3 K^-1. The solver advances physical enthalpy and inverts it
with safeguarded Newton/bisection. `Capp` is only the derivative and stability
capacity; an invalid value now raises an error instead of silently reverting to
sensible heat capacity. A fixed offset carries T&C's post-cryosuction phase
departure without creating latent heat.

Status: equation, derivative, water-phase closure, energy closure, stiff-edge
bounds, and off-equilibrium state are tested.

### 19. ERW lateral aqueous and soil-gas transport

Aqueous Na, Al species, residual anion, and DIC use the same donor-volume
advection as other solutes. Soil-gas CO2 uses diffusive edge exchange, and
the closed lateral diffusion term cancels, so the lateral step conserves
inorganic carbon. Air pushed out or drawn in by lateral water is booked once,
by the `ADV` term of the next ERW step, which sees the same saturation change.
An earlier piston boundary booked that exchange a second time. Its atmosphere
residual is gone from the driver and the campaign check; lateral inorganic
carbon is checked through the DIC entry of the lateral mass residual.

Status: tested for species transport, saturation response, concentration
bounds, closed diffusion, and inorganic-carbon closure of the lateral step.

### 20. 2-D topology and parameter adapters

The topology partitions one square domain into main and margin cells while
retaining exact cell areas, face widths, half-cell distances, and a closed
outer perimeter. Parameter adapters map legacy field names without changing
values. They do not introduce a new physical equation.

Status: geometry and adapter contracts tested. The 60 m square with a 5 m
margin is a model domain, not the literal Kantola experimental plot geometry.

### 21. Prairie-margin parameterization

`PFT_Prairie_Margin.m` uses legacy photosynthesis, hydraulics, phenology,
allocation, turnover, and nutrient equations. It introduces parameter priors,
not new governing equations. Big bluestem rooting beyond 2 m is literature
plausible, but exact root percentiles, coneflower rooting, phenology thresholds,
hydraulic vulnerability, turnover, initial biomass, and tissue stoichiometry
are not Kantola measurements.

Status: experimental parameterization. It is not used by 1-D runs 10-12 and
must not be presented as site-validated without margin observations and a
separate sensitivity analysis.

## Dimensional analysis

| Variable | Expected SI Units | Actual Units provided in code |
|---|---:|---:|
| `Zbio` | m | mm, converted through 1 mm = 1 L m^-2 |
| `theta`, `theta_hyg`, exchange fractions | m^3 m^-3 or 1 | dimensionless |
| `W` active water | m^3 m^-2 | L m^-2, numerically equal to mm |
| aqueous `C_i` | mol m^-3 | umol L^-1, numerically equal to mmol m^-3 |
| areal aqueous `M_i` | mol m^-2 | umol m^-2 |
| `H`, cations, anions | mol m^-3 | umol L^-1 |
| `h`, `K1`, `K2` in carbonate ratios | mol m^-3 on one basis | mol L^-1 after `/1e6` conversion |
| DIC species | mol m^-3 | umol L^-1 |
| alkalinity | eq m^-3 | ueq L^-1 or umol charge L^-1 |
| CEC inventory | eq m^-2 | umol charge m^-2 after soil-depth conversion |
| mineral `Ksp`, `IAP`, `Omega` | dimensionless standard-state ratio | Davies activities at the step's ionic strength |
| calcite `K_CaCO3` | dimensionless | (umol L^-1)^2, Plummer and Busenberg (1982) at soil temperature |
| primary mineral mass | kg m^-2 | native mineral amount plus g m^-2 diagnostics |
| mineral amount | mol m^-2 | umol m^-2 internally |
| mineral surface area | m^2 m^-2 | m^2 mineral m^-2 ground |
| mineral rate constant | mol m^-2 mineral s^-1 or source-specific equivalent | converted to umol m^-2 ground d^-1 |
| calcite precipitation `alpha` | mol m^-3 soil s^-1 | 5e-10 mol dm^-3 soil s^-1, converted to 43200 umol m^-3 d^-1 |
| T&C daily `DeltaX` | kg or mol m^-2 d^-1 over retained day | g m^-2 daily deltas converted to umol m^-2 |
| N pools | kg N m^-2 | g N m^-2 |
| internal `B(31)`, `B(32)` | kg N kg^-1 dry soil | g N g^-1 dry soil, dimensionally equivalent |
| internal `B(12) + B(13)` | kg C kg^-1 dry soil | g C g^-1 dry soil, dimensionally equivalent |
| dry bulk density `rsd` | kg dry soil m^-3 | converted internally to g dry soil m^-3 |
| `DOC_v` | kg C m^-3 soil | g C m^-3 soil |
| active `Kc` | kg C m^-3 soil | 17.5 g C m^-3 soil |
| inactive `Kn` | kg N m^-3 soil | 82.5 g N m^-3 soil; not used in the rate |
| `kv`, `KNmaxn`, `KNmaxd` | s^-1 | d^-1 |
| N response factors | 1 | dimensionless |
| N process rates | kg N m^-2 s^-1 | g N m^-2 d^-1 |
| reported N2O event flux | kg N m^-2 s^-1 | kg N ha^-1 d^-1 after diagnostic conversion |
| fixation carbon allocation | kg C m^-2 s^-1 | g C m^-2 d^-1 |
| fixation cost | kg C kg^-1 N | g C g^-1 N |
| crown fraction | 1 | dimensionless plot-area fraction |
| Darcy head `h` | m | mm, divided by 1000 in edge gradient |
| hydraulic conductivity `K` | m s^-1 | mm h^-1, divided by 1000 for m h^-1 |
| lateral face area | m^2 | m^2 |
| lateral water flux `Q` | m^3 s^-1 | m^3 h^-1 |
| lateral solute concentration | mol m^-3 or kg m^-3 | pool-specific amount m^-3 |
| lateral solute transfer | mol or kg per step | pool-specific amount per cell per step |
| temperature `T` | K difference | degC relative to 0 degC; increments equal K |
| `cv`, `Capp` | J m^-3 K^-1 | J m^-3 K^-1 |
| `Cw` retention derivative | m^-1 | mm^-1 |
| `dtheta_liq/dT` | K^-1 | K^-1 after explicit mm conversion |
| thermal conductivity `lambda` | W m^-1 K^-1 | W m^-1 K^-1 |
| edge thermal conductance `G` | W K^-1 | W K^-1 |
| cell enthalpy density `H` | J m^-3 | J m^-3 |
| advective/conductive heat | J | J |
| soil-gas CO2 | mol m^-3 air | umol L^-1 air |
| gas diffusion transfer | mol per step | umol per cell per step |

No unresolved dimensional mismatch remains in the accepted post-legacy
equations. The previous areal interpretation of `Kc` and `Kn` was removed.
The concentration-as-activity approximation, empirical volatilization term and
prairie parameters are scientific scope limitations, not unit mismatches.

## Boundary conditions and assumptions

- Vertical T&C hydrology, energy, vegetation, and legacy soil C-N-P processes
  remain the base model.
- The 1-D domain has no lateral exchange. The 2-D perimeter is closed to
  lateral water, solute, soil gas, and conductive heat. Soil air exchanges
  with the atmosphere vertically, through the ERW step's `ADV` and surface
  flux.
- Pore-water chemistry uses Davies activities with ionic-strength-conditional
  constants.
- Background Na and residual anion hold their seed inventories through the
  legacy steady closure (`keyword_add = 1`); the leached part returns as a
  neutral salt (section 9b).
- Pore-water pH drives aqueous and mineral reactions. Reconstructed bulk pH
  drives T&C microbial response functions and remains the soil-extract
  diagnostic. These call paths must remain separate.
- Nitrification, volatilization and denitrification share explicit NH4 and NO3
  availability limits. A pH response cannot authorize more N loss than the
  available donor pool.
- The ERW step's NH4 and NO3 are copies of T&C's pools, stepped linearly
  through the day, and feed only its chemistry; below zero they are
  rounding and read as zero.
- Ammonia volatilization is the retained first-order NH4 diagnostic. Direct
  NH3 speciation and atmosphere-transfer resistance are outside this model.
- Nitrate leakage is drainage-driven mass flow with an inventory cap. pH
  changes it only through the nitrate-producing and nitrate-consuming rates.
- Primary mineral inventories are finite. Secondary carbonates and Al(OH)3 are
  distinct retained pools.
- Internal speciation, exchange, sorption, mineral transfer, biological
  immobilization, and plant-litter transfer cancel in element ledgers.
- Fertilizer, deposition, fixation, gas loss, leakage, lateral transfer,
  amendment, harvest, and export are explicit boundaries.
- A state clamp is not accepted as conservation. Source availability is
  enforced before update, and infeasible chemistry raises an error. The one
  exception is rounding: a pool the day's limited sinks empty ends at zero
  when the sum of its terms lands within 16 eps of it (per gram of soil in
  T&C, the band the donor limits treat as nothing), and the uptake cap stops
  at the headroom exactly. Anything beyond that band is an error.

## Test-to-mechanism map

| Mechanism | Focused tests | Acceptance identity |
|---|---|---|
| Active water | `test_biogeo_mobile_water_volume` | depth scaling and nonnegative volume |
| Appelo exchange and Al | `test_appelo_*`, `test_kantola_seed_charge` | charge, fraction, and finite Al closure |
| Aluminium leaching scale | `test_erw_al_leaching_scale` | mobile inventory on the solver Al basis, drained fraction closes the Al_tot balance |
| Carbonate saturation basis | `test_erw_omega_activity_basis` | ion activity product formed from the code own CO3; Davies z^2 identity |
| Bulk/pore pH | `test_pH_conversion`, `test_kantola_bulk_pH_calibration` | initialization, forward extraction and separate reaction call paths |
| Saturation and mineral kinetics | `test_erw_mineral_kinetics`, `test_erw_particle_distribution` | IAP stoichiometry and finite inventory |
| Calcite | `test_calcite_*` | DIC/Ca stoichiometry, event routing, primary/secondary separation |
| ERW/T&C coupling | `test_erw_exchange_pool_coupling`, `test_2d_driver_equivalence` | daily delta exactly once and 1-D/2-D identity |
| N transformations | `test_nitrification_*`, `test_denitrification_*`, `test_nitrogen_process_shapes` | source coefficients, first-order nitrate response, volumetric DOC limit, pH response and N2+N2O closure |
| Plant reserves | `test_crop_*`, `test_plant_*`, `test_biological_fixation_carbon_allocation` | plot-area reserve and export closure |
| Weathered P return | `test_erw_phosphorus_return`, `test_current_12_run_mineralogy_zero_p` | one transfer per window, areal unit match and zero release from P-free mineralogy |
| Pool availability | `test_pool_overdraw_guards`, `test_no3_immobilization_replay` | no donor overdraw or mass-creating clamp |
| Lateral flow and solutes | `test_lateral_routing` | water/solute closure and timestep convergence |
| Lateral heat and ERW gas | `test_lateral_energy_erw` | physical enthalpy, phase, IC, and elemental closure |
| Initialization boundaries | `test_erw_initial_boundary_defaults`, `test_erw_equilibrium_initialization` | no arbitrary source or unconstrained seed |
| Adapters and resume safety | `test_model_adapter_contracts`, `test_checkpoint_*` | exact mapping and source fingerprint coverage |

## Primary literature used

- Appelo, 1994, cation exchange and transport:
  https://agupubs.onlinelibrary.wiley.com/doi/10.1029/94WR01048
- Palandri and Kharaka, 2004, mineral kinetic constants:
  https://pubs.usgs.gov/publication/ofr20041068
- Hermanska et al., 2022, amphibole rate law (ferroactinolite, actinolite):
  https://doi.org/10.1016/j.chemgeo.2022.120807
- Plummer, Wigley, and Parkhurst, 1978, calcite dissolution kinetics:
  https://ajsonline.org/article/59968-the-kinetics-of-calcite-dissolution-in-co-2-water-systems-at-5-degrees-to-60-degrees-c-and-0-0-to-1-0-atm-co-2.pdf
- Kirk et al., 2015, reactive-transport model of soil calcite precipitation:
  https://doi.org/10.1016/j.gca.2015.05.017
- May et al., 1979, gibbsite equilibrium:
  https://doi.org/10.1016/0016-7037(79)90224-2
- Kanzaki et al., 2024, SCEPTER soil weathering chemistry:
  https://gmd.copernicus.org/articles/17/4515/2024/gmd-17-4515-2024.pdf
- Dickinson et al., 2002, first-order NH4 volatilization, nitrification and
  nitrate-denitrification reservoir losses:
  https://doi.org/10.1175/1520-0442(2002)015%3C0278:NCOCME%3E2.0.CO;2
- Li et al., 1992, volumetric double-Monod denitrifier-growth parameters:
  https://doi.org/10.1029/92JD00509
- PnET-DNDC parameter and denitrifier-state table:
  https://www.dndc.sr.unh.edu/papers/2000_JGR_PnETDNDC1.pdf
- Parton et al., 1996, soil N gas process representation:
  https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/96GB01455
- Val Martin et al., 2023, Energy Farm N-cycle and N2O response:
  https://gmd.copernicus.org/articles/16/5783/2023/gmd-16-5783-2023.pdf
- Parton et al., 2001, revised nitrification rate and NOx:N2O partition:
  https://doi.org/10.1029/2001JD900101
- Del Grosso et al., 2000, denitrification rate and N2:N2O ratio functions:
  https://doi.org/10.1029/1999GB001225
- Inatomi et al., 2019, meta-analysis of the nitrification N2O fraction:
  https://doi.org/10.1371/journal.pone.0219159
- Liu et al., 2010, soil pH control on N2O reductase and the N2:N2O ratio:
  https://academic.oup.com/femsec/article/72/3/407/566642
- Blanc-Betes et al., enhanced weathering and N2O field evidence:
  https://eprints.whiterose.ac.uk/id/eprint/169796/1/gcbb.12757.pdf
- Simek et al., 2002, soil-community-specific denitrification pH response:
  https://doi.org/10.1016/S0038-0717(02)00059-7
- Dall'Amico et al., 2011, conservative frozen-soil enthalpy:
  https://tc.copernicus.org/articles/5/469/2011/tc-5-469-2011.pdf
- Swenson et al., 2012, frozen-pore hydraulic impedance:
  https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2012MS000165
- Nippert et al., 2012, tallgrass root depth and water use:
  https://doi.org/10.1007/s11104-011-1112-4

## Bouchet validation and production gate

The package includes this audit. The coordinated controller,
`runtime/Base_1D/submit_all_12.sh`, creates this dependency graph under
MATLAB R2025a:

```text
R2025a ten-day baseline
    -> 1-D validation: 6 CPUs, 40 GB, four process workers
       control, basalt, calcite and repeated calcite for 315 days
    -> 2-D validation: array element 1, 6 CPUs, 32 GB

both validations pass
    -> 1-D production allocation: three exclusive one-CPU steps
    -> 2-D production array: elements 1-9, at most nine active
```

The two validations become eligible independently after the baseline and may
run concurrently. Both production allocations wait for both validation job
IDs. Once that barrier clears, runs 1 to 9 and runs 10 to 12 are all eligible.
Slurm may run all twelve together, but it may also stagger them according to
available resources. Eligibility is not evidence that a run started or
completed. Production acceptance still requires `COMPLETED 0:0` and the final
PASS line from each validation.

## Evidence still required before production acceptance

1. Four clean executions of the post-legacy equation test set from one source
   state.
2. A cold-start 92-day control/basalt/calcite gate with treatment equality,
   all-pool finiteness, exchange/charge/element closure, and a complete N
   ledger.
3. Fresh full 1-D trajectories crossing application day 314, with matched
   audits at days 180, 314, 400, and 679.
4. Comparison of bulk pH, base saturation, and exchangeable acidity with
   Beerling et al.; N2O and mineral N with Val Martin et al.; and GPP, ecosystem
   respiration, NEE, and biomass with Kantola et al.
5. A separate calibration and validation exercise before the prairie-margin
   parameterization is treated as site-specific evidence.

The zero basalt-derived phosphorus supply is accepted campaign scope, not an
open acceptance item. The Blue Ridge feedstock carries no P-bearing mineral
phase, which section 9a documents against Lewis et al. (2021) and which
Blanc-Betes et al. (2021) corroborate independently by describing Blueridge as
the lower-P, higher-alkalinity basalt. Results from these twelve runs therefore
carry no P-mediated response by construction. A phase-resolved feedstock
measurement would be required only if a basalt-derived P supply were later
claimed.
