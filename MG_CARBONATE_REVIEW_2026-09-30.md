# Mg carbonate sink, 30 September 2026

Claude session: `session_016HTNKJZpuzjtxMePkGbRBX`. Package m is package l
plus the change below. Package l is the comparison copy.

## What failed in l

The 35-year lower-set basalt run (run 11) stopped in all three SSPs in
2037 with the activity guard in `ERW_Calcite_Kinetics`: gamma_IAP above 1.
The Davies coefficient goes above 1 once I passes 1.94 mol/L, and the
equation is only meant for I up to about 0.5. Medians from the run 11
checkpoint (lower set, ssp126), as contributions to I:

| year | Ca | Mg | HCO3 | CO3 | pH | median I | max I |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2025 | 0.006 | 0.014 | 0.001 | 0.000 | 7.62 | 0.032 | 0.142 |
| 2030 | 0.001 | 0.028 | 0.003 | 0.000 | 8.05 | 0.047 | 0.104 |
| 2031 | 0.000 | 0.047 | 0.013 | 0.003 | 8.68 | 0.076 | 0.386 |
| 2033 | 0.000 | 0.133 | 0.040 | 0.029 | 9.33 | 0.220 | 0.556 |
| 2036 | 0.000 | 0.235 | 0.057 | 0.092 | 9.65 | 0.408 | 0.730 |
| 2037 | 0.000 | 0.229 | 0.044 | 0.112 | 9.82 | 0.407 | 1.668 |

The original and upper run 11 drift the same way (median I 0.45 to 0.54
by 2045) but had not yet crossed the guard. Controls and calcite runs
stay near 0.03. So every basalt 35-year trajectory after about 2031 is
affected, not only the one that stopped.

The cause is a missing Mg sink. Calcite precipitation takes the Ca out by
2031. After that, Mg and carbonate alkalinity from basalt weathering have
nowhere to go but drainage: by 2036 the pore water holds 0.12 mol/L Mg at
pH 9.8. In l, `K_MgCO3` was magnesite (log K -7.46, the PHREEQC minteq.v4
value) and `r_MgCO3` was 0. Magnesite rightly never forms at soil
temperature: it is kinetically inhibited and rarely forms below 60 C
(Harrison et al., 2019), and Hänchen et al. (2008) only got it at 120 C.
The Mg carbonate that does form there is the hydrated one, nesquehonite,
MgCO3·3H2O. At 25 C and 1 bar CO2 only nesquehonite precipitated in
Hänchen et al. (2008). In Harrison et al. (2013) it formed within 2 to 8 h
and "nesquehonite precipitation occurred at equilibrium". At the 2036
composition nesquehonite is about 80 to 150 times supersaturated.

## Change

1. Nesquehonite solubility, per ERW step, overriding the value passed in
   the same way `K_CaCO3` already is:

       log10 K(T) = -5.27 - dH_r/(R ln 10) * (1/T - 1/298.15),  dH_r = -17.17 kJ/mol

   Both numbers are from Harrison et al. (2019), Tables 2 and 3. The
   dissolution is exothermic, so K rises in cold soil. The line gives
   -5.05 at 5 C and -5.37 at 35 C against the measured -5.03 +/- 0.13 and
   -5.34 +/- 0.04 (3 sigma). The parameter file sets the same expression
   at the seed temperature.

2. Saturation on the activity basis already used for calcite:

       Omega = gamma_div^2 * Mg * CO3 / K

3. Rate law, the existing `carb_W` form, now with a nonzero coefficient:

       W = r * Zrw * (1 - Omega)            Omega > 1, precipitation
       W = s * M_s * (1 - Omega) / tau      Omega <= 1, dissolution of the solid M_s

   `r_MgCO3 = r_CaCO3`, the Kirk et al. (2015) soil-volume rate
   5e-10 mol dm^-3 s^-1. No soil rate for nesquehonite is published. `tau_MgCO3`
   stays 40 d.

4. Closed-water bound. `ERW_Calcite_Reaction_Limit` is written for any 1:1
   metal carbonate, so it is reused unchanged: the extent xi is cut back
   to where

       Omega(Mg + xi/V, DIC + xi/V, Alk + 2 xi/V) = 1,   xi >= -min(Mg, DIC) V

   with V = n Zrw s 1000 L m^-2. It runs after the calcite bound, on the
   DIC and alkalinity calcite leaves, in `ERW_function` and in the seed
   step of `ERW_MOD_PARAM_excel`.

5. Mass balance is unchanged. `W_MgCO3` was already booked everywhere it
   belongs:

       Mg_tot(t+dt)  = Mg_tot + (I_Mg + weathering + W_MgCO3 + dMg) dt
       MgCO3(t+dt)   = MgCO3 - W_MgCO3 dt
       IC_tot       += (... + W_MgCO3 + ...) dt
       dAlk_carb    += 2 W_MgCO3 dt

`ERW_Calcite_Reaction_Limit.m` itself is not modified. `ERW_function.m`
gains the solubility line and the second bound. `ERW_MOD_PARAM_excel.m`
changes two constants and gains the same bound for the seed step.

## Assumptions, and how each was checked

**Which phase.** Nesquehonite is metastable. At 25 C it turns into
dypingite within about 57 days, and hydromagnesite is more stable still
(Harrison et al., 2019; Königsberger et al., 1999). The model has one Mg
carbonate pool, so this uses the most soluble phase that actually forms
quickly. That gives the highest dissolved Mg, so it is the conservative
choice for the ionic strength limit.

**Solubility value.** The minteq.v4 (-4.67) and LLNL (-5.06) values are
close to Kline (1929). Harrison et al. (2019) attribute the difference to
aqueous speciation that the older work left out, and recommend their own
values. These are used here. Their Ksp includes a_H2O^3. That is taken as 1,
which is within 0.02 log units below I = 0.5.

**Free ions.** Harrison et al. derived Ksp with the MgCO3° and MgHCO3+ ion
pairs. The model has no ion pairs, for Mg or for Ca. With the phreeqc.dat
constants (log K 2.98 and 11.40), the pairs would hold about 40 % as much
Mg again as the free ion at the capped state. The model therefore caps
total dissolved Mg lower than a fully speciated code would, by about that
fraction. The same simplification is
already in the calcite treatment.

**Rate.** Borrowing the calcite coefficient is an assumption. The box
model below shows that the cap hardly depends on it.

**Operator splitting.** Calcite is bounded first, nesquehonite second, each
from the same start-of-step rates. When calcite has already taken the
carbonate, the Omega handed to the second bound is stale. The bisection
then returns zero extent, and the next step recomputes from the true
state. The test covers that case.

**Bulk-pH extract.** `Porewater_to_Bulk_pH` leaves solid carbonate out of
the 1:2.5 extract; its calcite term is set to 0. Nesquehonite is left out the
same way. In a real extract some of it would redissolve and raise the pH,
so once the pool is large the modelled bulk pH of the basalt soils is a
lower bound. (Package r, 1 October 2026, removes this limit: the paste
dissolves the pool beside the calcite, `PASTE_MG_CARBONATE_2026-10-01.md`.)

**Brucite.** At nesquehonite saturation the brucite ion activity product
reduces to

    log(a_Mg / a_H^2) = log K_nsq - log(K_H K_1 K_2) - log pCO2

which at 25 C is 12.88 - log pCO2. Brucite (log K 16.84, minteq.v4) can
only saturate below pCO2 = 10^-3.96 atm, about 110 ppm, which soil air
never reaches. The box model gives brucite SI -0.6 to -3.6 over 5 to 25 C.
In the drifted l state brucite was supersaturated (SI +0.4). With the
nesquehonite cap that cannot happen.

**When it switches on.** With calcite at saturation, nesquehonite needs

    Mg / Ca > K_nsq / K_cal

which is about 1600 at 25 C and 2000 at 10 C, and higher still when calcite
is undersaturated. Controls and calcite runs are nowhere near that.
Neither are the five-year basalt runs: Mg/Ca in the lower-set run 11 is
still only about 2 in 2025. In those runs `W_MgCO3` stays exactly 0 as in
l, because every workbook case seeds `MgCO3_in = 0`. `Omega_MgCO3` is not
in the saved list, so nothing saved changes. In the 35-year basalt runs
calcite has taken almost all the Ca by 2031, and that is when the drift
begins.

## Dimensional analysis

| Variable | Expected SI Units | Actual Units provided in the code |
|---|---|---|
| T_K | K | K |
| R | J mol^-1 K^-1 | J mol^-1 K^-1 (8.314) |
| dH_r | J mol^-1 | J mol^-1 (17170 in the exponent, sign folded in) |
| dH_r/R * (1/T - 1/298.15) | dimensionless | dimensionless |
| K_MgCO3 | dimensionless (activity product, 1 mol/kg standard state) | umol^2 L^-2, via conv_mol^2; cancels in Omega |
| Mg, CO3, DIC | mol m^-3 | umol L^-1 |
| Alk | mol_c m^-3 | umol_c L^-1 |
| gamma_IAP | dimensionless | dimensionless |
| Omega_MgCO3 | dimensionless | dimensionless |
| r_MgCO3 | mol m^-3 soil s^-1 | umol m^-3 soil d^-1 |
| Zrw | m | m |
| W_MgCO3 | mol m^-2 s^-1 | umol m^-2 d^-1 |
| MgCO3 (solid) | mol m^-2 | umol m^-2 |
| s, n | dimensionless | dimensionless |
| tau_MgCO3, dtw | s | d |
| V = n Zrw s 1000 | m^3 m^-2 | L m^-2 |
| xi = W dtw | mol m^-2 | umol m^-2 |
| k1, k2 | mol L^-1 | umol L^-1 |
| k_w | mol^2 L^-2 | umol^2 L^-2 |

Micromoles, litres and days are used throughout and cancel consistently.
mol/L is used for mol/kg, which is fine for dilute pore water. Two inherited
comments are wrong and are left alone: the `ERW_function` header gives
`W_MgCO3` in mol/m^2/s, and the seed comments say [umol/d] without the
per-area. Neither changes a number.

## Tests

`test_mg_carbonate_precipitation` is new and is in the audit list.

* The bound at the 2036 composition (Mg 0.12 mol/L, CO3 0.056 mol/L, pH 9.8,
  I 0.4). It stops at Omega = 1 within 1e-8. It removes Mg, DIC and
  alkalinity 1:1:2. A stale Omega above 1 on undersaturated water gives
  zero flux.
* `K_MgCO3` through the kernel at 5, 25 and 35 C, with magnesite passed in.
  Each lies inside Harrison's 3 sigma.
* A kernel replay from the stored step 22053 inputs. Undersaturated water
  with no solid gives exactly zero flux. It then runs 12 coupled steps
  from a drifted, Ca-poor state (Omega 156, I 0.456). Omega and I fall at
  every step and gamma_IAP stays at or below 1. Each flux is already inside
  its closed-water bound. Dissolved plus exchangeable plus solid Mg is
  conserved to 1e-9.
* Source pins for both callers and the parameter file.

The test fails on l, at the solubility check, and passes on m.
`test_calcite_event_smoke` now expects `r_MgCO3 == r_CaCO3`.

Local runs used GNU Octave 8.4 with small shims for `contains`, `count`,
`string` and `optimoptions` (lsqnonlin from the optim package). The kernel
replay reproduces the documented MATLAB state of the fixture (pH 11.29).
Every audit test was run on l and m. 51 of the 89 in m run under Octave,
including the new test. Every test shared with l gives the same result on
both. The rest need MATLAB-only functions (datetime, readtable, table,
string arrays, parpool), or the Bouchet baseline, and fail the same way on
l and m. The full audit runs in MATLAB R2025a as part of the Bouchet
validation.

Longer kernel replay from the same drifted state, with no inputs or
drainage. Calcite reaches saturation within the first hours. After that
nesquehonite does the work:

| day | Omega | Mg (mol/L) | I (mol/L) | pH | MgCO3 (mol/m2) |
|---:|---:|---:|---:|---:|---:|
| 0 | 157 | 0.120 | 0.456 | 11.58 | 0 |
| 1 | 109 | 0.103 | 0.386 | 11.48 | 1.6 |
| 2 | 84 | 0.091 | 0.342 | 11.37 | 2.9 |
| 5 | 47 | 0.069 | 0.255 | 10.99 | 5.3 |
| 10 | 23 | 0.050 | 0.178 | 10.40 | 7.4 |

Dissolved plus exchangeable plus solid Mg stays at 36907447.49 umol/m2 to
the last printed digit. The fixture has almost no soil CO2, which is why
the pH stays above 10. What matters is the direction, the conservation and
I back inside the Davies range within days. From an abrupt start like this
the approach takes weeks. In the runs the build-up took years, and the
box model below shows the sink keeps up with that. The same day at 5 C
starts at Omega 103 (nesquehonite is more soluble in the cold) and ends
at I 0.406.

## Box model

An open pore-water box with the model's constants (Plummer and Busenberg
1982 for K1, K2 and KH, Davies A = 0.509), the lower-set background ions
(An 0.018, K 0.008 mol/L) and no ion pairs.

At nesquehonite saturation:

| T (C) | pCO2 (%) | pH | Mg (mol/L) | I (mol/L) | SI brucite |
|---:|---:|---:|---:|---:|---:|
| 5 | 0.04 | 9.35 | 0.019 | 0.067 | -1.68 |
| 5 | 3 | 8.20 | 0.070 | 0.218 | -3.55 |
| 15 | 0.1 | 9.05 | 0.019 | 0.068 | -1.52 |
| 15 | 1 | 8.44 | 0.038 | 0.124 | -2.52 |
| 25 | 0.04 | 9.24 | 0.014 | 0.050 | -0.57 |
| 25 | 3 | 8.11 | 0.045 | 0.143 | -2.45 |

At steady state with all of the Mg supply precipitating and no drainage
(15 C), for the l drift rate of 1.2 mol m^-2 yr^-1 and ten times that:

| supply (mol m^-2 yr^-1) | rate factor | pCO2 (%) | Omega | Mg (mol/L) | I (mol/L) |
|---:|---:|---:|---:|---:|---:|
| 1.2 | 0.1 | 1 | 3.5 | 0.061 | 0.19 |
| 1.2 | 1 | 1 | 1.25 | 0.042 | 0.13 |
| 1.2 | 10 | 1 | 1.03 | 0.039 | 0.13 |
| 12 | 0.1 | 1 | 26 | 0.127 | 0.41 |
| 12 | 1 | 1 | 3.5 | 0.061 | 0.19 |

The box reproduces the drifted l state when the sink is left out (Mg 0.12
mol/L at pCO2 0.09 %: pH 9.72, I 0.43, Omega 83 against the model's 9.65
to 9.8 and 0.41 to 0.47). With the sink, I stays under 0.25 except when
the rate is cut tenfold and the supply is raised tenfold at the same time.

## Failure points

| Failure point | Check |
|---|---|
| Real precipitation is much slower than the calcite rate | Box model: rate x0.1 still keeps I at 0.19 at the l supply |
| The phase is wrong (dypingite or hydromagnesite control Mg) | Those are less soluble, so Mg and I would be lower still; nesquehonite is the upper bound |
| Ion pairs are missing | About 40 % more dissolved Mg in a fully speciated code; I remains below 0.3 |
| Stale Omega after the calcite bound | Zero extent, covered in the test |
| Temperature outside 5 to 35 C | The van't Hoff line extrapolates mildly: -5.00 at 0 C |

## Adversarial edge cases

1. Calcite exhausts the dissolved carbonate (`DIC_final = 0`). The Mg bound
   then has no inventory, `min(Mg, DIC) V = 0`, and the extent is 0. No
   division by DIC happens.
2. Dry soil, s close to 0, so V is close to 0. The calcite bound already
   rejects V <= 0 with `ERW:invalidCalciteLimit` before the Mg bound runs,
   exactly as in l. Evaporative concentration can still push I up through
   other ions. The gamma guard is kept so that case stops loudly rather
   than running with gamma above 1.
3. A solid pool that rounds a few ulps below zero. Dissolution is at most
   s M_s / 40 per day, which is far below M_s / dtw. The pool therefore stays
   positive and the dissolution cap never binds, so a negative dissolution
   flux cannot reach the direction check.

## What the reruns should show

* Controls, calcite runs and the five-year basalt runs match l bit for
  bit in every saved pool and flux, with `MgCO3` at 0 throughout.
* 35-year basalt runs match l until about 2031. After that Mg is held near
  0.01 to 0.07 mol/L, pH near 8.1 to 9.4 and I under about 0.25, and the
  `MgCO3` pool grows. Carbon in that pool is solid inorganic carbon. Two
  equivalents of alkalinity per mole leave solution rather than drain.
  CO2 removal accounting on these runs has to count `MgCO3` beside `CaCO3`.

## Sources

* Harrison, A.L., Mavromatis, V., Oelkers, E.H., Bénézeth, P., 2019.
  Solubility of the hydrated Mg-carbonates nesquehonite and dypingite from
  5 to 35 C. Chemical Geology 504, 123-135. doi:10.1016/j.chemgeo.2018.11.003
* Harrison, A.L., Power, I.M., Dipple, G.M., 2013. Accelerated carbonation
  of brucite in mine tailings for carbon sequestration. Environmental
  Science & Technology 47, 126-134. doi:10.1021/es3012854
* Hänchen, M., Prigiobbe, V., Baciocchi, R., Mazzotti, M., 2008.
  Precipitation in the Mg-carbonate system: effects of temperature and CO2
  pressure. Chemical Engineering Science 63, 1012-1028.
  doi:10.1016/j.ces.2007.09.052
* Königsberger, E., Königsberger, L.-C., Gamsjäger, H., 1999. Low-temperature
  thermodynamic model for the system Na2CO3-MgCO3-CaCO3-H2O. Geochimica et
  Cosmochimica Acta 63, 3105-3119. doi:10.1016/S0016-7037(99)00238-0
* Kirk, G.J.D., Versteegen, A., Ritz, K., Milodowski, A.E., 2015. A simple
  reactive-transport model of calcite precipitation in soils and other
  porous media. Geochimica et Cosmochimica Acta 165, 108-122.
  doi:10.1016/j.gca.2015.05.017
* Plummer, L.N., Busenberg, E., 1982. The solubilities of calcite,
  aragonite and vaterite in CO2-H2O solutions between 0 and 90 C.
  Geochimica et Cosmochimica Acta 46, 1011-1040.
  doi:10.1016/0016-7037(82)90056-4
* PHREEQC databases minteq.v4.dat and phreeqc.dat (Parkhurst and Appelo,
  2013), for the brucite, magnesite and ion-pair constants.
