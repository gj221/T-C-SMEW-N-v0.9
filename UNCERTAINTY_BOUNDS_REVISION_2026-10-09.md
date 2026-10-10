# Uncertainty bounds revision (package ac)

Package ac (10 October 2026) implements section 2 of this note in
`runtime/Uncertainty_1D/Uncertainty_Members.m`, pinned by
`test_uncertainty_campaign`; `runtime/Uncertainty_1D/README.md` carries the
member table and `RUNBOOK_2026-09-24.md`, section 0ac, the 1-D pre-check of
section 3 (`precheck_1d.m` at the package root prints the checks). The note
below is the record of the decision as written on 9 October.

---

# Uncertainty bounds for the 108-run rerun (decided 9 Oct 2026, to apply next week)

Greg's instruction, 9 Oct 2026, after the first 5-year validation figures
(`Downloads\kantola_figures_5yr\validation_5yr_*.pdf`, package u results):

> Alter uncertainty bounds
> 1. Much reduce Knmaxd
> 2. Alter alpha to ensure uncertainty bounds cover upper pH bounds and
>    overlap with observations, rather than the uncertainty region drop
>    below the main line
>
> determine the more suitable bounds based on literature and heuristics,
> but so we don't get really unreasonable uncertainty regions but that the
> bounds can be justified. Make sure there isn't any unrealistic overshoot
> with uncertainty regions.

The 108 runs are the three parameter sets (original, lower, upper) over the
36 configurations (5 and 35 years, 1-D and 2-D, standard and 5 t/ha
calcite arms, three SSPs). They are to be rerun once the 35-year figures
have been looked at. Members live in
`runtime/Uncertainty_1D/Uncertainty_Members.m`; `test_uncertainty_campaign`
pins every value and must change with them.

## 1. What the first run actually did (1-D, runs 10-12, 2016-2020)

Read from `Downloads\kantola_extract_5yr` (`scratchpad\member_readout.m`).
The figure's band is the envelope of the three sets; the CSV's
"lower"/"upper" columns are min/max, not the member names. By member:

| control run 10 | original | lower set | upper set |
|---|---|---|---|
| N2O 2016-2020, kg N/ha | 24.0 | 6.4 | 149 |
| N2O in 2016 (soybean, no fertiliser) | 1.7 | 0.7 | 36 |
| peak daily N2O, g N/ha/day | 81 | 35 | 401 (2016) |
| NO3 leaching, kg N/ha/yr (2017-2019) | 61-92 | 8-35 | 194-249 (145 in 2016) |
| NEE 2016, g C/m2 | -258 | -348 | +794 |
| NEE 2016-2020 | -617 | -1288 | +462 |
| bulk pH, 5-yr mean / August 2020 | 5.86 / 5.73 | 6.04 / 6.00 | 5.07 / 5.13 |
| basalt effect, August 2020 | +0.39 | +0.21 | +0.74 |

Observations: chambers 1.6 (2017) and 7.0 (2018) kg N2O-N/ha over the
campaign windows; control pH 6.07-6.31 in August; basalt effect +0.47 by
2020; tower NEE -523, -416, +79, -292 g C/m2 for 2017-2020.

Diagnosis:

- The upper set is dominated by a spin-up transient, not by KNmaxd or
  alpha. All twelve decomposition literals sit at the top of Table S2
  (Vm x5.3, Vpc x4, Vpl x34, Km x3.8, rem/rep x4.2, pepf x5.2, mrb x3,
  mrf x6, fd x2, Vdf x3) on pools spun up with the original values, so the
  first year mineralises about a decade of litter: +865 g C/m2 of net
  respiration, 145 kg N/ha of nitrate leached and 36 kg N2O-N/ha emitted
  before any fertiliser, and the nitrification acidity takes the control
  from 6.1 to 4.9 within months. Everything after 2016 runs on that
  damaged state. KNmaxd 0.01 can contribute at most x2 (the rate is first
  order in KNmaxd, capped by nitrate availability); the other x3 and the
  pH collapse come from the literals.
- The lower set is mild despite moving the same literals 5-80x the other
  way (Vm x0.012, Vpc x0.024): decomposition slows, the system becomes
  sink-limited, N2O falls to 0.27x and pH stays flat at 6.04 (no acid
  drift because little nitrate is produced or leached).
- So the band is asymmetric by construction: the only member above the
  central pH line is the lower set (+0.18), and its basalt effect is the
  smallest (+0.21), because alpha 1.3 sits in it. Alpha 6, which gives the
  largest basalt effect (+0.74), sits in the acidified member where it
  cannot lift the band. That is why "the uncertainty region drops below
  the main line".
- Table S2's ranges are the cross-site literature priors of the MEND-type
  microbial model (Wang et al. 2013, Waring and Averill 2013), spanning
  two to three orders of magnitude. They describe ignorance across soils,
  not the uncertainty that remains after the model has been calibrated to
  this field's tower NEE, respiration and chambers. Taking every one of
  them to its end at once (a corner of a 16-dimensional box) is a
  parameter set no soil has, and it is not spun up.

## 2. Proposed bounds

Principle: the central values stay; each range is the smallest one that
can be defended from data or literature; the two moved sets stay
spin-up-consistent (no first-year burst); and the members are assigned so
that the band straddles the central line for the quantity the paper
reports, the bulk-pH response to basalt.

### 2a. KNmaxd (Greg: "much reduce")

| | lower | original | upper |
|---|---|---|---|
| now | 0.002 | 0.0048 | 0.01 |
| proposed | 0.002 | 0.0048 | 0.0064 |

Justification: 0.0048 d-1 is calibrated to N2O/obs = 1 on the multi-year
cumulative, but it hides two compensating per-year errors recorded in
`Biogeochemistry_Parameter.m`: 2017 modelled 2.4x the chambers, 2018
0.75x. The per-year fits therefore bracket the constant at 0.0048/2.4 =
0.002 and 0.0048/0.75 = 0.0064. Both ends have literature support:
Hou et al. (2024) measure first-order soil denitrification constants of
0.002-0.003 d-1 in the upper 40 cm; CropSyst's reference constant is
0.01 d-1 at 15 C (Heinen 2003; Heinen 2006 reviews the first-order
family), so 0.0064 is comfortably inside the modelled range and 0.01 is
dropped as the top of a cross-model bracket rather than a value fitted
to this site. 0.0064 was already run in the package-t sweep grid
(kn0p0064), so its exchange behaviour is known. Expected effect on the
band: N2O at the upper end about 1.3x the central instead of 6x once the
literals are tamed.

### 2b. Alpha (Greg: "alter alpha so the band covers the upper pH")

| | lower | original | upper |
|---|---|---|---|
| now | 1.3 | 3.4 | 6 |
| proposed | 1.3 | 3.4 | 6 (constants slid by the campaign rule, unchanged) |

Alpha's range is defensible as it stands; what must change is which
member carries which end (section 2d). Anchors: 1.3 is the fitted
mesocosm value (26 Sep); 3.4 is Appelo's (1994) proton-loading slope and
SCEPTER's default (Kanzaki et al.); the field's own pH-versus-base-
saturation pairs (Beerling 2024 / Pihlap, 0.063-0.080 pH per BS point,
against the model's 0.056 at 3.4) correspond to alpha 4.1-5.8 at the
seed's f_H = 0.20, so 6 brackets the field slope with a margin and is
the value under which the basalt effect reaches +0.57 to +0.74, above
the field's +0.47. Narrowing the top to 5.8 would be data-exact but
gains nothing and forfeits the sweep's existing alpha-6 runs.

What alpha cannot do: it sets pH per base-saturation point for the
control's acid drift and for the basalt gain alike (0.035 / 0.056 /
0.082 per point at 1.3 / 3.4 / 6). A high alpha therefore raises the
basalt plot and lowers the control plot in the same member. The band's
upper edge for the control is set by the nitrogen cycle (slower
mineralisation, more denitrification), not by alpha; the basalt plot's
upper edge needs high alpha in that same slow-nitrogen member. The
remaining gap to the observations after the change (control drift of
0.2-0.3 and a basalt effect about 0.8x the field's at alpha 3.4) is the
model's known bias, to be reported as such, not absorbed into the band
([[accept-residual-ph-discrepancy-decision]]).

### 2c. The twelve biogeochemistry literals

Replace the Table S2 ends with a factor-of-two bracket around the shipped
values, clipped to Table S2 where the bracket would leave it. A factor
of two is the usual post-calibration bracket for microbial kinetic
constants and keeps the potential decomposition within 2x of the state
the pools were spun up with, which the first run shows is the condition
for no burst (x4-34 burst; x0.01-0.4 did not).

| parameter | shipped | lower | upper | Table S2 |
|---|---|---|---|---|
| Vm | 100 | 50 | 200 | 1.2-528 |
| Vpc | 200 | 100 | 400 | 4.8-792 |
| Km | 0.25 | 0.125 | 0.5 | 0.057-0.95 |
| rem | 0.018 | 0.009 | 0.036 | 0.0075-0.075 |
| rep | 0.018 | 0.009 | 0.036 | 0.0075-0.075 |
| pepf | 0.006 | 0.0031 (S2 floor) | 0.012 | 0.0031-0.031 |
| mrb | 0.010 (legacy) | 0.005 | 0.020 | 0.0012-0.03 |
| mrf | 0.005 (legacy) | 0.0025 | 0.010 | 0.0012-0.03 |
| fd | 0.4 | 0.3 | 0.5 | 0.2-0.8 |
| Vpl | 23 | 11.5 | 46 | 4.8-792 |
| Vdb | 0.04 | 0.02 | 0.06 (S2 ceiling) | 0.0024-0.06 |
| Vdf | 0.02 | 0.01 | 0.04 | 0.0024-0.06 |

fd is a fraction, so +/-0.1 rather than x2. If the 1-D check (section 3)
still shows a first-year respiration excess above about 150 g C/m2 or an
unfertilised-year N2O above 5 kg N/ha, fall back to x/÷1.5.

The rigorous alternative is to spin each moved set up with its own
literals (`Run_BG_Spinup.m`, three phases, one 1-D spin-up per member)
and start the members from their own `Data_Kantola_BG` state. That makes
even wide ranges spin-up-consistent, but it is a new build (member-
specific initial states carried through the provenance checks) and does
not remove the objection that Table S2's ends are not this site's
uncertainty. Not recommended for next week's rerun.

### 2d. Dissolution factor and saturation exponent

Keep 0.8 / 1 / 1.2 and s^2 / s / sqrt(s) (the legacy analysis, Zhang et
al. SI). The CDR band they give is 1.3-2.1 t CO2/ha around 1.6 (+/-25 %),
with no overshoot; the distance to Beerling's 4-11 t CO2/ha is the
known export/mass-balance gap, not a bound problem.

### 2e. Member assignment: by direction of the basalt pH response

Now both moved sets are the ends of the parameter ranges, so parameters
that push pH opposite ways share a member and cancel, and the band does
not straddle the central line. Assign each parameter's end by the sign
of its effect on the basalt pH response instead:

| parameter | "low response" member | "high response" member | reason |
|---|---|---|---|
| twelve literals | x2 (upper) | x0.5 (lower) | faster mineralisation -> more nitrification acidity and nitrate export -> lower pH |
| KNmaxd | 0.002 | 0.0064 | denitrification consumes acidity and nitrate |
| alpha | 1.3 | 6 | pH per base-saturation point |
| dissolution factor | 0.8 | 1.2 | release |
| saturation exponent | 2 (s^2) | 0.5 (sqrt s) | wetter-limited weathering |

Relative to the current code this flips only the literals' side. The
"high response" member is then the slow-nitrogen, high-alpha, fast-
weathering set: control near 6.0-6.1 (no drift, as the present lower set)
and a basalt effect of roughly +0.6, which puts the basalt plot's upper
edge at about 6.5-6.7 and overlaps the 2019-2020 observations. The "low
response" member gives a control around 5.5-5.7 and a basalt effect near
+0.15. N2O is then not aligned with the members (fast mineralisation
pairs with slow denitrification and vice versa), so its band will be
narrow, roughly 10-45 kg N/ha over five years; that is acceptable for a
band whose purpose is the pH response, and it is inside the Midwest maize
range (Hoben et al. 2011: emission factors 0.6-1.5 % of applied N; the
field's 1.6-7 kg N/ha/yr).

If Greg prefers to keep the README's "ends of the parameter ranges"
meaning, apply 2a-2d only: the band then straddles the central line for
N2O and NEE and the control pH (lower set about +0.15 above, upper set
about -0.3 below), but the basalt plot's upper edge stays near 6.1-6.2
because high alpha remains in the acid member.

## 3. Check before the 108 (the "no overshoot" guarantee)

Run the two new members in 1-D first (`TC_UC_DIMS=1d`, runs 10-12, five
years: six one-CPU jobs of about two hours), as the controllers already
allow, and accept them only if, for the control:

- N2O 2016 (unfertilised) below 5 kg N/ha and 2016-2020 within 8-50;
- peak daily N2O below 150 g N/ha/day;
- NO3 leaching below 100 kg N/ha in every year;
- NEE 2016 within +/-150 g C/m2 of the original (-258);
- bulk pH: one member above and one below the central line in every
  August; neither below 5.4 nor above 6.4; the basalt effect between
  +0.1 and +0.8 by August 2020;
- calcite 50 t/ha unchanged at the carbonate ceiling (8.1).

Only then submit the 2-D and 35-year arrays. `test_uncertainty_campaign`
(values, 12-significant-figure constants, longest result name under 255
characters) and the README table go with the member change; all of it is
inside the validation stamp, so it is a new package.

## Sources

- Hoben et al. 2011, GCB 17:1140-1152, nonlinear N2O response in Midwest
  maize: https://lter.kbs.msu.edu/citations/3130
- Heinen 2006, Geoderma 133:444-463, simplified denitrification models:
  https://research.wur.nl/en/publications/simplified-denitrification-models-overview-and-properties/
- Hou et al. 2024, doi:10.3389/fpls.2024.1482292 (first-order soil
  denitrification constants); Confalonieri et al. 2006, CropSyst
  (0.01 d-1 at 15 C), both already cited in the ensemble README.
- Kanzaki et al., SCEPTER (Appelo exchange, alpha default 3.4):
  https://gmd.copernicus.org/articles/17/4515/2024/gmd-17-4515-2024-supplement.pdf
- Zhang et al., GCB, T&C-SMEW supplement Table S2 (the literal ranges);
  Beerling et al. 2024 PNAS Fig. 1E and Dataset S06 (pH, CDR);
  Kantola et al. 2023 GCB (tower NEE).
