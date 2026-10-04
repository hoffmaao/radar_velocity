# radar_velocity

Vertical deformation from repeat-pass radar interferometry, built on the
OPR `+multipass` chain. Target application: the 2022_Antarctica_Ground
accum3 repeat passes at Windless Bight, McMurdo.

Structured to mirror `../fabric_anisotropy`: a dependency-free numerics
package (`+vdef`, the analogue of `+ptt`) plus a drop-in OPR processing
module (`opr_vvel/`, the analogue of `opr_fabric/`) that adapts it to the
toolbox's file and product conventions.

## Status

Working and validated:

- `+vdef` - the full numerical chain, end-to-end verified in Octave.
- `scripts/synthetic_vertical_velocity.m` - synthesises a repeat-pass SLC
  pair from a known strain rate, runs the whole chain, and asserts
  recovery. Currently recovers S1 to 0.3% and S2 to 1.0% of truth.
- `scripts/test_surface_reference.m` - regression test for the surface
  reference bin: a trace whose `Surface` is NaN, or whose reference twtt
  falls outside the fast-time axis, must be dropped rather than clamped to
  the first or last bin. Real `pass.surface` carries NaN wherever the
  surface tracker failed, and a clamped trace comes back fully finite and
  entirely wrong.
- `opr_vvel/` - the OPR adapter: `vvel.m`, `vvel_task.m`,
  `vvel_defaults.m`, `vvel_load_multipass.m`, the local/server drivers,
  and `test/test_vvel_task.m`, which builds a synthetic multipass product
  (with a disabled pass in the middle, to cover the pass-index vs
  data-slice mapping) and recovers S1 to 0.6% and S2 to 3.2% of truth in
  Octave in ~20 s. See `opr_vvel/README.md`.

Run against the real products: all five EAGER_2022 repeat-pass products,
in the main, sequential and all-pairs pairings, at two along-track block
sizes. Coherence is good (0.95 at 50 m, 0.8 at 300 m) and the chain
completes cleanly. The whole server chain is reproducible from one entry
point; see `opr_vvel/server/README.md`.

Not yet done:

- `docker-compose.yml` - the MATLAB container.

### Fixed: the tide-proportional artefact (2026-08-04)

The tidal analysis originally carried a tide-proportional artefact: two
`multipass` builds of the SAME traverse leg gave anti-correlated strain
(r = -0.944), with a difference proportional to the tide (about 57 mm of
apparent column displacement per metre of heave).

**Mechanism, confirmed from source and from the data.** `multipass`
comp_mode 3 motion-compensates each pass's FCS z-motion by advancing it
`ref_z/(c/2)` in fast time (multipass.m:512-522, envelope and phase), and
`pass.surface` is never updated for the shift (`layers.twtt_ref`
explicitly is, line 667-669). On grounded ice `ref_z` is platform motion
and the compensation is right. On the floating shelf the platform and the
surface ride the tide together, so `ref_z` is essentially the tide and
the compensation displaces the returns by an amount that never was a
range change.

**Whether it reaches a given product depends on that product's
calibration (2026-08-05).** `multipass`'s own
`param.multipass.coregistration_time_shift` applies a per-pass fast-time
shift derived from the comp_mode 2 coregistration stage, and it removes
this misalignment as a side effect. Across the five EAGER products the
field splits them cleanly in two:

| product | `coregistration_time_shift` | surface offset measured in `data` |
|---|---|---|
| `EAGER_2022` | all zeros | full residual, 0.97-1.07 x `-(ref_z_sec - ref_z_ref)/(c/2)` |
| `EAGER_2022_GL4` | all zeros | full residual |
| `EAGER_2022_GL1` | nonzero, up to 2 bins | aligned, within +/-0.9 ns of zero |
| `EAGER_2022_GL2` | nonzero | aligned |
| `EAGER_2022_GL3` | nonzero | aligned |

So this correction is a **repair for uncoregistered products**, and on a
properly coregistered one it must measure ~0 and do nothing. That split
is not geography: it is exactly the split seen in the tidal response,
where GL1/GL2/GL3 show a hinge and `EAGER_2022`/GL4 do not.

**The fix** is `vdef.coalignPair`, applied per pair in `vvel_task` before
the interferogram: measure the shift empirically from the data - per
along-track window, since it varies along track (see the retraction
below) - and remove it, envelope and carrier. Empirical rather than the
deterministic inverse because how much survives depends on
`coregistration_time_shift`, which `ref_z` alone does not reveal.
Recorded in every product as `dtau_bulk` / `dtau_bulk_pred` /
`coalign_quality` / `coalign_peak_ratio` / `coalign_applied`, plus - since
coalignment went per-column - the applied profile and per-window
diagnostics (`dtau_bulk_profile`, `dtau_bulk_win`, `coalign_x_win`,
`coalign_quality_win`, `coalign_n_win`, `coalign_n_win_ok`).

**The estimator, replaced 2026-08-05.** It is now the normalised
cross-correlation of the two slices' trace-averaged POWER profiles in a
surface window, FFT-upsampled 32x. The previous cross-spectrum group
delay was validated only on a synthetic that contained no surface return
at all - depth-decaying white noise, whose trace-averaged profile has no
bin-scale structure - and on the real products it ran 1.2-1.3x the true
residual on the uncoregistered products and returned up to 8.5 ns of pure
noise on the coregistered ones, where the truth is ~0, all at quality
0.965-1.000. Since coalignment is applied to every pair, that noise was
corrupting the *correctly* calibrated products. The envelope correlation
reproduces the truth at ratio 0.97-1.07 on the former and stays within
+/-0.9 ns of zero on the latter. The synthetic now carries a realistic
band-limited surface return, without which it cannot exercise the
estimator that runs on real data.

The `coalign_max_lag` bound is load-bearing rather than a formality: the
trace-averaged surface profile carries a range sidelobe at +/-5.4 bins at
0.37-0.55 of the main peak, and on one pair (GL1 / 20221211_07) that
sidelobe outranked the true peak and gave -19 ns against a true 0.8 ns. A
peak-dominance test does not separate these cases - that pair's
peak-to-sidelobe ratio was 1.8 while a correctly measured pair sat at
1.03 - so the bound is physical (3 bins, above the tidal range over c/2
and well inside the sidelobe), and a peak found on the bound is rejected
rather than clamped.

**What the leg-1 comparison says now**
(`scripts/figures/leg1_merge_check.m`): the two builds' tidal-response
profiles went from uncorrelated (0.125) to correlated in shape (0.807),
but an offset in absolute r remains. That residual is an ANALYSIS
limitation, not a processing one: the two builds reference different
epochs (2022-12-09 vs 2022-12-12), and a plain correlation is invariant
to the constant strain offset but not to the secular trend, which
aliases in through the sample covariance of time with tide - a quantity
that depends on which pairs a build happens to contain. The
reference-invariant quantity is the tide admittance of the joint
strain = a + b*t + c*tide fit, `vdef.fitTideAdmittance`: re-referencing
shifts strain, t and tide by constants that the intercept absorbs, so
the trend b and admittance c cannot move. `scripts/test_tide_admittance.m`
demonstrates the aliasing deterministically (two noiseless builds of the
same truth differ in plain r by up to 0.23 while their admittances are
bit-equal) and is the unit test. Both analysis scripts now read the
admittance and its trend-removed partial correlation; plain r survives
only in printed tables for continuity.

### RETRACTED: the tidal flexure hinge was the artefact (2026-08-05)

The project's headline result - vertical column strain correlating with
the tide and reversing sign along track at a flexure hinge - **does not
survive correct coalignment**. It was the residual misalignment.

**Why a scalar correction manufactured it.** Both
`coregistration_time_shift` and the first version of `vdef.coalignPair`
remove a line MEAN. What survives is proportional to `(a(x) - 1)`, where
`a(x)` is the local surface tidal admittance normalised so the line mean
is 1. That residual changes sign exactly where `a(x) = 1` - a position
set by the arbitrary normalisation of the survey, not by the ice.
Measured: `a(x)` crosses 1 at 2.5-3.0 km in every product, and the
apparent hinge sat within ~0.5 km of that crossing in four of the five
products; the exception is `EAGER_2022` at 1.50 km, which is also the
only wholly uncalibrated product, so its residual is the largest and
least well described by the simple `(a(x) - 1)` form. Real
flexure predicts something different - bending strain follows the
CURVATURE of the deflection, so it changes sign at an inflection of
`a(x)`, and `a(x)` is concave-down at every point in the surveyed window.
The data matched the artefact prediction and contradicted the flexure
one.

**What per-column coalignment did to it.** The correlation between the
measured admittance and the artefact predictor
(`scripts/diagnostics/artefact_vs_signal.m`) went from -0.84…-0.97 with
every product at p < 0.005, to +0.34 / -0.63 / +0.29 / +0.55 / -0.64 with
none significant. The change-point detector now finds NO hinge in any of
the five products, where with the scalar correction all five showed one
at contrast 0.96-1.60. The surviving admittance is -33 to +40 µε/m with
no along-track structure.

The four products that had "agreed" on the hinge to ~570 m agreed because
they share a survey geometry and therefore the same `a(x) = 1` crossing -
not because they each saw the same ice.

**The strain rates themselves** are in
`scripts/figures/strain_rates.m`, which reports the secular term and the
tide admittance of the joint fit over a systematic floor measured from
the two builds of leg 1 rather than assumed from a noise model. Every
product's line mean sits INSIDE that floor for both quantities, at every
depth. The floor scales as 1/z, from 148 µε/m and 4.3e-2 /yr at 50 m to
29 µε/m and 8.8e-3 /yr at 250 m.

**Measured independently by ApRES (2026-08-06; corrected 2026-08-17;
position-matched 2026-08-19).** Phase-sensitive radar was deployed at
Windless Bight in the same weeks (`~/projects/EAGER_ApRES`; four sites
carry pair products, and the site of record GA04 has 277 half-hour pairs
over 5.8 days).

The first comparison (`scripts/diagnostics/apres_comparison.py`) chained
the per-pair displacement profiles end to end and fitted the cumulative
series with this project's estimator. It gave +3.79 mm per metre of tide
at 100 m. **That number is RETRACTED as an analysis artefact.** The
retraction rests on an internal inconsistency rather than on a vote among
methods: the chained analysis and the rate-method check below use the
SAME instrument at the SAME site over the SAME window and give
incompatible answers, which a bug explains and physics does not. The
chaining method also has specific known flaws already recorded as leads -
it integrates per-pair displacements as a random walk, applies iid OLS to
errors that are autocorrelated by construction, and may double-count.
Neither of those depends on any radar number, so nothing below changes
the retraction. The bug has not yet been localised - candidates are sign
handling in the chain, cumulative drift correlated with the tide, and the
shallow grid points the vsr fit excludes - so the chained script is kept,
marked RETRACTED and with its CSV export removed, for the pending bug
hunt.

The comparison of record is `scripts/diagnostics/apres_rate_check.py`,
which needs almost no processing of ours: if strain = A*tide then the
strain RATE of a half-hour pair is A*d(tide)/dt, so regressing the ApRES
project's own `vsr_per_year_medfilt` against the CATS2008 tide rate gives
A directly - no chaining, no cumulative drift, no depth-grid handling.

The numbers, over the top 100 m in mm per metre of tide:

- ApRES rate method (GA04): **-1.24 +/- 0.04** (R2 = 0.75)
- radar network, pooled over the four calibrated lines: **-1.46 +/- 0.55**
- thin-plate flexure model: about **-1.2**
- radar network AT GA04's OWN along-track position: **+2.22 +/- 0.86**
  (500 m bins) or **+0.30 +/- 0.43** (1 km bins)

**The ApRES site positions were recovered on 2026-08-19** and are in
`data/gis/eastwind_2022_2023_apres_xy.txt` (EPSG:3031 metres). That file
lists exactly the four sites that have pair products - GA01, GA04, GA05,
GA10 - so it also settles which of the twelve GA sites in
`data/gis/McM_GNSS_ApRES.txt` are ApRES rather than GNSS. GA04 projects
to 4.70 km along track and 0.07 km off the line, i.e. into the
grounding-line approach rather than the flat floating section the earlier
write-up assumed.

That makes the first three numbers above a comparison of DIFFERENT
quantities. The radar profile runs -3.7 to +2.5 mm/m along the line, so
its line mean is not what a point measurement should be compared with,
and the close -1.24 against -1.46 was an averaging coincidence. Like for
like at 4.70 km the two are opposite in sign and about 4 sigma apart. The
pooled line mean is not retracted - it is a different quantity, and it
still carries the strain-rate and error-budget figures.

**An unexplained coincidence, stated rather than buried.** Locating the
site also cost the retraction above the sign argument it used to lean on.
Against the LINE MEAN the retracted +3.79 was the odd one out of four;
position-matched it is 2-2 - rate method (-1.24) and flexure model (-1.2)
negative, chained (+3.79) and radar at GA04 (+2.22) positive. So the
retracted value and the radar at the site share a sign. That is not
currently understood, and it does NOT rehabilitate the chained profile:
the radar's far-end MAGNITUDE is itself unresolved (+2.22 +/- 0.86 at
500 m bins against +0.30 +/- 0.43 in the 1 km bin containing the site,
per `scripts/diagnostics/far_end_check.m`, below), and the chaining flaws
above are methodological rather than inferred from any disagreement. It
is an open question worth revisiting if the chaining bug is ever
localised.

How much the far-end radar bin can carry is tested before it is drawn, in
`scripts/diagnostics/far_end_check.m`, and the verdict is deliberately
split. The SIGN is robust: all four lines are independently positive at
4.65-4.75 km, a jackknife dropping each line in turn leaves +1.94 to
+2.64, and the radar-free GPS-curvature prediction also rises toward the
grounding line. The MAGNITUDE is not resolved: +2.22 +/- 0.86 at 500 m
bins against +0.30 +/- 0.43 at 1 km, because the adjacent 4.25 km blocks
are negative, and the far blocks are thinner than mid-line ones (valid
samples -22%, fitted depth span -21%). The disagreement with ApRES
survives either binning (4.0 and 3.6 sigma), so it is not a bin artefact.
`scripts/figures/convergence.py` therefore draws both stack rows and
compares every estimate at the ApRES site rather than pooling.

Caveats: only GA04 is usable - GA01 and GA05 have bed picks that drift 82
and 159 m, and GA05 and GA10 have only a handful of pairs each.

Advection is not the issue: measured from the per-pass GPS, the ice moves
0.3-1.2 m between passes, against 500 m blocks - so over a 3-day baseline
this Eulerian measurement samples effectively the same column an ApRES
would follow.

**The evidence, in one figure**: `scripts/figures/tidal_evidence.m` draws
the whole argument - the apparent hinge under a scalar coalignment, its
disappearance under per-column, the GPS a(x) profile showing the scalar
residual must flip sign where a(x) = 1, and the collapse of the
correlation with the artefact predictor. Note the scale in panel (a): the
apparent signal reached 40 mm per metre of tide against an 8.3 mm floor,
which is why it was convincing.

**What does survive**: real tidal flexure IS present, measured from the
GPS alone with no radar (`scripts/diagnostics/gps_flexure.m`). Regressing
each pass's along-track `ref_z` on the tide gives `a(x)` falling
monotonically by a factor ~6 along the line; GL3 and GL4 agree to 0.02
despite different reference passes and pass sets, so it is physical
rather than a GPS baseline ramp. Grounding is toward the north/northeast
end. The radar simply does not resolve the column-strain signature of it
in this dataset. (The first estimator regressed on each pass's own line
mean; see 'The tide is CATS2008' under the elasticity section for why
that is biased and what replaced it.)

**Interpretation caveats that remain**: the five products are four legs
of ONE line (120-145 m apart, walked ~24 min apart), not independent
lines, so they were never five independent tests of anything.

## Target data

The 2022_Antarctica_Ground season contains a genuine repeat-pass
experiment at Windless Bight, McMurdo - and **the multipass comp_mode 3
products already exist on the CReSIS servers**, so the `combine_passes` /
`multipass` stages do not need re-running. Read from
`mem1.cresis.ku.edu:/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass`
on 2026-08-03:

| product | passes | main pass | line | dates |
|---|---|---|---|---|
| `EAGER_2022` | 13 | 1 (20221209_03) | 4.9 km | 12-09 to 12-12 |
| `EAGER_2022_GL1` | 14 | - | ~4.8 km | 12-09 to 12-12 |
| `EAGER_2022_GL2` | 15 | 5 (20221210_05) | 4.91 km | 12-07 to 12-12 |
| `EAGER_2022_GL3` | 13 | 11 (20221211_09) | 4.82 km | 12-09 to 12-12 |
| `EAGER_2022_GL4` | 14 | 9 (20221211_02) | 4.79 km | 12-09 to 12-12 |

Each is a `<pass_name>_multipass03.mat` of 1.0-1.5 GB holding
`data` as 4250 x ~1930 x Npass single complex.

Geometry, from the products themselves: accum3 at **fc = 750 MHz**,
**fs = 300 MHz** (dt = 3.333 ns), 4250 fast-time bins spanning -2.070 to
12.093 us with the surface at ~0.002 us, so the record reaches
**~1020 m** below the surface. Along-track sampling is **2.5 m** over a
~4.8 km line (~1930 range lines), which is why `block_size` is set to 200
here rather than the OPR per-frame default of 1000.

The passes run 2022-12-09 to 2022-12-12 at 21-25 min each, with
sub-hourly to few-hourly sampling on 12-10 and 12-11 and gaps out to
~3 days, so the tidal band (M2 = 12.42 h, K1 = 23.93 h) is resolvable -
relevant because Windless Bight is floating and tidal flexure strains the
column. Passes are walked in both directions depending on the product,
which is why the repeat interval is reconstructed per along-track column
rather than taken as one number per pair.

**Cross-track baselines** decide which products are usable. GL3 and GL4
are clean: every pass sits within ~4 m of the main_pass, most within 1 m.
GL2 is not - 20221210_06/07/08 reach 45 m, 146 m and 22 m, far enough
that topographic phase (which surface referencing does not remove) will
contaminate those pairs. Start with GL3 and GL4.

All the underlying segments have `CSARP_sar` frames, so the products can
be rebuilt if a different pass selection or main pass is wanted.

## Method (`+vdef`)

Depth `d` is measured DOWNWARD from the snow surface (the opposite of
`+ptt`, which uses height above the bed). SI units throughout; only the
reported strain rates and velocities are converted to per-year.

- `vdef.constants`, `vdef.defaultParams` - physical constants and the
  firn column parameterisation.
- `vdef.firnColumn` - single-stage Herron-Langway density recast in depth
  and pinned at `rho_sfc` / `rho_bco`, plus refractive index (Kovacs or
  Looyenga) and the vertical twtt table.
- `vdef.depthFromTwtt` - twtt below the surface to depth and local `n`.
- `vdef.coalignPair` - measures the residual fast-time misalignment of a
  pair, per along-track window, by normalised cross-correlation of the two
  slices' trace-averaged power profiles in a surface window, and removes
  it per column, envelope and carrier, before the interferogram. This is
  the fix for the tide-proportional artefact; see 'Fixed: the
  tide-proportional artefact' above for the mechanism and the estimator
  rationale, and the retraction section for why it must be per column.
- `vdef.multilook` - boxcar interferogram and coherence from a coregistered
  SLC pair. Cross product per pixel, averaged after - never the reverse.
- `vdef.differentialRange` - interferogram phase to `dtau(twtt, x)`,
  referenced to zero just below the surface return. Unwraps along fast
  time outward from the surface bin rather than with a 2-D unwrapper,
  because `dtau` is smooth in depth. Bridges short incoherent gaps and
  invalidates everything below a long one, so it fails toward "no data"
  instead of accumulating confident nonsense.
- `vdef.blockAverage` - coherence-weighted along-track block averaging.
  This is where the precision comes from: the per-pixel signal over a
  few-day baseline is a few hundredths of a radian, well under the
  per-pixel phase noise.
- `vdef.verticalDisplacement` - `dtau` to vertical displacement relative
  to the surface, `dh = dtau*c/(2*n(d))`, using the LOCAL index at the
  reflector depth. Optional densification correction, off by default
  (negligible at days, not at seasons).
- `vdef.invertStrainRate` - Legendre inversion for the strain-rate
  profile, reporting `c0/c1/c2`, `S1`, `S2`, `epszz_mean` and `p_quad` in
  the same schema as the earlier EGIG 2011-2012 vertical-strain products
  so the two are directly comparable.
- `vdef.fitTideAdmittance` - the reference-invariant joint fit
  strain = a + b*t + c*tide per along-track block. The trend `b` and tide
  admittance `c` cannot depend on the reference-pass choice, unlike a
  plain correlation with tide; this is the acceptance metric of the tidal
  analysis. Unit test: `scripts/test_tide_admittance.m`.
- `vdef.invertNetwork` - per-epoch displacement from every pass pair at
  once, under a sum(x) = 0 datum (no privileged reference epoch) with
  robust rejection of inconsistent pairs; the fit residuals are the
  closure errors. Unit test: `scripts/test_invert_network.m`.
- `vdef.surfaceAdmittance` - the surface tidal admittance `a(x)` from
  repeat-pass platform heights: per-block regression on an EXTERNAL tide,
  with the line-mean residual as a pass gate and a common-mode nuisance
  column. The line mean itself is not a usable regressor (see the
  elasticity section). Unit test: `scripts/test_surface_admittance.m`.
- `vdef.surveyLines` - the four independent profiles, in one place.
  The multipass directory holds FIVE products and the fifth, `EAGER_2022`,
  is not a fifth line: it is a second build of GL1. The two centre lines
  are a median 0.5 m apart where the true neighbours are 140, 267 and
  402 m away, the column counts and end points are identical, and GL1's
  pass list is a strict superset, its thirteen day segments plus
  20221209_01. Every table, map and pooled statistic that listed all five
  therefore counted leg 1 twice; `scripts/figures/leg1_merge_check.m` said
  so in its header long ago and the product lists were never updated.
  They are now, and the duplicate is kept only where two builds of the
  same ice are the point: the systematic-error floor in
  `scripts/figures/strain_rates.m` and `scripts/figures/tidal_evidence.m`,
  the between-build matrix, and the coalignment regression test. Those
  call `[lines, dup] = vdef.surveyLines()` and say why. Worth knowing
  before quoting either build: at 100 m GL1 comes out consistent with
  zero while its second build comes out at -2.9 to -3.9 mm per m of tide
  at every block length, on the same track, differing by one pass and by
  which pass is the reference.

- `vdef.complexBlocks` and `vdef.tidalStack` - the COHERENT ALL-PAIRS
  tidal estimator. Every pair's multilooked interferogram is referenced
  to the surface bin by phase and block-averaged as a complex field,
  never unwrapped; the tidal column response a(z, x) is then the trial
  response at which all pairs, counter-rotated by the model phase
  `phase_sign * 4 pi fc n(z)/c * a * dtide`, add coherently. Every pair
  constrains every other through the shared model, the response comes
  out at every depth, and the estimate is independent of the
  unwrap -> dh -> regression chain. Its peak coherence `F` says how
  tide-locked the phase is (0.94-1.00 on these lines: the residual
  mis-registration, seen directly); its 2-D scan takes each pair's
  MEASURED residual misalignment as a second regressor and reports the
  collinearity `rcol` per block, because with `|rcol|` near 1 the two
  are one regressor and the controlled estimate is a ridge, not a
  number. A self-test injects a known response as a phase rotation and
  requires it back as an exact shift, which is what pins the phase
  sign. Driver: `scripts/diagnostics/tidal_stack.m`, figure:
  `scripts/figures/tidal_stack_figure.m`. Unit test:
  `scripts/test_tidal_stack.m`.
- `vdef.beamFlexure` - the grounding zone as an Euler-Bernoulli beam of
  varying thickness on a hydrostatic foundation,
  `d2/dx2[D(x) w''] = rho_w g (A0 - w)` after Holdsworth (1969), solved by
  central differences with the beam clamped landward and flat at the
  far-field seaward boundary. Nondimensionalised by the flexural length and
  row-equilibrated, because a biharmonic operator on a few hundred nodes
  otherwise loses most of its digits.
- `vdef.invertElasticModulus` - least-squares inversion of an observed
  flexure profile for the effective Young's modulus `E*`, after Elgart,
  Minchew and Meyer (2025). See 'Effective elasticity' below. Unit test:
  `scripts/test_beam_flexure.m`.
- `vdef.forwardDisplacement`, `vdef.legendreBasis` - forward model and
  basis.

### Effective elasticity from tidal flexure

`vdef.invertElasticModulus` fits the beam above to a measured profile of
tidal surface deflection and returns the single effective elastic
parameter `E*` that minimises the misfit. The driver is
`scripts/diagnostics/elastic_modulus.m`, which inverts the GPS surface
tidal admittance `a(x)` of each calibrated leg against the tracked-bed
thickness `h(x)` - the OPR layer_tracker bed run over the four main-pass
segments, validated to ~2 m against the ApRES bed at GA10 - with the
BedMachine pseudo-layer and two constant ApRES/radar thicknesses kept as
bracket cases. It prints and returns fits;
`scripts/figures/flexure_inversion.m` calls it and draws them, so the
inversion has one implementation and the plot another and neither can
drift. Its four panels are the observation and the fit, the misfit surface
over `E*` and clamp position, the misfit profiled over clamp position in
units of the data variance, and `E*` against assumed thickness.
`scripts/figures/elasticity_results.m` is the results companion: the
fitted beams, the tracked-bed `h(x)` against BedMachine and the ApRES bed
depths, and a forest of `E*` against the published estimates (Vaughan
1995; Sayag and Worster 2013; Elgart and others 2025; laboratory ice).

**The tide is CATS2008, not the line mean.** `a(x)` is a regression of
block heights on the tide across passes, and the regressor matters.
The first estimator used each pass's own line-mean height. Against the
CATS2008 prediction at each pass mid-time
(`scripts/diagnostics/cats2008_tide.m`) the line means carry a 13-15 cm
rms non-tidal residual, about half of it a height offset UNIFORM along
the line - a per-pass platform/GPS error. A uniform offset enters every
block and the line-mean regressor with coefficient one, so the slope is
pulled toward one and the profile compressed toward flat, which the beam
fit reads as a stiffer beam: on synthetic data with this geometry 14 cm
rms of such offsets biased `E*` by +58%. `vdef.surfaceAdmittance` now
regresses on the external tide, uses the line-mean residual as a PASS
GATE (GL2's first pass, 20221207_03, sat 1.13 m off the tide and alone
produced its negative-admittance block and its 16 GPa) and as a
common-mode nuisance column, and the same tide and gate are handed to
the strain chain (`scripts/test_surface_admittance.m`). One helper,
`scripts/diagnostics/pass_tide.m`, is where every script that regresses
on the tide gets it - the flexure driver, the strain loader, the tidal
response map, the strain-rate profiles, the line sections and the
per-line tidal deformation figures - so all of them use the same tide
and drop the same passes. `scripts/figures/strain_flexure.m` draws the
englacial strain admittance per line against the joint fit's beam, the
companion of the surface panel in `flexure_inversion.m`.
`scripts/diagnostics/tidal_depth_profile.m` takes the englacial response
to a FUNCTION OF DEPTH: the project's Legendre expansion of a column
profile (`vdef.invertStrainRate`) at order 3 on every pair's stored
dh(z), the coefficients network-inverted and tide-fitted one by one, and
the thin-plate shape `A (z_n z - z^2/2)` fitted to them in coefficient
space to read the neutral-plane depth off the data. The plate shape is
exactly three Legendre terms - `m_2 = -A H^2/12` alone carries the
amplitude, `z_n = H/2 - (H/6)(adm_1/adm_2)` is a ratio of two fitted
numbers, and every order above 2 is identically zero for any plate - so
the neutral plane is measurable BLOCK BY BLOCK and the order-3 term is a
shape test that refitting cannot absorb. Fit quality is judged against a
basis-free run of the same chain, with dh interpolated onto a depth grid
and every depth network-inverted and tide-fitted on its own. On two legs the expansion
earns its order, the residual against the pointwise profile falling to
order 3 and beyond, while on the other two it degrades past order 1.

**The bending amplitudes it reports are not evidence of flexure**, and
the test that shows this is the project's own systematic floor: run the
same chain on `EAGER_2022` and `EAGER_2022_GL1`, which are the same leg
with the first uncalibrated, and the UNCALIBRATED build wins every
quality criterion - a column response of 6.4 mm/m against 1.7, a bending
amplitude of 15 sigma against 6.5, and a neutral plane marching
monotonically along the line at 148-176 m. Coherence, plate shape, a
tight `z_n` and a clean order sweep are reproduced in full by a build
known to carry the coalignment artefact, and the difference between the
two builds exceeds the calibrated build's entire signal. The formal
errors do not see this. Nothing from this driver should be read as
flexure until it reproduces between independent builds of the same ice. What no
regressor removes is the offsets' random projection onto the tide,
which shifts a whole profile by a constant; the driver therefore prints
a LEAVE-ONE-PASS-OUT JACKKNIFE of `E*` beside the formal interval, and
that is the error to quote. The three legs are the same passes walked
minutes apart, so their agreeing with each other does not test it.

The primary fit is JOINT with the englacial strain admittance - the
radar's own dh(100 m) per metre of tide, loaded by
`scripts/diagnostics/load_strain_admittance.m` through the same network
inversion and reference-invariant fit as the tidal response map. The beam
predicts it as `amp * nu/(1-nu) * (z_n*zr - zr^2/2) * w''(x)` with the
SAME shared amplitude as the deflection fit, and that coupling is the
point: the line-mean normalised `a(x)` constrains only the beam's shape,
while the strain admittance is absolute and its curvature scale goes as
`D^(-1/2)` - an amplitude equation the surface expression cannot supply,
which is what closes the `E*`-clamp trade-off on a window that never
sees the far field (asserted in `test_beam_flexure.m`, check 11). Three
things the joint fit is honest about by construction: it REFUSES to run
without shape sigmas, because with unweighted shape rows the shared
amplitude is set by the strain and the constraint silently cancels; it
reports per-dataset reduced chi-squareds against the INPUT sigmas, so
englacial strain beyond thin-plate bending shows up as `X2s` above 1
rather than vanishing into a pooled variance; and with `opts.rescale` it
re-weights the two datasets from their own residuals, because at face
value the strain rows misfit at 2-4 against 0.1-1 for the shape and were
carrying twenty-odd times the weight their scatter supports (on the real
lines the re-weighting leaves `E*` unchanged and brings the joint
interval down to the shape one - at its real precision the strain adds
little here). The interval itself is never taken narrower than the input
sigmas imply: residuals smaller than the sigmas are, on these lines, the
signature of errors shared by every block, and scaling the interval by
them made it two to three times narrower than the pass jackknife
(`test_beam_flexure.m`, checks 14 and 15).

LOCAL `E*(x)` is available through the inverter's `E_patch` mode (the
searched modulus applies only inside a window, via the exact
`D = E h^3` equivalence) and drawn by `scripts/figures/elasticity_map.m`,
which runs a constant-truth control through the identical pipeline
beside the data. Flexure is nonlocal, so nothing is resolved below the
~1.2 km flexural length, and the constraint dies seaward where moment
and curvature both vanish - on this line the map is meaningful over
roughly the first two kilometres and honestly unconstrained beyond,
which the per-patch intervals and the control both show
(`test_beam_flexure.m`, check 12).

Three parameters, two of them searched. `E*` and the landward boundary
`x0` are found on a grid and then by a pattern search - the same role
MATLAB's `patternsearch` plays in Elgart and others, written out so that
it needs no toolbox and behaves identically in Octave. The far-field
amplitude enters linearly and is eliminated in closed form at every trial,
which is what lets a profile with an arbitrary overall scale be inverted:
only the SHAPE is being fitted, so `a(x)` may be per metre of tide, per
metre of line mean, or in any other units. The pattern search refines off the lattice but stays INSIDE the
grid it started from, so an ill-conditioned patch ends at the nearest edge
instead of walking off into moduli the grid never proposed; `R.interior`
then reports that edge - judged from the refined optimum, not from the
grid profile alone - and `elasticity_map` draws such a patch as unresolved
rather than as its boundary value (`test_beam_flexure.m`, check 13).

Four things about the result that are not caveats but part of it:

- **`E*` is only ever `E*h^3`.** Flexure constrains the rigidity, so a 10%
  thickness error is a 33% modulus error. The driver runs three thickness
  cases and reports the achieved exponent, which comes out at -3.00.
- **It is not the Young's modulus of ice.** Surface flexure does not
  separate bending of the ice from bending of the bed beneath it, or from
  the anelastic part of the response at tidal frequency. Elgart and others
  find 0.6-9 GPa across three Ross sites, against ~9 GPa in the laboratory.
- **This survey is a partial window.** The line spans ~4.8 km inside the
  flexure zone and reaches neither flat end, so `x0` sits outside the data
  and trades off against `E*`. On synthetic data with that exact geometry
  the fit lands tens of percent off, with an interval seven times wider
  than the same method achieves on a full profile. A TIGHT interval on this
  geometry is the thing to distrust.
- **Block means are modelled as block means.** `a(x)` is a mean over a
  500 m block and the deflection is curved on the scale of a flexural
  length, so the block mean is not the value at the block centre. Left
  uncorrected that bias runs an order of magnitude above the formal error
  on `a(x)`, and it does not announce itself as a bad fit. `opts.avg_width`
  is what turns it off.

### Sign conventions

`dh < 0` means the column between the surface and the reflector has
SHORTENED (vertical compression). The interferogram is `sec .* conj(ref)`,
matching `polarimetric.m` and the fabric chain, so `phase_sign` is fixed
at -1 rather than auto-detected. `opr_vvel/vvel_defaults.m` defines the
parameter and carries the full rationale.

### What surface referencing removes

Everything common to the whole column: antenna height change, rigid-body
tidal heave, bulk timing drift, and the unwrapping constant. What survives
is the differential vertical motion within the column - the vertical
strain signal - plus tidal flexure, which at Windless Bight is real
signal rather than error.

## Input product

`vvel_task` consumes the `+multipass` comp_mode 3 output,
`CSARP_multipass/<pass_name><midfix>_multipass03.mat`, which holds:

- `data` - Nt x Nx x Npass complex SLCs, all coregistered onto
  `pass(baseline_main_idx).time` and the main pass along-track axis
- `pass` - per-pass geometry (`ref_z`, `ref_y`, `along_track`, `layers`,
  `gps_time`, `wfs`, `surface`); `data`/`ref_data` are stripped at save
- `ref` - the baseline main pass
- `param_combine_passes`, `param_multipass`

`data` holds only the ENABLED passes while `pass` stays complete, so with
any `false` in `pass_en_mask` the slice index and the pass index diverge:
`data(:,:,k)` is `pass(pass_en_idxs(k))`. `vvel_load_multipass` resolves
that mapping; everything user-facing is a pass index.

`ref` is a copy of the main pass taken before the coregistration loop
and can still carry a full SLC image, so `opr_vvel` deliberately does not
load it - the main pass's own vectors in `pass` are the common axes.

Surface flattening is NOT applied upstream, so `vdef.differentialRange`
does all the surface referencing. See "Things the product does not tell
you" in `opr_vvel/README.md` for why, and for what recording the flag in
the output protects against.

## Running the validation

```sh
docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work/scripts \
  gnuoctave/octave:latest octave --no-gui synthetic_vertical_velocity.m
```

Figures land in `figs/`. EVERY figure and every fitted product goes
there, whichever script drew it and whichever machine it ran on:
`vdef.figureDir` derives the path from the position of the `+vdef`
package, so the outputs land beside the code that made them, in the
working copy locally and in the deployed copy on the processing server.
The directory is gitignored. Before this there were two absolute scratch
paths on the server hardcoded into a dozen scripts, and a figure's home
depended on which script drew it. Set `RADAR_VELOCITY_FIGS` to redirect
a run without editing anything. The fitted `.mat` products
(`flexure_fit_cats.mat`, `tidal_stack.mat`, `tidal_depth_profile.mat`)
live there too, because they are the exact inputs the figures were drawn
from and several figure scripts reload them to redraw without refitting.

The surface-reference, tide-admittance,
network-inversion, surface-admittance, tidal-stack and beam-flexure tests
run the same way, with `test_surface_reference.m`, `test_tide_admittance.m`,
`test_invert_network.m`, `test_surface_admittance.m`, `test_tidal_stack.m`
or `test_beam_flexure.m` in place of the script name; the OPR adapter's
end-to-end test has its own command in `opr_vvel/README.md`. CI
(`.github/workflows/tests.yml`) runs all eight
entrypoints in Octave on every push.
