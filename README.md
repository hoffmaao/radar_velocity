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
  Octave in ~16 s. See `opr_vvel/README.md`.

Run against the real products: all five EAGER_2022 repeat-pass products,
both pairings, at two along-track block sizes. Coherence is good (0.95 at
50 m, 0.8 at 300 m) and the chain completes cleanly.

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
the interferogram: measure the bulk shift empirically from the data and
remove it, envelope and carrier. Empirical rather than the deterministic
inverse because how much survives depends on
`coregistration_time_shift`, which `ref_z` alone does not reveal.
Recorded in every product as `dtau_bulk` / `dtau_bulk_pred` /
`coalign_quality` / `coalign_peak_ratio` / `coalign_applied`.

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
apparent hinge sat within 0.5 km of that crossing every time. Real
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

**Measured independently by ApRES (2026-08-06).** Phase-sensitive radar was
deployed at Windless Bight in the same weeks (`~/projects/EAGER_ApRES`,
site GA04, 277 chained pairs over 5.8 days). Its per-pair displacement
profiles were fitted with THIS project's estimator - `dh = a + b*t +
c*tide`, the same model as `vdef.fitTideAdmittance` - so the two
instruments are compared with one estimator rather than through assumed
conversions. Reproduced by `scripts/diagnostics/apres_comparison.py`.

Result at 100 m: **+3.79 +/- 0.10 mm per metre of tide**, rising from
+1.1 mm at 25 m to about +4.9 mm at 145 m. Three things follow, and two of
them corrected earlier claims in this README:

- **The thin-plate flexure model on the figure was wrong** - it predicted
  -1.2 mm, opposite in sign and three times too small. The measured
  profile is now drawn on `scripts/figures/strain_rates.m` alongside it,
  because a measurement from a second instrument beats a model whose
  flexure wavelength was guessed and whose curvature goes as 1/L^2.
- **The gap is 2.2x, not 7x.** The 8.3 mm systematic floor is only just
  above the real signal, so this method is far closer to useful than the
  model-based estimate suggested.
- **ApRES resolves it 87x more finely** (1-sigma 0.095 mm against our
  8.3 mm floor), which is the scale of improvement the InSAR would need.

Caveats: only GA04 is usable (GA01 and GA05 have bed picks that drift 82
and 159 m); the ApRES grid starts at 24 m so it says nothing about the top
24 m; and no ApRES site coordinates exist anywhere, because GPS was off
for the whole deployment, so "the same site" means the same few-km area
rather than a known offset.

Advection is not the issue: measured from the per-pass GPS, the ice moves
0.3-1.2 m between passes, against 500 m blocks - so over a 3-day baseline
this Eulerian measurement samples effectively the same column an ApRES
would follow.

**What does survive**: real tidal flexure IS present, measured from the
GPS alone with no radar (`scripts/diagnostics/gps_flexure.m`). Regressing
each pass's along-track `ref_z` on its own line mean gives `a(x)` falling
monotonically by a factor ~4 along the line; GL3 and GL4 agree to 0.02
despite different reference passes and pass sets, so it is physical
rather than a GPS baseline ramp. Grounding is toward the north/northeast
end. The radar simply does not resolve the column-strain signature of it
in this dataset.

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
- `vdef.coalignPair` - measures the residual bulk fast-time shift between
  a pair by cross-spectrum group delay in a surface window and removes it,
  envelope and carrier, before the interferogram. This is the fix for the
  tide-proportional artefact; see 'Fixed: the tide-proportional artefact'
  above for the mechanism and the estimator rationale.
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
- `vdef.forwardDisplacement`, `vdef.legendreBasis` - forward model and
  basis.

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

Figures land in `figs/`. The surface-reference regression test runs the
same way, with `test_surface_reference.m` in place of the script name; the
OPR adapter's end-to-end test has its own command in
`opr_vvel/README.md`.
