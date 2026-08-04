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

### Known issue: the tidal analysis is not yet trustworthy

`scripts/figures/tidal_deformation.m` correlates the vertical column
strain against the tide, and the correlation appeared to reverse sign
along track - a flexure hinge. **Do not rely on that result yet.**

Two `multipass` builds of the SAME traverse leg, with the same passes,
give per-block strain series that are ANTI-correlated (r = -0.944), and
the difference between them is proportional to the tide (r = 0.968, 569
microstrain per metre of heave - about 57 mm of apparent displacement over
a 100 m column per metre of tide). Pass selection is not the cause:
dropping the one pass that differs between the two builds changes nothing.
The tide input itself is identical between builds up to a constant offset.

Leading hypothesis: `map.Surface` is taken from `pass.surface`, which is
the FCS surface twtt recorded BEFORE `multipass` motion-compensates each
pass's FCS z-motion by a time shift of `ref_z/(c/2)`. `ref_z` is platform
height relative to the main pass, which on a floating shelf is essentially
the tide. The surface reference bin and the coregistered data would then
be offset from each other in proportion to the tidal heave, with a sign
set by which pass is the reference - which matches every feature of the
signature.

Fix to try: derive `Surface` by tracking the surface return in the
coregistered `data` itself, per pass slice, rather than trusting
`pass.surface`. `scripts/figures/leg1_merge_check.m` is the regression
test for it: two builds of the same leg must agree.

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
matching `polarimetric.m` and the fabric chain, so `phase_sign = -1` (the
matched-filter convention settled in the fabric project). There is no
coregistration field in the multipass product to re-derive the sign from,
so it is a fixed parameter here rather than auto-detected.

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

Surface flattening is NOT applied upstream: `surf_flatten_en` is a
hardcoded `false` local in `multipass.m` (~line 369), not a parameter, so
`vdef.differentialRange` does all the surface referencing. The output
records the flag anyway, so a future toolbox version that exposes it
cannot silently change what these products mean.

## Running the validation

```sh
docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work/scripts \
  gnuoctave/octave:latest octave --no-gui synthetic_vertical_velocity.m
```

Figures land in `figs/`. The surface-reference regression test runs the
same way, with `test_surface_reference.m` in place of the script name; the
OPR adapter's end-to-end test has its own command in
`opr_vvel/README.md`.
