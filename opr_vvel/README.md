# opr_vvel: vertical-deformation module for the OPR toolbox

Drop-in processing module for the OPR toolbox
(https://gitlab.com/openpolarradar/opr) that turns an existing
`CSARP_multipass` comp_mode 3 product into profiles of the vertical
strain rate `eps_zz` versus depth and along-track position, using the
`+vdef` package (one directory up).

## Pipeline

1. `combine_passes.m` (existing): gathers the SAR images of every pass
   over a line into one `CSARP_multipass/<pass_name>.mat`.
2. `multipass.m` comp_mode 3 (existing): coregisters every pass onto the
   baseline main - resampling in along-track, motion-compensating the
   FCS z-motion, matching fast-time axes - and saves
   `<pass_name><midfix>_multipass03.mat` with the complex images stacked
   in `data`.
3. `vvel.m` / `vvel_task.m` (this module): `vvel_task.m` is a thin OPR
   adapter (product loading, index bookkeeping, output/figure conventions)
   around the pure numerical chain in `+vdef` - `vdef.multilook` ->
   `vdef.differentialRange` -> `vdef.blockAverage` ->
   `vdef.verticalDisplacement` -> `vdef.invertStrainRate` - which reads
   its options directly off the `param.vvel` struct and is equally
   callable from standalone scripts and tests. Per pass PAIR,
   - coalignment (`vdef.coalignPair`): the residual fast-time
     misalignment between the two slices - which varies along track - is
     measured per along-track window by normalised cross-correlation of
     the trace-averaged power profiles around the surface return, and
     removed per column, envelope and carrier. On floating ice
     `multipass`'s z-motion compensation treats tidal heave as platform
     motion and misaligns the pair in proportion to the tide; this step
     undoes what survives of that shift (see the project README for the
     full mechanism and history),
   - a multilooked interferogram and coherence from the two coregistered
     SLC slices, cross product formed per pixel and only then averaged,
   - the interferogram phase unwrapped along fast time outward from the
     surface bin and converted to `dtau(twtt, x) = t_sec - t_ref`,
     referenced to zero just below the surface return,
   - coherence-weighted averaging into along-track blocks,
   - `dtau` mapped onto depth through the firn column and converted into
     the vertical displacement of each reflector relative to the surface,
     then into a relative vertical velocity using that block's own repeat
     interval,
   - a Legendre inversion for `eps_zz(d)`, reporting `c0/c1/c2`, `S1`,
     `S2`, `epszz_mean` and `p_quad` in the same schema as the earlier
     EGIG 2011-2012 vertical-strain products,
   - output `CSARP_<out_path>/<pass_name><midfix>_vvel_<ref>_<sec>.mat`
     plus three overview images (`_coh`, `_dh`, `_epszz`).

## Why this is not a cluster job

Every other OPR processing step is per frame, so one cluster task per
frame is the natural decomposition. A multipass product is instead a
single 1-1.5 GB file holding every pass in the experiment, and loading it
dominates the run. One task per pair would reload the whole stack for
every pair. `vvel.m` therefore loads it once and runs the pairs in
process, which is also how `multipass.m` itself is invoked. `vvel_task.m`
still takes one pair and still returns a success flag, so it can be
wrapped in a cluster task later if that ever becomes worthwhile.

## The two index spaces

`multipass.m` concatenates only the ENABLED passes into `data` while
leaving `pass` complete, so with any `false` in `pass_en_mask` the slice
index and the pass index diverge:

    data(:,:,k)  is  pass(pass_en_idxs(k)),
    pass_en_idxs = find(param_multipass.multipass.pass_en_mask)

`param.vvel.pairs` and `param.load.pair` are always PASS indices.
`vvel_load_multipass.m` resolves the mapping and refuses to guess if the
mask and the stack disagree; `test/test_vvel_task.m` builds its synthetic
product with a disabled pass in the middle specifically to cover this.

## Things the product does not tell you

- **The repeat interval.** `multipass` resamples the images, `ref_y`,
  `ref_z`, the layers and the attitude onto the main pass along-track axis,
  but leaves each pass's `gps_time` on its own axis. `vvel_task` maps both
  passes' `gps_time` through their own `along_track` vectors to get the
  interval per column, uses each block's own mean for that block's
  velocity, and records the along-track spread. On the EAGER 2022 lines
  the passes are walked in ~22 min and in both directions, so this is not
  a constant.
- **The phase sign.** The multipass product carries nothing equivalent to
  the polarimetric product's coregistration `row_offset`s, so `phase_sign`
  is FIXED at -1 rather than estimated: changing it asserts a different
  convention, it does not fit one. `vvel_defaults.m` defines the parameter
  and carries the full rationale.
- **Whether the surface phase was already normalised.** `surf_flatten_en`
  is a hardcoded `false` local in `multipass.m` (line ~369), not a
  parameter, so no upstream flattening has happened and
  `vdef.differentialRange` does all the surface referencing. The output
  records the flag anyway, so a future toolbox version that exposes it
  cannot silently change what these products mean.

## What surface referencing does not remove

Referencing `dtau` to zero just below the surface return removes
everything common to the whole column; "What surface referencing removes"
in the project README lists what that covers and what survives.

It does NOT remove topographic phase from a cross-track baseline, which is
not common to the column. `vvel_task` reports `baseline_y` and warns above
`param.vvel.max_baseline` (default 10 m). This matters in practice:
`EAGER_2022_GL2` contains passes tens to >100 m off the main pass track.

## Uncertainties in the product

`dtau_std_blk`, and the `dh_std_blk` / `v_std_blk` derived from it, is the
1-sigma of the block MEAN - a real error bar on `dh_blk` / `v_blk` - not
the scatter of the samples inside the block. The within-block weighted
scatter is divided by the square root of the effective number of
independent along-track samples behind the mean. Adjacent columns are not
independent, because `vdef.multilook` has already run a boxcar of
`mlook_window(2)` range lines along track, so the (Kish) weighted sample
count is deflated by `param.vvel.cols_per_look` - the along-track
counterpart of the `bins_per_look` correction `vdef.invertStrainRate`
applies in fast time.

The undeflated scatter survives in the product as `dtau_scatter_blk` and
the effective count as `neff_blk`. It is the scatter, not the error bar,
that sets the relative weights in the Legendre fit, so `S1`, `S2`,
`epszz_mean` and `eps_zz` do not depend on how many columns happened to
clear the coverage test.

`S1` and `S2` describe `eps_zz` over the normalisation depth `norm_depth`,
whereas `epszz_mean` is the mean over the range actually fitted,
`fit_top_depth` to `fit_bot_depth`. The two agree only when the fit spans
exactly `0..norm_depth`, which the production settings do not.

## A fixed artefact worth knowing about

Products produced before 2026-08-04 carried a tide-proportional strain
artefact: `multipass`'s z-motion compensation advances each pass by
`ref_z/(c/2)`, which on floating ice treats tidal heave as platform
motion and misaligns the pair slices by up to two range bins, varying
along track. The coalignment step above is the fix, and
`test/test_vvel_task.m` injects exactly this artefact - along-track-
varying per-pass delays - and asserts recovery of the full per-column
profile, not just its mean; a first, scalar generation of the fix removed
only the line mean, and its along-track residual manufactured a false
flexure hinge (see 'RETRACTED' in the project README). Two details worth
remembering: interferometric coherence never revealed the misalignment
(the returns are layered rather than pure speckle, so phase corrupted
while coherence stayed high), and the products archived as
`CSARP_vvel*_precoalign` on the server scratch preserve the pre-fix state
for comparison. Full mechanism, measurements, and remaining analysis
caveats: 'Fixed: the tide-proportional artefact' in the project README.

## Running the test

The end-to-end test builds a synthetic multipass product from a known
strain rate and checks the recovery, the index mapping, the per-block
repeat intervals and the baseline diagnostics. It runs in Octave:

```sh
docker run --rm --platform linux/amd64 -v "$PWD":/work \
  -w /work/opr_vvel/test gnuoctave/octave:latest \
  octave --no-gui test_vvel_task.m
```

Roughly 20 s, and it currently recovers S1 to 0.6% and S2 to 3.2% of
truth.

## Deployment on the CReSIS servers

- `vvel.m`, `vvel_task.m`, `vvel_defaults.m`, `vvel_load_multipass.m` ->
  `opr/matlab/processing/` (or keep on your personal path);
  `run_vvel.m` -> your `run_opr` repo (`gRadar.path_override`), edited per
  experiment.
- `+vdef` (from the project root) must be on the MATLAB path, e.g. copy to
  `opr/matlab/+vdef` or your `run_opr` repo.
- No param spreadsheet worksheet is needed or wanted: a repeat-pass
  experiment is named by its `pass_name` and spans several segments, so no
  single `day_seg` row owns it. Drive it from `run_vvel.m`.
- `server/run_vvel_scratch.m` batches the five EAGER 2022 repeat-pass
  products, reading them in place from the shared season tree through
  `param.vvel.in_fn` and writing `CSARP_vvel` (and `CSARP_vvel_seq` for
  the sequential pairing) to a user scratch tree through the `test/stubs`
  opr_* shims. See its header for the table and the `matlab -batch` launch
  line.
- `local/run_local_pair.m` runs a single pair against a product copied to
  a laptop, in MATLAB or Octave.
