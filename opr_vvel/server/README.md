# Server workflow: reproducing every product from the archive

Everything below runs on mem1 (`ssh hoffmana_sta@mem1.cresis.ku.edu`,
VPN required). The goal of this file is that any product in the chain can
be regenerated from the archived inputs and PROVEN identical, and that a
processing variant (different master pass, different pairing) can be
built without touching the standard products.

## The whole thing, one command

Deploy (stamps the code version into the tree - bare rsync does not):

    opr_vvel/server/deploy.sh

Then, on the server, detached so it survives a dropped connection:

    setsid nohup <code>/opr_vvel/server/reproduce_all.sh \
        < /dev/null > /dev/null 2>&1 &

It refuses to start unless the environment matches
`environment.pin` (OPR toolbox SHA including the headless docking patch,
MATLAB version), then rebuilds multipass for GL1-GL4, verifies the
rebuilds bit-identical to the archive, rebuilds every vvel variant the
standing analyses read (main / sequential / all-pairs / 2.5 km blocks /
master-override), and reruns the standing diagnostics and figures - the
exact list is the stage header of `reproduce_all.sh`; scripts outside it
(the elasticity chain, the maps, the movie) are run manually. Per-stage logs
land in `.../vvel/logs/repro_*.log`, the summary in `reproduce_all.log`,
and `reproduce_all_done` appears at the end. Idempotent: finished outputs
are skipped, so rerunning after an interruption resumes. Every vvel
product records `code_version`, `input_fn`/`input_bytes`/`input_mtime`
and `matlab_version`, so any file can answer "what built me" on its own.

EAGER_2022 is the standing exception: its archived combine_passes input
was replaced after the product was built (10 passes on disk against the
product's 13), so stage 1 cannot include it and `run_multipass_scratch.m`
refuses it with that explanation. Its DERIVED products still rebuild from
the archived multipass03 file.

`scripts/figures/tidal_evidence.m` is NOT part of the reproduction
pipeline. It is the historical retraction figure and requires the
preserved legacy `CSARP_vvel_v2` products (the retired scalar-coalignment
generation), which the current per-column code cannot and should not
rebuild. Run it manually when those products are present.

## The chain

```
combine_passes (archived)          .../CSARP_multipass/<name>.mat
  └─ multipass comp_mode 3         run_multipass_scratch.m  -> <name>_multipass03.mat
       └─ vvel (this repo)         run_vvel_scratch.m       -> CSARP_vvel*/<name>_vvel_i_j.mat
            └─ analysis            scripts/figures, scripts/diagnostics
```

## Reproducing a multipass product

One product per MATLAB session - `multipass.m` is a SCRIPT and leaks its
workspace, so a second product in the same session inherits the first
one's `coregistration_time_shift` (whose length may not even match):

    /opt/sw/matlab/2024b/bin/matlab -batch \
      "product='EAGER_2022_GL3'; run('<code>/opr_vvel/server/run_multipass_scratch.m')"

Then prove it:

    /opt/sw/matlab/2024b/bin/matlab -batch \
      "PRODUCT='EAGER_2022_GL3'; run('<code>/scripts/diagnostics/verify_multipass_rerun.m')"

`verify_multipass_rerun` requires the `data` cube to be bit-identical to
the archive, tolerates new `pass` fields the archive predates and
projection rounding below 1e-6 m, and errors on anything else. All four
GL products passed this on 2026-08-04 (the only toolbox commit touching
`+multipass` since the archived products, 6fcec75f, adds three `pass`
fields and nothing in the comp_mode 3 data path).

## Where the calibration numbers come from

The frozen `coregistration_time_shift` and `equalization` vectors in
`run_multipass_scratch.m` are copied verbatim from
`mem1:~/scripts/multipass_eastwind/run_multipass_EAGER.m`, which built
the archived products. They are the outputs of the full three-run
calibration workflow:

1. comp_mode 2 - coregistration estimate, prints per-pass time shifts
2. comp_mode 1 - equalization estimate, prints per-pass complex gains
3. comp_mode 3 - the differential product, consuming both vectors

The printed estimates are pasted between stages by hand. The vectors are
tied to the pass set AND to the master: any changed pass set must
re-derive both. `EAGER_2022` carries neither vector (wholly uncalibrated)
and `EAGER_2022_GL4` carries equalization but no time shifts (commented
out in the original); both facts are reproduced, not repaired.

## Processing variants

Master-pass override (standard build untouched; output suffixed `_mNN`):

    matlab -batch "product='EAGER_2022_GL3'; master_override=6; run('.../run_multipass_scratch.m')"

Keeping the frozen vectors under a master change is first-order correct:
the master enters `ref_z` as one constant across all passes, and a shift
common to every pass cancels in any interferometric pair. What a master
change DOES exercise - and what `scripts/diagnostics/master_sensitivity.m`
measures - is the resampling grid, the surface reference, and any
master-keyed residual.

vvel on a variant product:

    matlab -batch "product_tbl_override={'EAGER_2022_GL3_m06','all','vvel'}; \
      mp_dir_override='<scratch>/CSARP_multipass'; out_suffix='_netm06'; \
      run('.../run_vvel_scratch.m')"

Pairings: `'main'` (every pass against the product's master),
`'sequential'` (consecutive passes), `'all'` (every unordered pair - what
`vdef.invertNetwork` consumes).

## Long jobs and connectivity

`reproduce_all.sh` is already detached-safe when launched under
`setsid nohup ... < /dev/null`; poll for `reproduce_all_done` rather than
holding a connection open. The same pattern applies to any one-off job.
`pkill -f <pattern>` over ssh matches the remote shell issuing it - kill
by exact pid instead. (The ad-hoc run_chain*.sh scripts from Aug 2026 are
superseded by reproduce_all.sh and deleted.)

## Standing QC, in order of what they catch

| check | script | catches |
|---|---|---|
| rerun identity | `diagnostics/verify_multipass_rerun.m` | toolbox drift, wrong settings |
| loop closure | `diagnostics/closure.m` | inconsistent pairs (secular) |
| artefact tracking | `diagnostics/artefact_vs_signal.m` | tide-proportional processing residuals |
| between-build | `diagnostics/between_build.m` | reproducibility across products |
| master sensitivity | `diagnostics/master_sensitivity.m` | master-keyed systematics |
