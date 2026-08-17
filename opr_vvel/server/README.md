# Server workflow: reproducing every product from the archive

Everything below runs on mem1 (`ssh hoffmana_sta@mem1.cresis.ku.edu`,
VPN required). The goal of this file is that any product in the chain can
be regenerated from the archived inputs and PROVEN identical, and that a
processing variant (different master pass, different pairing) can be
built without touching the standard products.

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

Launch chains detached so they survive a dropped ssh session:

    setsid nohup <code>/run_chain.sh < /dev/null > /dev/null 2>&1 &

`run_chain.sh` writes per-step logs to `.../vvel/logs/` and touches
`chain_done` when finished; poll for that file rather than holding a
connection open. `pkill -f <pattern>` over ssh matches the remote shell
issuing it - kill by exact pid instead.

## Standing QC, in order of what they catch

| check | script | catches |
|---|---|---|
| rerun identity | `diagnostics/verify_multipass_rerun.m` | toolbox drift, wrong settings |
| loop closure | `diagnostics/closure.m` | inconsistent pairs (secular) |
| artefact tracking | `diagnostics/artefact_vs_signal.m` | tide-proportional processing residuals |
| between-build | `diagnostics/between_build.m` | reproducibility across products |
| master sensitivity | `diagnostics/master_sensitivity.m` | master-keyed systematics |
