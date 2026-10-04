#!/bin/bash
# Rebuild the four EAGER GL multipass products WITHOUT z-motion
# compensation (see run_multipass_scratch.m, zmotion_off). Per product the
# calibration is re-derived in three stages - shift estimate, equalization,
# product - each in a fresh MATLAB session, because multipass.m is a script
# and its workspace must not leak between stages or products. A failed
# stage stops that product; the others still run.
#
# Requires the multipass.m zmotion_comp_en switch in the server's OPR
# checkout. Outputs go to .../CSARP_multipass_nozc; the standard products
# are not touched. EAGER_2022 is not rebuildable (see the driver) and is a
# duplicate build of GL1 in any case.
#
# Launch detached on mem1:
#   nohup bash /kucresis/scratch/hoffmana_sta/vvel/code_ls_review/opr_vvel/server/run_multipass_nozc.sh \
#     > /kucresis/scratch/hoffmana_sta/vvel/multipass_nozc.log 2>&1 &
set -u
CODE=/kucresis/scratch/hoffmana_sta/vvel/code_ls_review

for p in EAGER_2022_GL1 EAGER_2022_GL2 EAGER_2022_GL3 EAGER_2022_GL4; do
  for s in coreg equalize product; do
    echo "=== $p $s start $(date) ==="
    /opt/sw/matlab/2024b/bin/matlab -batch \
      "product='$p'; zmotion_off=true; stage='$s'; run('$CODE/opr_vvel/server/run_multipass_scratch.m')"
    rc=$?
    echo "=== $p $s exit=$rc $(date) ==="
    if [ $rc -ne 0 ]; then echo "*** $p stopped at stage $s"; break; fi
  done
done
echo "=== all done $(date) ==="
