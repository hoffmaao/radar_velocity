#!/bin/bash
# Rerun multipass comp_mode 3 for the four EAGER GL products, one fresh
# MATLAB session per product - multipass.m is a script, not a function,
# and per-product state (e.g. coregistration_time_shift length) must not
# leak from one product into the next.
#
# Launch detached on mem1:
#   nohup bash /kucresis/scratch/hoffmana_sta/vvel/code/opr_vvel/server/run_multipass_all.sh \
#     > /kucresis/scratch/hoffmana_sta/vvel/multipass_rerun.log 2>&1 &
set -u
CODE=/kucresis/scratch/hoffmana_sta/vvel/code

for p in EAGER_2022_GL1 EAGER_2022_GL2 EAGER_2022_GL3 EAGER_2022_GL4; do
  echo "=== $p start $(date) ==="
  /opt/sw/matlab/2024b/bin/matlab -batch \
    "product='$p'; run('$CODE/opr_vvel/server/run_multipass_scratch.m')"
  echo "=== $p exit=$? $(date) ==="
done
echo "=== all done $(date) ==="
