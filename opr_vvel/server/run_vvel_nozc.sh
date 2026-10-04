#!/bin/bash
# All-pairs vvel network over the multipass products built WITHOUT z-motion
# compensation (run_multipass_nozc.sh). Writes the two network builds the
# standing analyses read, under _nozc names so the compensated builds stay
# for comparison:
#   CSARP_vvel_net_nozc       500 m blocks (200 columns)
#   CSARP_vvel_net_fine_nozc  125 m blocks (50 columns)
# EAGER_2022 is absent: it duplicates GL1 and cannot be rebuilt.
#
# Launch detached on mem1 once run_multipass_nozc.sh has finished:
#   nohup bash /kucresis/scratch/hoffmana_sta/vvel/code_ls_review/opr_vvel/server/run_vvel_nozc.sh \
#     > /kucresis/scratch/hoffmana_sta/vvel/vvel_nozc.log 2>&1 &
set -u
CODE=/kucresis/scratch/hoffmana_sta/vvel/code_ls_review
MP=/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_multipass_nozc
TBL="{'EAGER_2022_GL1','all','vvel';'EAGER_2022_GL2','all','vvel';'EAGER_2022_GL3','all','vvel';'EAGER_2022_GL4','all','vvel'}"

for p in EAGER_2022_GL1 EAGER_2022_GL2 EAGER_2022_GL3 EAGER_2022_GL4; do
  if [ ! -f "$MP/${p}_multipass03.mat" ]; then echo "*** missing $MP/${p}_multipass03.mat"; exit 1; fi
done
for spec in "_net_nozc:" "_net_fine_nozc:block_size_override=50;"; do
  sfx="${spec%%:*}"; extra="${spec#*:}"
  echo "=== vvel$sfx start $(date) ==="
  /opt/sw/matlab/2024b/bin/matlab -batch \
    "product_tbl_override=$TBL; mp_dir_override='$MP'; $extra out_suffix='$sfx'; run('$CODE/opr_vvel/server/run_vvel_scratch.m')"
  echo "=== vvel$sfx exit=$? $(date) ==="
done
echo "=== all done $(date) ==="
