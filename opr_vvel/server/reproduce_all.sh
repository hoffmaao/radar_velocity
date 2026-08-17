#!/bin/bash
# REPRODUCE_ALL - the one entry point that rebuilds and verifies everything
# this project derives from the archived inputs.
#
#   ssh mem1 'setsid nohup <code>/opr_vvel/server/reproduce_all.sh \
#       < /dev/null > /dev/null 2>&1 &'
#
# Stages, each logged to .../vvel/logs/repro_<stage>.log:
#   0 env        refuse to run unless the pinned environment matches
#   1 multipass  rebuild GL1-GL4 from archived combine_passes + frozen
#                calibration (EAGER_2022 is NOT rebuildable - its archived
#                input was replaced after the build; its DERIVED stages
#                below still run from the archived multipass03 product)
#   2 verify     require the rebuilds bit-identical to the archive
#   3 vvel       every product variant the standing analyses read:
#                  main pairing -> CSARP_vvel_v3      (5 products)
#                  sequential   -> CSARP_vvel_seq_v3  (GL3, GL4)
#                  all pairs    -> CSARP_vvel_net     (5 products)
#                  all @ 2.5 km -> CSARP_vvel_net1k   (5 products)
#                  master-override all-pairs -> CSARP_vvel_netm (GL1-4)
#   4 diagnostics closure, artefact tracking, between-build, master
#                sensitivity, line means
#   5 figures    strain_rates (network), error_budget, tidal_evidence
#
# Idempotent: existing outputs are skipped (vvel via rerun_only, multipass
# via its own exists-check), so a rerun after an interruption resumes.
# FORCE_MULTIPASS=1 rebuilds stage 1 from scratch. The vvel stages have no
# force flag on purpose - delete the CSARP_vvel* directory you want rebuilt,
# which makes destruction explicit.
#
# Detached-safe: run under setsid; progress in the logs; touches
# .../vvel/logs/reproduce_all_done at the end with PASS/FAIL per stage.
set -uo pipefail

CODE=$(cd "$(dirname "$0")/../.." && pwd)
LOG=/kucresis/scratch/hoffmana_sta/vvel/logs
SCR=/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_multipass
PIN="$CODE/opr_vvel/server/environment.pin"
mkdir -p "$LOG"
rm -f "$LOG/reproduce_all_done"
SUMMARY="$LOG/reproduce_all.log"
: > "$SUMMARY"

note() { echo "$(date -u +%H:%M:%SZ) $*" | tee -a "$SUMMARY"; }
fail=0

# ---- stage 0: environment ------------------------------------------------
# shellcheck disable=SC1090
source "$PIN"
note "code version: $(cat "$CODE/code_version.txt" 2>/dev/null || echo 'UNSTAMPED - deploy with deploy.sh')"
sha=$(cd "$OPR_TOOLBOX_DIR" && git rev-parse HEAD)
dirty=$(cd "$OPR_TOOLBOX_DIR" && git status --porcelain | head -1)
mlver=$("$MATLAB_BIN" -batch "fprintf('%s',version)" 2>/dev/null | tail -1)
if [[ "$sha" != "$OPR_TOOLBOX_SHA" || -n "$dirty" ]]; then
  note "env FAIL: OPR toolbox is $sha (dirty='$dirty'), pinned $OPR_TOOLBOX_SHA"
  note "Update environment.pin deliberately if the toolbox moved on purpose."
  touch "$LOG/reproduce_all_done"; exit 1
fi
case "$mlver" in
  "$MATLAB_VERSION"*) note "env OK: toolbox $sha, MATLAB $mlver" ;;
  *) note "env FAIL: MATLAB $mlver, pinned $MATLAB_VERSION"; touch "$LOG/reproduce_all_done"; exit 1 ;;
esac

M="$MATLAB_BIN"
RUN_MP="$CODE/opr_vvel/server/run_multipass_scratch.m"
RUN_VV="$CODE/opr_vvel/server/run_vvel_scratch.m"

# ---- stage 1: multipass --------------------------------------------------
FORCE=${FORCE_MULTIPASS:-0}
for p in EAGER_2022_GL1 EAGER_2022_GL2 EAGER_2022_GL3 EAGER_2022_GL4; do
  $M -batch "product='$p'; force_rerun=$([ "$FORCE" = 1 ] && echo true || echo false); run('$RUN_MP')" \
    > "$LOG/repro_mp_$p.log" 2>&1 \
    && note "multipass $p OK" || { note "multipass $p FAIL"; fail=1; }
done
# master-override builds for the sensitivity diagnostic
for spec in "EAGER_2022_GL1:7" "EAGER_2022_GL2:8" "EAGER_2022_GL3:6" "EAGER_2022_GL4:7"; do
  p="${spec%%:*}"; m="${spec##*:}"
  $M -batch "product='$p'; master_override=$m; run('$RUN_MP')" \
    > "$LOG/repro_mp_${p}_m0${m}.log" 2>&1 \
    && note "multipass $p m0$m OK" || { note "multipass $p m0$m FAIL"; fail=1; }
done

# ---- stage 2: verify -----------------------------------------------------
for p in EAGER_2022_GL1 EAGER_2022_GL2 EAGER_2022_GL3 EAGER_2022_GL4; do
  $M -batch "PRODUCT='$p'; run('$CODE/scripts/diagnostics/verify_multipass_rerun.m')" \
    > "$LOG/repro_verify_$p.log" 2>&1 \
    && note "verify $p OK" || { note "verify $p FAIL"; fail=1; }
done

# ---- stage 3: vvel variants ---------------------------------------------
T5="{'EAGER_2022','PAIR','vvel';'EAGER_2022_GL1','PAIR','vvel';'EAGER_2022_GL2','PAIR','vvel';'EAGER_2022_GL3','PAIR','vvel';'EAGER_2022_GL4','PAIR','vvel'}"
vv() { # name, table, extra
  local name=$1 tbl=$2 extra=$3
  $M -batch "product_tbl_override=$tbl; $extra run('$RUN_VV')" \
    > "$LOG/repro_vvel_$name.log" 2>&1 \
    && note "vvel $name OK" || { note "vvel $name FAIL"; fail=1; }
}
vv v3    "${T5//PAIR/main}"       "out_suffix='_v3';"
vv seq   "{'EAGER_2022_GL3','sequential','vvel';'EAGER_2022_GL4','sequential','vvel'}" "out_suffix='_seq_v3';"
vv net   "${T5//PAIR/all}"        "out_suffix='_net';"
vv net1k "${T5//PAIR/all}"        "block_size_override=1000; out_suffix='_net1k';"
vv netm  "{'EAGER_2022_GL1_m07','all','vvel';'EAGER_2022_GL2_m08','all','vvel';'EAGER_2022_GL3_m06','all','vvel';'EAGER_2022_GL4_m07','all','vvel'}" \
         "mp_dir_override='$SCR'; out_suffix='_netm';"

# ---- stages 4+5: diagnostics and figures --------------------------------
for s in diagnostics/closure diagnostics/artefact_vs_signal diagnostics/between_build \
         diagnostics/master_sensitivity diagnostics/line_means \
         figures/strain_rates figures/error_budget figures/tidal_evidence; do
  n=$(basename "$s")
  $M -batch "run('$CODE/scripts/$s.m')" > "$LOG/repro_$n.log" 2>&1 \
    && note "$n OK" || { note "$n FAIL"; fail=1; }
done

note "=== reproduce_all finished, fail=$fail ==="
touch "$LOG/reproduce_all_done"
exit $fail
