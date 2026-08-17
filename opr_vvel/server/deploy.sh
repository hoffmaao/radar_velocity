#!/bin/bash
# Deploy the repo to mem1, stamping the code version as it goes.
#
# This is THE deployment path. The deployed tree is not a git checkout, so
# products built there cannot learn their own version from git; this script
# writes code_version.txt at the deployed root and vvel_code_version.m
# reads it into every product. Deploying by bare rsync skips the stamp and
# products then record 'unknown' - visible, but useless for provenance.
#
# Refuses a dirty working tree unless FORCE_DIRTY=1, because a product
# stamped with a SHA that does not match the code that built it is worse
# than one stamped unknown.
set -euo pipefail
cd "$(dirname "$0")/../.."

HOST=hoffmana_sta@mem1.cresis.ku.edu
DEST=/kucresis/scratch/hoffmana_sta/vvel/code

DESC=$(git describe --always --dirty --tags)
if [[ "$DESC" == *-dirty && "${FORCE_DIRTY:-0}" != 1 ]]; then
  echo "Working tree is dirty ($DESC). Commit first, or FORCE_DIRTY=1 to override." >&2
  exit 1
fi
STAMP="git-$DESC (branch $(git rev-parse --abbrev-ref HEAD), deployed $(date -u +%Y-%m-%dT%H:%M:%SZ))"
echo "$STAMP" > /tmp/code_version.txt
echo "deploying: $STAMP"

rsync -a --delete +vdef opr_vvel "$HOST:$DEST/"
rsync -a scripts "$HOST:$DEST/"
scp -q /tmp/code_version.txt "$HOST:$DEST/code_version.txt"
ssh "$HOST" "chmod +x $DEST/opr_vvel/server/reproduce_all.sh 2>/dev/null || true"
echo "deployed."
