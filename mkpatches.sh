#!/usr/bin/env bash

set -euo pipefail
# set -x

METIS_DIR="${1:-$HOME/projects/repos/METIS}"
GKLIB_DIR="${2:-$HOME/projects/repos/GKlib}"
DEST_DIR="${3:-"$(pwd)"}/normalization"

test -d "$METIS_DIR"
test -d "$GKLIB_DIR"
test -d "$DEST_DIR"

echo "metis dir: $METIS_DIR"
echo "gklib dir: $GKLIB_DIR"
echo "dest dir:  $DEST_DIR"


cd "$METIS_DIR"
git diff origin/master > "$DEST_DIR/metis-normalized.patch"

cd "$GKLIB_DIR"
git diff origin/master > "$DEST_DIR/gklib-normalized.patch"
