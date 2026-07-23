#!/usr/bin/env bash
# Build the ddbj-gff wheel from the local sibling repo into this repo's root so
# `docker build -f Dockerfile.slim` can COPY it into the image (../gff_submission
# is outside the Docker build context, so it cannot be referenced directly).
#
# The wheel is NOT committed (see .gitignore). ddbj-gff is under active
# development and its version stays 0.1.0, so regenerate the wheel whenever
# ddbj-gff changes — this always bundles the current source.
#
# Usage:
#   scripts/build-ddbj-gff-wheel.sh          # build from ../gff_submission
#   DDBJ_GFF_DIR=/path/to/gff scripts/build-ddbj-gff-wheel.sh
#   PYTHON=python3.12 scripts/build-ddbj-gff-wheel.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GFF_DIR="${DDBJ_GFF_DIR:-$REPO_ROOT/../gff_submission}"
PYTHON="${PYTHON:-python3}"

if [ ! -f "$GFF_DIR/pyproject.toml" ]; then
    echo "error: ddbj-gff source not found at '$GFF_DIR'" >&2
    echo "       set DDBJ_GFF_DIR to your gff_submission checkout." >&2
    exit 1
fi

echo "[build-ddbj-gff-wheel] building wheel from: $GFF_DIR"
rm -f "$REPO_ROOT"/ddbj_gff-*.whl
"$PYTHON" -m pip wheel --no-deps -w "$REPO_ROOT" "$GFF_DIR"
echo "[build-ddbj-gff-wheel] wrote:"
ls -1 "$REPO_ROOT"/ddbj_gff-*.whl
