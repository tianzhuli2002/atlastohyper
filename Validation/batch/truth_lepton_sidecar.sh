#!/bin/bash
# Worker for one chunk (HTCondor or interactive): truth_lepton_sidecar.sh <chunk_index> <OUT> [REPO]
set -euo pipefail
INDEX=$1; OUT=$2; REPO=${3:-/home/leet/FF/atlastohyper}
LIST=$(printf "%s/lists/chunk_%04d.txt" "$OUT" "$INDEX")
PART=$(printf "%s/parts/part_%04d.npz" "$OUT" "$INDEX")
mkdir -p "$OUT/parts"
if [[ -f "$PART" && -f "${PART%.npz}.manifest.json" ]]; then echo "exists: $PART"; exit 0; fi
python3 "$REPO/Validation/build_truth_lepton_sidecar.py" build "$LIST" --output "$PART" --step-size "200 MB"
python3 "$REPO/Validation/build_truth_lepton_sidecar.py" validate "$PART" --output "${PART%.npz}.validation.json"
