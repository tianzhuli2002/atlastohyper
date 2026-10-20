#!/bin/bash
# CPU/I-O only truth-lepton sidecar, one array task per chunk of input files.
# Usage:
#   split -l 50 -d -a 4 --additional-suffix=.txt $FILELIST $OUT/lists/chunk_   # 50 ROOT files per task
#   NCHUNK=$(ls $OUT/lists/chunk_*.txt | wc -l)
#   sbatch --array=0-$((NCHUNK-1))%40 --export=ALL,OUT=$OUT Validation/batch/truth_lepton_sidecar.sbatch
#SBATCH --job-name=truthlep
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.log
set -euo pipefail
: "${OUT:?set OUT to the sidecar work directory}"
REPO=${REPO:-/home/leet/FF/atlastohyper}
LIST=$(printf "%s/lists/chunk_%04d.txt" "$OUT" "$SLURM_ARRAY_TASK_ID")
PART=$(printf "%s/parts/part_%04d.npz" "$OUT" "$SLURM_ARRAY_TASK_ID")
mkdir -p "$OUT/parts"
if [[ -f "$PART" && -f "${PART%.npz}.manifest.json" ]]; then
  echo "exists: $PART"; exit 0
fi
python3 "$REPO/Validation/build_truth_lepton_sidecar.py" build "$LIST" --output "$PART" --step-size "200 MB"
python3 "$REPO/Validation/build_truth_lepton_sidecar.py" validate "$PART" --output "${PART%.npz}.validation.json" || {
  echo "validation failed for $PART"; exit 3; }
