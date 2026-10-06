#!/usr/bin/env bash
# Score the 10 example Alus as deletions, over all six sequence augmentations
# (shifts of 0, +1 and -1 bp, each forward and reverse complemented).
set -euo pipefail

# Paths - edit for your environment
SUPREMO_DIR="${SUPREMO_DIR:-/pollard/home/szhang20/akita_variant_scoring}"
FASTA="${SUPREMO_DIR}/data/hg38.fa"

EXAMPLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT="${EXAMPLE_DIR}/sample_input_10alus.txt"
OUT_DIR="${EXAMPLE_DIR}/output"
NAME="sample_10alus"

mkdir -p "${OUT_DIR}"

# SuPreMo.py is invoked by a relative path from its own checkout
cd "${SUPREMO_DIR}"

python scripts/SuPreMo.py "${INPUT}" \
    --file "${NAME}" \
    --dir "${OUT_DIR}" \
    --get_Akita_scores \
    --get_maps \
    --genome hg38 \
    --fa "${FASTA}" \
    --shift_by 0 1 -1 \
    --revcomp add_revcomp

echo "Wrote ${OUT_DIR}/${NAME}_scores and ${OUT_DIR}/${NAME}_maps.npy"
