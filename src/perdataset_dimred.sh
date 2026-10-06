#!/usr/bin/env bash
#
# PCA + tSNE for every dataset separately, from the all-batches expression matrix.
#
# Two steps:
#   1. split the matrix into one matrix + sampleinfo per dataset (split_matrix_by_dataset.sh)
#   2. run PCA + tSNE on each of them (dimred_perdataset.py)
#
# Usage: bash src/perdataset_dimred.sh <matrix.tsv> <sampleinfo.tsv> [out_dir]
#   matrix.tsv      all-batches gene x sample matrix, e.g.
#                   results/all-batches_meta_filtered/expression_matrix/gene_tpm_norm.tsv
#   sampleinfo.tsv  the matching sample metadata
#   out_dir         default ./perdataset_dimred
#
# Plots land in <out_dir>/plots/<dataset>/{pca,tsne}/2d-{pca,tsne}_<variable>.{png,pdf}.
# The split under <out_dir>/split is reused if it is already there, and can be deleted
# once the plots exist.
#
# Overridable: SPLIT_COL (dataset_batch, per batch | dataset_short_name, per lab),
# PYTHON and MARKER_SIZE.
#
# Needs pandas, numpy, matplotlib, seaborn, scikit-learn.
#
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: bash src/perdataset_dimred.sh <matrix.tsv> <sampleinfo.tsv> [out_dir]" >&2
    exit 1
fi

MATRIX=$1
SAMPLEINFO=$2
OUT=${3:-$PWD/perdataset_dimred}
SRC=$(cd "$(dirname "$0")" && pwd)
PYTHON=${PYTHON:-python}

for f in "$MATRIX" "$SAMPLEINFO"; do
    [ -r "$f" ] || { echo "Error: cannot read $f" >&2; exit 1; }
done

"$PYTHON" -c "import pandas, numpy, matplotlib, seaborn, sklearn" 2>/dev/null || {
    echo "Error: $PYTHON is missing pandas, numpy, matplotlib, seaborn or scikit-learn." >&2
    exit 1; }

SPLIT="$OUT/split"
PLOTS="$OUT/plots"
mkdir -p "$OUT"

echo "split : $SPLIT"
echo "plots : $PLOTS"
echo

if [ -d "$SPLIT" ] && [ -n "$(ls -A "$SPLIT" 2>/dev/null)" ]; then
    echo "[1/2] split already present, reusing it"
else
    echo "[1/2] splitting by ${SPLIT_COL:-dataset_batch}"
    bash "$SRC/split_matrix_by_dataset.sh" "$MATRIX" "$SAMPLEINFO" "$SPLIT" \
        "${SPLIT_COL:-dataset_batch}"
fi

echo
echo "[2/2] PCA + tSNE per dataset"
ok=0; failed=()
for d in "$SPLIT"/*/; do
    ds=$(basename "$d")
    [ -r "$d/matrix.tsv" ] || { echo "skipping $ds: no matrix.tsv"; continue; }
    echo "==================== $ds ===================="
    # one dataset failing (too few samples, a missing column) must not stop the rest
    if "$PYTHON" "$SRC/dimred_perdataset.py" \
            --matrix "$d/matrix.tsv" --sampleinfo "$d/sampleinfo.tsv" \
            --label "$ds" --out-dir "$PLOTS" \
            --marker-size "${MARKER_SIZE:-12}"; then
        ok=$((ok + 1))
    else
        failed+=("$ds")
        echo "!!! FAILED: $ds (continuing)"
    fi
done

echo
echo "done. ok=$ok failed=${#failed[@]}"
[ "${#failed[@]}" -gt 0 ] && printf 'failed: %s\n' "${failed[*]}"
echo "plots under $PLOTS/<dataset>/{pca,tsne}/"
