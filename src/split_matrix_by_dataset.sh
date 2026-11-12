#!/usr/bin/env bash
#
# Splits an all-batches expression matrix into one matrix + sampleinfo per dataset, so
# src/dimred_perdataset.py can run per dataset without one quantification run per dataset.
#
# Valid because gene_tpm_norm is normalised per sample: every column sums to 1e6, so a
# subset of columns is exactly the TPM that dataset would give on its own. This would NOT
# hold for TMM, whose factors are fitted across the whole sample set.
#
# Usage: bash src/split_matrix_by_dataset.sh <matrix.tsv> <sampleinfo.tsv> <out_dir> [split_col]
#   matrix.tsv      gene x sample matrix (gene_id, gene_name, samples...), e.g.
#                   results/all-batches_meta_filtered/expression_matrix/gene_tpm_norm.tsv
#   sampleinfo.tsv  sample metadata; needs a sample_name column plus the split column
#   out_dir         written as <out_dir>/<group>/{matrix.tsv,sampleinfo.tsv}
#   split_col       granularity, default dataset_short_name (18 groups, one per dataset).
#                   dataset_batch gives 26 groups instead, splitting the datasets that
#                   were sequenced in more than one batch.
#
# Writes about as much as the source matrix, and can be deleted once the plots exist.
#
set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: bash src/split_matrix_by_dataset.sh <matrix.tsv> <sampleinfo.tsv> <out_dir> [split_col]" >&2
    exit 1
fi

MATRIX=$1
SAMPLEINFO=$2
OUT=$3
SPLIT_COL=${4:-dataset_short_name}

for f in "$MATRIX" "$SAMPLEINFO"; do
    [ -r "$f" ] || { echo "Error: cannot read $f" >&2; exit 1; }
done

# gawk and mawk both stream this comfortably; BSD awk hits its open-file limit on 26 groups
AWK=$(command -v gawk || command -v mawk || command -v awk)

# checked here rather than inside the awks below, whose exit still runs their END block
for col in "$SPLIT_COL" sample_name; do
    "$AWK" -F'\t' -v name="$col" \
        'NR==1 { for (i=1; i<=NF; i++) if ($i == name) found=1; exit !found }' "$SAMPLEINFO" \
        || { echo "Error: column $col not in $SAMPLEINFO" >&2; exit 1; }
done

mkdir -p "$OUT"
echo "matrix : $MATRIX"
echo "out    : $OUT"
echo "split  : $SPLIT_COL"

# one directory per group, and its slice of the sampleinfo
"$AWK" -v OUT="$OUT" -v SC="$SPLIT_COL" 'BEGIN{FS=OFS="\t"}
NR==1 { for (i=1; i<=NF; i++) if ($i == SC) c=i
        hdr=$0; next }
{ g=$c; if (g=="" || g=="NA") next
  if (!(g in seen)) { seen[g]=1
                      system("mkdir -p \"" OUT "/" g "\"")
                      print hdr > (OUT "/" g "/sampleinfo.tsv") }
  print > (OUT "/" g "/sampleinfo.tsv") }
END { print "groups : " length(seen) }' "$SAMPLEINFO"

# per-group matrix, one streaming pass over the big file
"$AWK" -v OUT="$OUT" -v SC="$SPLIT_COL" 'BEGIN{FS=OFS="\t"}
# first file: sample_name -> group
NR==FNR {
  if (FNR==1) { for (i=1; i<=NF; i++) { if ($i==SC) c=i; if ($i=="sample_name") sn=i }
                next }
  grp[$sn]=$c; next
}
# matrix header: work out which columns belong to which group
FNR==1 {
  for (i=3; i<=NF; i++) {
    g=grp[$i]
    if (g=="" || g=="NA") { unmapped++; continue }
    if (!(g in seen)) { seen[g]=1; ord[++ng]=g; f[g]=OUT "/" g "/matrix.tsv" }
    col[g,++n[g]]=i
  }
  if (unmapped) print "warning: " unmapped " matrix columns had no group, skipped" > "/dev/stderr"
}
{
  for (k=1; k<=ng; k++) { g=ord[k]; line=$1 OFS $2
    for (j=1; j<=n[g]; j++) line=line OFS $(col[g,j])
    print line > f[g] }
}
END { for (k=1; k<=ng; k++) { g=ord[k]; close(f[g]); print "  " g "\t" n[g] " samples" } }
' "$SAMPLEINFO" "$MATRIX"
