#!/usr/bin/env bash
#
# make_ng80_tables.sh
#
# Generates the two NG80 tables that boxplots_fig2.R and diversity_scatterplots.R read:
#
#   tables/genes_contributing_to_percentage_reads.tsv
#       NG80 over PROTEIN-CODING genes only. Read as `ng_only_mrna` and used for the
#       "_pc" columns (boxplots_fig2.R, section "NG80 protein coding").
#
#   tables/genes_contributing_to_percentage_reads_no_spike_ins.tsv
#       NG80 over all genes EXCEPT the spike-ins (SIRV / ERCC).
#
# Both are produced by src/ng80.R, which always writes
# 'genes_contributing_to_percentage_reads.tsv' into the current directory, so each run
# happens in its own temporary directory and the result is moved into tables/.
#
# Usage:
#   bash src/make_ng80_tables.sh <count_matrix.tsv> <gene_info.tsv>
#
#   count_matrix.tsv  Gene-level matrix, UNFILTERED (all ~61.6k genes), columns:
#                     gene_id, gene_name, <one column per sample>.
#                     snakeDA writes it to
#                     results/<dataset>/expression_matrix/gene_tpm_norm.tsv
#                     Note: NOT the filterByExpr_* version, which drops genes and would
#                     make the NG80 counts meaningless.
#
#   gene_info.tsv     GENCODE annotation with 'gene_id' and 'gene_type' columns
#                     (e.g. GENCODEv39_biotype_gene_info.tsv).
#
set -euo pipefail

if [ "$#" -ne 2 ]; then
    sed -n '2,33p' "$0" >&2
    exit 1
fi

MATRIX=$1
GENE_INFO=$2
REPO=$(cd "$(dirname "$0")/.." && pwd)

for f in "$MATRIX" "$GENE_INFO"; do
    [ -r "$f" ] || { echo "Error: cannot read $f" >&2; exit 1; }
done

mkdir -p "$REPO/tables"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# Pull a named column out of a TSV.
column_by_name() {
    awk -F'\t' -v name="$2" '
        NR==1 { for (i=1; i<=NF; i++) if ($i == name) col=i
                if (!col) { print "Error: column " name " not found" > "/dev/stderr"; exit 1 }
                next }
        { print $col }' "$1"
}

n_genes=$(( $(wc -l < "$MATRIX") - 1 ))
echo "Count matrix: $MATRIX ($n_genes genes)"

# ---------------------------------------------------------------------------
# 1. NG80 over protein-coding genes only
#    filter_gene_ids.py removes the IDs it is given, so the list is every gene
#    whose gene_type is NOT protein_coding.
# ---------------------------------------------------------------------------
python3 "$REPO/src/filter_rows.py" "$GENE_INFO" gene_type '!=' protein_coding \
    | column_by_name /dev/stdin gene_id > "$WORK/ids_not_protein_coding.txt"
echo "  non-protein-coding genes to drop: $(wc -l < "$WORK/ids_not_protein_coding.txt")"

python3 "$REPO/src/filter_gene_ids.py" "$MATRIX" "$WORK/ids_not_protein_coding.txt" \
    > "$WORK/matrix_protein_coding.tsv"
echo "  protein-coding matrix: $(( $(wc -l < "$WORK/matrix_protein_coding.tsv") - 1 )) genes"

( cd "$WORK" && Rscript "$REPO/src/ng80.R" "$WORK/matrix_protein_coding.tsv" )
mv "$WORK/genes_contributing_to_percentage_reads.tsv" \
   "$REPO/tables/genes_contributing_to_percentage_reads.tsv"

# ---------------------------------------------------------------------------
# 2. NG80 over all genes except the spike-ins
# ---------------------------------------------------------------------------
awk -F'\t' 'NR>1 && ($1 ~ /SIRV|ERCC/ || $2 ~ /SIRV|ERCC/) { print $1 }' "$MATRIX" \
    > "$WORK/ids_spike_ins.txt"
echo "  spike-in genes to drop: $(wc -l < "$WORK/ids_spike_ins.txt")"

python3 "$REPO/src/filter_gene_ids.py" "$MATRIX" "$WORK/ids_spike_ins.txt" \
    > "$WORK/matrix_no_spike_ins.tsv"
echo "  no-spike-in matrix: $(( $(wc -l < "$WORK/matrix_no_spike_ins.tsv") - 1 )) genes"

( cd "$WORK" && Rscript "$REPO/src/ng80.R" "$WORK/matrix_no_spike_ins.tsv" )
mv "$WORK/genes_contributing_to_percentage_reads.tsv" \
   "$REPO/tables/genes_contributing_to_percentage_reads_no_spike_ins.tsv"

# ---------------------------------------------------------------------------
echo
echo "Done. Samples per table:"
for t in genes_contributing_to_percentage_reads genes_contributing_to_percentage_reads_no_spike_ins; do
    printf '  %-52s %s\n' "tables/$t.tsv" "$(( $(wc -l < "$REPO/tables/$t.tsv") - 1 ))"
done
