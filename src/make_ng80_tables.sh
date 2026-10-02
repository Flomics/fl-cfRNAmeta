#!/usr/bin/env bash
#
# Generates the two NG80 tables read by boxplots_fig2.R and diversity_scatterplots.R:
#   tables/genes_contributing_to_percentage_reads.tsv                protein-coding genes only ("_pc")
#   tables/genes_contributing_to_percentage_reads_no_spike_ins.tsv   all genes except SIRV/ERCC
#
# Usage: bash src/make_ng80_tables.sh <count_matrix.tsv> <gene_info.tsv> <spikein_ids.tsv>
#   count_matrix.tsv  unfiltered RAW COUNTS matrix (gene_id, gene_name, samples...),
#                     e.g. results/<dataset>/expression_matrix/gene_raw_counts.tsv
#                     NG80 is defined on raw counts: TPM divides by gene length, which
#                     reshapes the per-sample distribution and changes the metric.
#                     Not the filterByExpr_* version either, which drops genes.
#   gene_info.tsv     GENCODE annotation with gene_id and gene_type columns.
#   spikein_ids.tsv   spike-in gene_ids, e.g. snakeDA's
#                     resources/gene_list/<gtf_id>/spikein_gene_ids.tsv (99 ERCC+SIRV).
#                     Do not match on /ERCC/ instead: that also catches the human
#                     ERCC1-8 DNA-repair genes.
#
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: bash src/make_ng80_tables.sh <count_matrix.tsv> <gene_info.tsv> <spikein_ids.tsv>" >&2
    exit 1
fi

MATRIX=$1
GENE_INFO=$2
SPIKEINS=$3
REPO=$(cd "$(dirname "$0")/.." && pwd)

for f in "$MATRIX" "$GENE_INFO" "$SPIKEINS"; do
    [ -r "$f" ] || { echo "Error: cannot read $f" >&2; exit 1; }
done

mkdir -p "$REPO/tables"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

column_by_name() {
    awk -F'\t' -v name="$2" '
        NR==1 { for (i=1; i<=NF; i++) if ($i == name) col=i
                if (!col) { print "Error: column " name " not found" > "/dev/stderr"; exit 1 }
                next }
        { print $col }' "$1"
}

# ng80.R always writes genes_contributing_to_percentage_reads.tsv into the current
# directory, so each run happens in WORK and the result is moved into tables/.
run_ng80() {
    ( cd "$WORK" && Rscript "$REPO/src/ng80.R" "$1" )
    mv "$WORK/genes_contributing_to_percentage_reads.tsv" "$REPO/tables/$2"
    echo "  $2: $(( $(wc -l < "$REPO/tables/$2") - 1 )) samples"
}

# protein-coding only: filter_gene_ids.py drops the IDs it is given
python3 "$REPO/src/filter_rows.py" "$GENE_INFO" gene_type '!=' protein_coding \
    | column_by_name /dev/stdin gene_id > "$WORK/drop_not_protein_coding.txt"
python3 "$REPO/src/filter_gene_ids.py" "$MATRIX" "$WORK/drop_not_protein_coding.txt" \
    > "$WORK/matrix_protein_coding.tsv"
run_ng80 "$WORK/matrix_protein_coding.tsv" genes_contributing_to_percentage_reads.tsv

# all genes except spike-ins
column_by_name "$SPIKEINS" gene_id > "$WORK/drop_spike_ins.txt"
python3 "$REPO/src/filter_gene_ids.py" "$MATRIX" "$WORK/drop_spike_ins.txt" \
    > "$WORK/matrix_no_spike_ins.tsv"
run_ng80 "$WORK/matrix_no_spike_ins.tsv" genes_contributing_to_percentage_reads_no_spike_ins.tsv
