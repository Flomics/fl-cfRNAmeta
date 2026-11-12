#!/bin/bash

BAM="$1"
BASENAME=$(basename "$BAM")
SAMPLE_ID="${BASENAME%%.*}"

# Count total uniquely mapped reads (-F 4 = exclude unmapped, -q 255 = MAPQ 255 only)
TOTAL=$(samtools view -F 4 -q 255 "$BAM" | wc -l)

# Count uniquely mapped spliced reads (CIGAR contains 'N')
SPLICED=$(samtools view -F 4 -q 255 "$BAM" | awk '$6 ~ /N/' | wc -l)

if [ "$TOTAL" -eq 0 ]; then
  PERCENT="0.00"
else
  PERCENT=$(awk -v s="$SPLICED" -v t="$TOTAL" 'BEGIN { printf "%.2f", (s / t) * 100 }')
fi

OUTFILE="${SAMPLE_ID}_spliced_reads_summary.tsv"
echo -e "sample_id\tspliced_reads\tpercent_spliced_reads" > "$OUTFILE"
echo -e "${SAMPLE_ID}\t${SPLICED}\t${PERCENT}" >> "$OUTFILE"
