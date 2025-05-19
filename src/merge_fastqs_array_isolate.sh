#!/bin/bash
#SBATCH --job-name=merge_fastqs
#SBATCH --output=logs/merge_fastqs_%A_%a.out
#SBATCH --error=logs/merge_fastqs_%A_%a.err
#SBATCH --array=1-250
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --no-requeue
#SBATCH -p genoa64
#SBATCH --qos='long'

set -euo pipefail

_term() {
         echo "Caught SIGTERM signal!"
         kill -s SIGTERM $pid
         wait $pid
 }


mkdir -p logs merged_fastqs

ISOLATE_LIST="isolate_list.txt"
MAPPING_FILE="isolate_to_runs.txt"
OUTPUT_DIR="merged_fastqs"

ISOLATE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ISOLATE_LIST")

echo "Merging for isolate group: $ISOLATE"

RUNS=$(awk -v isolate="$ISOLATE" '$1 == isolate {print $2}' "$MAPPING_FILE")

FILES_1=""
FILES_2=""

for RUN in $RUNS; do
    FILE_1="${RUN}_1.fastq.gz"
    FILE_2="${RUN}_2.fastq.gz"

    if [[ -e "$FILE_1" && -e "$FILE_2" ]]; then
        FILES_1="$FILES_1 $FILE_1"
        FILES_2="$FILES_2 $FILE_2"
    else
        echo "Warning: Missing files for $RUN, skipping..."
    fi
done

if [[ -n "$FILES_1" && -n "$FILES_2" ]]; then
    cat $FILES_1 > "${OUTPUT_DIR}/${ISOLATE}_1.fastq.gz"
    cat $FILES_2 > "${OUTPUT_DIR}/${ISOLATE}_2.fastq.gz"
    echo "Finished merging $ISOLATE."
else
    echo "No files to merge for $ISOLATE."
fi
