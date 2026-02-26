#!/bin/bash
#SBATCH --job-name=build_diploid_index
#SBATCH --output=build_diploid_index_%j.out
#SBATCH --error=build_diploid_index_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=10G
#SBATCH --partition=weka

set -euo pipefail

export PATH=/mnt/altnas/work/ishawnia/empvsag/process-chipseq/samtools-1.19:$PATH

REFS="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs"
FA="$REFS/diploid.fa"

echo "$(date) - Starting BWA index..."
bwa index $FA
echo "$(date) - BWA index done."

echo "$(date) - Building samtools faidx..."
samtools faidx $FA
echo "$(date) - samtools faidx done."

echo "$(date) - Generating chromosome lengths file..."
cut -f1,2 ${FA}.fai > $REFS/diploid.len
echo "$(date) - All done."