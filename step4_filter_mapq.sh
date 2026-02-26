#!/bin/bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --partition=weka2
#SBATCH --cpus-per-task=4

export PATH=/mnt/altnas/work/ishawnia/empvsag/process-chipseq/samtools-1.19:$PATH

set -euo pipefail

configFile=${1:?"Usage: $0 config.json"}
[ -f "$configFile" ] || { echo "Config not found: $configFile"; exit 1; }

topDir=$(jq -r '.topDir' "$configFile")

TOTAL_BAM="$topDir/tmp/output.sorted.bam"
outputdir="$topDir/aligned"
base="$outputdir/output"

[ -f "$TOTAL_BAM" ] || { echo "Total BAM not found: $TOTAL_BAM"; exit 1; }

echo "$(date) Filtering total BAM by MAPQ..."

total=$(samtools view        -@ 4 -c       "$TOTAL_BAM")
mapq30_count=$(samtools view -@ 4 -c -q 30 "$TOTAL_BAM")
mapq10_count=$(samtools view -@ 4 -c -q 10 "$TOTAL_BAM")
mapq1_count=$(samtools view  -@ 4 -c -q 1  "$TOTAL_BAM")
mapq0_count=$((total - mapq1_count))

echo "  Total reads: $total"
echo "  MAPQ >= 30:  $mapq30_count"
echo "  MAPQ >= 10:  $mapq10_count"
echo "  MAPQ = 0:    $mapq0_count"

echo "$(date) Writing output_mapq30.bam..."
samtools view -@ 4 -b -q 30 "$TOTAL_BAM" > "${base}_mapq30.bam"

echo "$(date) Writing output_mapq10.bam..."
samtools view -@ 4 -b -q 10 "$TOTAL_BAM" > "${base}_mapq10.bam"

echo "$(date) Writing output_mapq0.bam (MAPQ=0 reads only)..."
samtools view -@ 4 -b -q 1 -U "${base}_mapq0.bam" "$TOTAL_BAM" > /dev/null

echo "$(date) Indexing filtered BAMs..."
samtools index "${base}_mapq30.bam"
samtools index "${base}_mapq10.bam"
samtools index "${base}_mapq0.bam"

echo "$(date) MAPQ filtering complete."
