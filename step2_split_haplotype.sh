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

ALIGNED_BAM="$topDir/tmp/output.sorted.bam"
outputdir="$topDir/aligned"

[ -f "$ALIGNED_BAM" ] || { echo "Aligned BAM not found: $ALIGNED_BAM"; exit 1; }

HAP1_RAW="$outputdir/hap1.bam"
HAP2_RAW="$outputdir/hap2.bam"
HAP1_RENAMED="$outputdir/hap1_renamed.bam"
HAP2_RENAMED="$outputdir/hap2_renamed.bam"

echo "$(date) Extracting hap1 reads (chrom names matching _h1)..."
{ samtools view -@ 2 -H "$ALIGNED_BAM"; samtools view -@ 2 "$ALIGNED_BAM" | awk '$3 ~ /_h1/'; } \
    | samtools view -bS -o "$HAP1_RAW"

echo "$(date) Extracting hap2 reads (chrom names matching _h2)..."
{ samtools view -@ 2 -H "$ALIGNED_BAM"; samtools view -@ 2 "$ALIGNED_BAM" | awk '$3 ~ /_h2/'; } \
    | samtools view -bS -o "$HAP2_RAW"

echo "$(date) Stripping _h1 tags from chromosome names..."
samtools view -@ 2 -h "$HAP1_RAW" | sed 's/_h1//g' | samtools view -bS - > "$HAP1_RENAMED"

echo "$(date) Stripping _h2 tags from chromosome names..."
samtools view -@ 2 -h "$HAP2_RAW" | sed 's/_h2//g' | samtools view -bS - > "$HAP2_RENAMED"

echo "$(date) Haplotype splitting complete."
