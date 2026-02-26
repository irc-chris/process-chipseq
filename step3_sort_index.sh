#!/bin/bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --partition=weka2
#SBATCH --cpus-per-task=8

export PATH=/mnt/altnas/work/ishawnia/empvsag/process-chipseq/samtools-1.19:$PATH

set -euo pipefail

configFile=${1:?"Usage: $0 config.json"}
[ -f "$configFile" ] || { echo "Config not found: $configFile"; exit 1; }

topDir=$(jq -r '.topDir' "$configFile")
mode=$(jq -r '.mode // "haploid"' "$configFile")

outputdir="$topDir/aligned"

if [ "$mode" = "diploid" ]; then
    HAP1_IN="$outputdir/hap1_renamed.bam"
    HAP2_IN="$outputdir/hap2_renamed.bam"
    HAP1_SORTED="$outputdir/hap1_renamed.sorted.bam"
    HAP2_SORTED="$outputdir/hap2_renamed.sorted.bam"

    [ -f "$HAP1_IN" ] || { echo "hap1 BAM not found: $HAP1_IN"; exit 1; }
    [ -f "$HAP2_IN" ] || { echo "hap2 BAM not found: $HAP2_IN"; exit 1; }

    echo "$(date) Sorting hap1..."
    samtools sort -@ 8 -o "$HAP1_SORTED" "$HAP1_IN"
    echo "$(date) Sorting hap2..."
    samtools sort -@ 8 -o "$HAP2_SORTED" "$HAP2_IN"

    echo "$(date) Indexing hap1..."
    samtools index "$HAP1_SORTED"
    echo "$(date) Indexing hap2..."
    samtools index "$HAP2_SORTED"
else
    ALIGNED_IN="$topDir/tmp/output.sorted.bam"
    ALIGNED_SORTED="$outputdir/aligned.sorted.bam"

    [ -f "$ALIGNED_IN" ] || { echo "Aligned BAM not found: $ALIGNED_IN"; exit 1; }

    echo "$(date) Sorting aligned BAM..."
    samtools sort -@ 8 -o "$ALIGNED_SORTED" "$ALIGNED_IN"
    echo "$(date) Indexing..."
    samtools index "$ALIGNED_SORTED"
fi

echo "$(date) Sort and index complete."
