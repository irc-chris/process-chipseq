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
outputdir="$topDir/aligned"

HAP1_IN="$outputdir/hap1.bam"
HAP2_IN="$outputdir/hap2.bam"
HAP1_SORTED="$outputdir/hap1.sorted.bam"
HAP2_SORTED="$outputdir/hap2.sorted.bam"

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

echo "$(date) Sort and index complete."
