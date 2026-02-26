#!/bin/bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --partition=weka2
#SBATCH --cpus-per-task=2

export PATH=/mnt/altnas/work/ishawnia/empvsag/process-chipseq/samtools-1.19:$PATH

set -euo pipefail

configFile=${1:?"Usage: $0 config.json"}
[ -f "$configFile" ] || { echo "Config not found: $configFile"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "bedtools not in PATH"; exit 1; }

topDir=$(jq -r '.topDir'    "$configFile")
bedfile1=$(jq -r '.bedfile1' "$configFile")
bedfile2=$(jq -r '.bedfile2' "$configFile")

[ -f "$bedfile1" ] || { echo "bedfile1 not found: $bedfile1"; exit 1; }
[ -f "$bedfile2" ] || { echo "bedfile2 not found: $bedfile2"; exit 1; }

outputdir="$topDir/aligned"
HAP1_BAM="$outputdir/hap1.sorted.bam"
HAP2_BAM="$outputdir/hap2.sorted.bam"

[ -f "$HAP1_BAM" ] || { echo "hap1 sorted BAM not found: $HAP1_BAM"; exit 1; }
[ -f "$HAP2_BAM" ] || { echo "hap2 sorted BAM not found: $HAP2_BAM"; exit 1; }

# In diploid mode hap1/hap2 BAMs come from two separate alignments.
# In haploid mode they come from splitting a single alignment via WhatsHap HP tags.
# Either way, hap1 is counted against bedfile1 and hap2 against bedfile2.

echo "$(date) Getting counts: hap1 vs bedfile1..."
bedtools multicov -bams "$HAP1_BAM" -bed "$bedfile1" > "$topDir/hap1_counts.tsv"
echo "Counts written to $topDir/hap1_counts.tsv"

echo "$(date) Getting counts: hap2 vs bedfile2..."
bedtools multicov -bams "$HAP2_BAM" -bed "$bedfile2" > "$topDir/hap2_counts.tsv"
echo "Counts written to $topDir/hap2_counts.tsv"

echo "$(date) Counts complete."
