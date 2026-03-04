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

topDir=$(jq -r '.topDir'     "$configFile")
mode=$(jq -r '.mode // "haploid"' "$configFile")
bedfile1=$(jq -r '.bedfile1' "$configFile")
bedfile2=$(jq -r '.bedfile2' "$configFile")

[ -f "$bedfile1" ] || { echo "bedfile1 not found: $bedfile1"; exit 1; }
[ -f "$bedfile2" ] || { echo "bedfile2 not found: $bedfile2"; exit 1; }

outputdir="$topDir/${mode}-aligned"
HAP1_BAM="$outputdir/hap1.sorted.bam"
HAP2_BAM="$outputdir/hap2.sorted.bam"

[ -f "$HAP1_BAM" ] || { echo "hap1 sorted BAM not found: $HAP1_BAM"; exit 1; }
[ -f "$HAP2_BAM" ] || { echo "hap2 sorted BAM not found: $HAP2_BAM"; exit 1; }

# MAPQ threshold: >=30 for haploid, >=10 for diploid
if [ "$mode" = "haploid" ]; then
    mapq=30
else
    mapq=10
fi

echo -e "CHR\tPOS1\tPOS2\thaplotype1" > "$outputdir/hap1_counts.tsv"
echo -e "CHR\tPOS1\tPOS2\thaplotype2" > "$outputdir/hap2_counts.tsv"

HAP1_FILTERED="$outputdir/hap1_mapq${mapq}.bam"
HAP2_FILTERED="$outputdir/hap2_mapq${mapq}.bam"

# --- Unfiltered counts ---
echo "$(date) Getting counts: hap1 vs bedfile1..."
bedtools multicov -bams "$HAP1_BAM" -bed "$bedfile1" >> "$outputdir/hap1_counts.tsv"
echo "Counts written to $outputdir/hap1_counts.tsv"

echo "$(date) Getting counts: hap2 vs bedfile2..."
bedtools multicov -bams "$HAP2_BAM" -bed "$bedfile2" > "$outputdir/hap2_counts.tsv"
echo "Counts written to $outputdir/hap2_counts.tsv"

# --- Filter hap BAMs by MAPQ (>= threshold) ---
echo "$(date) Filtering hap BAMs at MAPQ >= ${mapq}..."
samtools view -@ 2 -b -q "$mapq" "$HAP1_BAM" -o "$HAP1_FILTERED"
samtools index "$HAP1_FILTERED"
echo "Filtered BAM written to $HAP1_FILTERED"

samtools view -@ 2 -b -q "$mapq" "$HAP2_BAM" -o "$HAP2_FILTERED"
samtools index "$HAP2_FILTERED"
echo "Filtered BAM written to $HAP2_FILTERED"

# --- Filtered counts ---
echo "$(date) Getting filtered counts: hap1_mapq${mapq} vs bedfile1..."
bedtools multicov -bams "$HAP1_FILTERED" -bed "$bedfile1" > "$outputdir/hap1_mapq${mapq}_counts.tsv"
echo "Counts written to $outputdir/hap1_mapq${mapq}_counts.tsv"

echo "$(date) Getting filtered counts: hap2_mapq${mapq} vs bedfile2..."
bedtools multicov -bams "$HAP2_FILTERED" -bed "$bedfile2" > "$outputdir/hap2_mapq${mapq}_counts.tsv"
echo "Counts written to $outputdir/hap2_mapq${mapq}_counts.tsv"

echo "$(date) Counts complete."
