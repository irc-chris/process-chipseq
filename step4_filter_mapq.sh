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
mode=$(jq -r '.mode // "haploid"' "$configFile")

outputdir="$topDir/aligned"

if [ "$mode" = "diploid" ]; then
    BAMS=("$outputdir/hap1_renamed.sorted.bam" "$outputdir/hap2_renamed.sorted.bam")
else
    BAMS=("$outputdir/aligned.sorted.bam")
fi

for HAP_BAM in "${BAMS[@]}"; do
    [ -f "$HAP_BAM" ] || { echo "BAM not found: $HAP_BAM"; exit 1; }

    base="${HAP_BAM%.bam}"
    hap_name=$(basename "$base")

    echo "$(date) --- Filtering $hap_name ---"
    total=$(samtools view       -@ 4 -c       "$HAP_BAM")
    mapq30_count=$(samtools view -@ 4 -c -q 30 "$HAP_BAM")
    mapq10_count=$(samtools view -@ 4 -c -q 10 "$HAP_BAM")
    mapq1_count=$(samtools view  -@ 4 -c -q 1  "$HAP_BAM")
    mapq0_count=$((total - mapq1_count))

    echo "  Total reads: $total"
    echo "  MAPQ >= 30:  $mapq30_count"
    echo "  MAPQ >= 10:  $mapq10_count"
    echo "  MAPQ = 0:    $mapq0_count"

    echo "  Writing ${hap_name}_mapq30.bam..."
    samtools view -@ 4 -b -q 30 "$HAP_BAM" > "${base}_mapq30.bam"

    echo "  Writing ${hap_name}_mapq10.bam..."
    samtools view -@ 4 -b -q 10 "$HAP_BAM" > "${base}_mapq10.bam"

    echo "  Writing ${hap_name}_mapq0.bam (MAPQ=0 reads only)..."
    samtools view -@ 4 -b -q 1 -U "${base}_mapq0.bam" "$HAP_BAM" > /dev/null

    echo "  Indexing..."
    samtools index "${base}_mapq30.bam"
    samtools index "${base}_mapq10.bam"
    samtools index "${base}_mapq0.bam"
done

echo "$(date) MAPQ filtering complete."
