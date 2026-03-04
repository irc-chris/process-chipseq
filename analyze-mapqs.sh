#!/bin/bash
#SBATCH --job-name=mapq_summary
#SBATCH --output=mapq_summary_%j.out
#SBATCH --error=mapq_summary_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=weka

export PATH=/mnt/altnas/work/ishawnia/empvsag/process-chipseq/samtools-1.19:$PATH

set -euo pipefail

bams=(
    # "/mnt/altnas/work/ragini/ENCODE-paper/brainstormanalysis/analysis-v2/phasechipseq/gm12878/CTCF_ENCSR000DRZ/CTCF_ENCSR000DRZ-untagged.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/DRZ/diploid-aligned/tmp/output.sorted.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/DRZ/haploid-aligned/untagged.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/SRR18530773/dip-aligned/aligned.sorted.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/SRR18530773/hap-aligned/untagged.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/DRZ/diploid-aligned/hap1.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/DRZ/diploid-aligned/hap2.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/DRZ/haploid-aligned/hap1.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/DRZ/haploid-aligned/hap2.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/SRR18530773/dip-aligned/hap1_renamed.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/SRR18530773/dip-aligned/hap2_renamed.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/SRR18530773/hap-aligned/hap1.bam"
    # "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/SRR18530773/hap-aligned/hap2.bam"
    "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/DRZ/haploid-aligned/tmp/output.sorted.bam"
    "/mnt/altnas/work/ishawnia/empvsag/process-chipseq/SRR18530773/tmp/output.sorted.bam"
)

printf "%-55s %12s %12s %12s %12s %12s %12s\n" \
    "BAM" "total" "m=0" "0<m<=10" "10<m<=30" "m>30" "m>10"
printf '%0.s-' {1..130}; echo

for bam in "${bams[@]}"; do
    name=$(basename "$bam" .bam)
    total=$(samtools view -c      "$bam")
    ge1=$(samtools view  -c -q 1  "$bam")
    ge11=$(samtools view -c -q 11 "$bam")
    ge31=$(samtools view -c -q 31 "$bam")

    m0=$(( total - ge1 ))
    m1_10=$(( ge1 - ge11 ))
    m11_30=$(( ge11 - ge31 ))
    mgt30=$ge31
    mgt10=$ge11

    printf "%-55s %12s %12s %12s %12s %12s %12s\n" \
        "$name" "$total" "$m0" "$m1_10" "$m11_30" "$mgt30" "$mgt10"
done
