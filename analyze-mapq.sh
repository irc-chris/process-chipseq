#!/bin/bash
#SBATCH --job-name=analyze_mapq
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=108:00:00
#SBATCH --mem=20G
#SBATCH --partition=weka2

OUTDIR=/mnt/altnas/work/ishawnia/empvsag/process-chipseq/filtering

# for bam in /mnt/altnas/work/ishawnia/empvsag/process-chipseq/DRZ/tmp/DRZ-output.bam \
#            /mnt/altnas/work/ishawnia/empvsag/process-chipseq/DRZ/tmp/DRZ-hap1.bam \
#            /mnt/altnas/work/ishawnia/empvsag/process-chipseq/DRZ/tmp/DRZ-hap2.bam \
for bam in /mnt/altnas/work/ragini/ENCODE-paper/brainstormanalysis/analysis-v2/phasechipseq/gm12878/CTCF_ENCSR000DRZ/CTCF_ENCSR000DRZ-haplotype1.bam \
           /mnt/altnas/work/ragini/ENCODE-paper/brainstormanalysis/analysis-v2/phasechipseq/gm12878/CTCF_ENCSR000DRZ/CTCF_ENCSR000DRZ-haplotype2.bam \
           /mnt/altnas/work/ishawnia/empvsag/process-chipseq/SRR18530773/aligned/hap1.bam \
           /mnt/altnas/work/ishawnia/empvsag/process-chipseq/SRR18530773/aligned/hap2.bam \
           /mnt/altnas/work/ishawnia/empvsag/process-chipseq/SRR18530773/aligned/aligned.sorted.bam; do

    echo "=== $bam ==="
    # Use just the basename for output, always written to OUTDIR
    base=$OUTDIR/$(basename ${bam%.bam})
    total=$(samtools view -c "$bam")
    mapq30=$(samtools view -c -q 30 "$bam")
    mapq10=$(samtools view -c -q 10 "$bam")
    mapq1=$(samtools view -c -q 1 "$bam")
    mapq0=$((total - mapq1))
    echo "Total reads: $total"
    echo "mapq >= 30:  $mapq30"
    echo "mapq >= 10:  $mapq10"
    echo "mapq = 0:    $mapq0"

    samtools view -b -q 30 "$bam" > "${base}_mapq30.bam"
    samtools view -b -q 10 "$bam" > "${base}_mapq10.bam"
    samtools view -b -q 1 -U "${base}_mapq0.bam" "$bam" > /dev/null
done