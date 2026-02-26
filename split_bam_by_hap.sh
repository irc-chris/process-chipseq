#!/bin/bash
#SBATCH --job-name=process_bams
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=108:00:00
#SBATCH --mem=20G
#SBATCH --partition=weka2

INPUT=$1
HAP1=${1%.bam}_hap1.bam
HAP2=${1%.bam}_hap2.bam
HAP1_r=${1%.bam}_hap1_renamed.bam
HAP2_r=${1%.bam}_hap2_renamed.bam

echo "Splitting by haplotype..."
{ samtools view -H "$INPUT"; samtools view "$INPUT" | awk '$3 ~ /_h1/'; } \
    | samtools view -bS -o "$HAP1"

{ samtools view -H "$INPUT"; samtools view "$INPUT" | awk '$3 ~ /_h2/'; } \
    | samtools view -bS -o "$HAP2"

echo "Removing haplotag..."
samtools view -h "$HAP1" | sed 's/_h1//g' | samtools view -bS - > "$HAP1_r"
samtools view -h "$HAP2" | sed 's/_h2//g' | samtools view -bS - > "$HAP2_r"
