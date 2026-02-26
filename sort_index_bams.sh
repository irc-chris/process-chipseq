#!/bin/bash
#SBATCH --job-name=process_bams
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=108:00:00
#SBATCH --mem=20G
#SBATCH --partition=weka2

INPUT=$1
HAP1_r=${1%.bam}_hap1_renamed.bam
HAP2_r=${1%.bam}_hap2_renamed.bam
HAP1_s=${1%.bam}_hap1_renamed.sorted.bam
HAP2_s=${1%.bam}_hap2_renamed.sorted.bam

echo "Sorting..."
samtools sort -o "$HAP1_s" "$HAP1_r"
samtools sort -o "$HAP2_s" "$HAP2_r"

echo "Indexing..."
samtools index "$HAP1_s"
samtools index "$HAP2_s"