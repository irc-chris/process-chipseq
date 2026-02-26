#!/bin/bash
#SBATCH --job-name=filter-bams
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=108:00:00
#SBATCH --mem=20G
#SBATCH --partition=weka2

samtools view -b -q 30 $1 > ${1}-filtered.bam
