#!/bin/bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=30G
#SBATCH --partition=weka
#SBATCH --cpus-per-task=21

export PATH=/mnt/altnas/work/ishawnia/empvsag/process-chipseq/samtools-1.19:$PATH

set -euo pipefail

configFile=${1:?"Usage: $0 config.json"}
[ -f "$configFile" ] || { echo "Config not found: $configFile"; exit 1; }

genomeID=$(jq -r '.genomeID'        "$configFile")
topDir=$(jq -r   '.topDir'          "$configFile")
singleend=$(jq -r '.singleend // 0' "$configFile")

case $genomeID in
    hg19)   refSeq="/gpfs0/juicer2/references/hg19/Homo_sapiens_assembly19.fasta" ;;
    mm9)    refSeq="/gpfs0/juicer2/references/Mus_musculus_assembly9_norandom.fasta" ;;
    GRCh38) refSeq="/gpfs0/juicer2/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" ;;
    GMh1)   refSeq="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs/hap1_corrected.fa" ;;
    GMh2)   refSeq="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs/hap2_corrected.fa" ;;
    GMdip)  refSeq="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs/diploid.fa" ;;
    *)  echo "Unknown genomeID: $genomeID"; exit 1 ;;
esac

[ -e "$refSeq" ] || { echo "Reference not found: $refSeq"; exit 1; }

tmpdir="$topDir/tmp"
mkdir -p "$topDir/aligned" "$tmpdir"

ALIGNED_BAM="$tmpdir/output.sorted.bam"

echo "$(date) Aligning to $genomeID ($refSeq)..."

if [ "$singleend" -eq 0 ]; then
    echo "Combining paired-end FASTQs..."
    cat "$topDir/fastq/"*R1*.fastq > "$tmpdir/.tmp_R1.fastq"
    cat "$topDir/fastq/"*R2*.fastq > "$tmpdir/.tmp_R2.fastq"
    bwa mem -t 20 -SP5M "$refSeq" "$tmpdir/.tmp_R1.fastq" "$tmpdir/.tmp_R2.fastq" \
        | samtools sort -@ 4 -o "$ALIGNED_BAM"
else
    echo "Combining single-end FASTQs..."
    cat "$topDir/fastq/"*.fastq > "$tmpdir/.tmp_R1.fastq"
    bwa mem -t 20 -SP5M "$refSeq" "$tmpdir/.tmp_R1.fastq" \
        | samtools sort -@ 4 -o "$ALIGNED_BAM"
fi

echo "$(date) Alignment complete: $ALIGNED_BAM"
