#!/bin/bash
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --partition=weka2
#SBATCH --cpus-per-task=4

export PATH=/mnt/altnas/work/ishawnia/empvsag/process-chipseq/samtools-1.19:$PATH

set -euo pipefail

configFile=${1:?"Usage: $0 config.json"}
[ -f "$configFile" ] || { echo "Config not found: $configFile"; exit 1; }
command -v whatshap >/dev/null 2>&1 || { echo "whatshap not in PATH"; exit 1; }

topDir=$(jq -r '.topDir'     "$configFile")
genomeID=$(jq -r '.genomeID' "$configFile")
vcf=$(jq -r '.vcf'           "$configFile")

[ -f "$vcf" ] || { echo "VCF not found: $vcf"; exit 1; }

case $genomeID in
    hg19)   refSeq="/gpfs0/juicer2/references/hg19/Homo_sapiens_assembly19.fasta" ;;
    mm9)    refSeq="/gpfs0/juicer2/references/Mus_musculus_assembly9_norandom.fasta" ;;
    GRCh38) refSeq="/gpfs0/juicer2/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" ;;
    GMh1)   refSeq="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs/hap1_corrected.fa" ;;
    GMh2)   refSeq="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs/hap2_corrected.fa" ;;
    *)  echo "Unknown genomeID: $genomeID"; exit 1 ;;
esac

[ -e "$refSeq" ] || { echo "Reference not found: $refSeq"; exit 1; }

tmpdir="$topDir/tmp"
outputdir="$topDir/aligned"
ALIGNED_BAM="$tmpdir/output.sorted.bam"

[ -f "$ALIGNED_BAM" ] || { echo "Aligned BAM not found: $ALIGNED_BAM"; exit 1; }

HAPLOTAGGED_BAM="$tmpdir/output.haplotagged.bam"
HAP1_OUT="$outputdir/hap1.bam"
HAP2_OUT="$outputdir/hap2.bam"

# WhatsHap requires a bgzipped, tabix-indexed VCF.
# If the provided VCF is uncompressed, bgzip/index a copy in tmpdir.
if [[ "$vcf" != *.gz ]]; then
    echo "$(date) bgzipping and indexing VCF..."
    bgzip -c "$vcf" > "$tmpdir/phased.vcf.gz"
    tabix -p vcf "$tmpdir/phased.vcf.gz"
    vcf="$tmpdir/phased.vcf.gz"
elif [ ! -f "${vcf}.tbi" ]; then
    echo "$(date) Copying and indexing VCF..."
    cp "$vcf" "$tmpdir/phased.vcf.gz"
    tabix -p vcf "$tmpdir/phased.vcf.gz"
    vcf="$tmpdir/phased.vcf.gz"
fi

echo "$(date) Running WhatsHap haplotag..."
whatshap haplotag \
    --output "$HAPLOTAGGED_BAM" \
    --reference "$refSeq" \
    "$vcf" \
    "$ALIGNED_BAM"

echo "$(date) Indexing haplotagged BAM..."
samtools index "$HAPLOTAGGED_BAM"

echo "$(date) Splitting by HP tag..."
samtools view -@ 2 -b -d HP:1 "$HAPLOTAGGED_BAM" > "$HAP1_OUT"
samtools view -@ 2 -b -d HP:2 "$HAPLOTAGGED_BAM" > "$HAP2_OUT"

echo "$(date) WhatsHap haplotagging and splitting complete."
