#!/bin/bash
#SBATCH --job-name=process_chipseq
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=108:00:00
#SBATCH --mem=30G
#SBATCH --partition=weka
#SBATCH --cpus-per-task=21
export PATH=/mnt/altnas/work/ishawnia/empvsag/process-chipseq/samtools-1.19:$PATH

set -euo pipefail

# This script takes in fastqs and outputs processed ChIP-Seq
# data (signal tracks, peak calls, etc.). 
# 
#
# BWA-mem is used for alignment and MACS2 for peakcalling/signal
# track generation.
#
#

# ---- Config parsing ----
# Usage: provide a JSON config file as first argument
#   sbatch my_process_chipseq.sh config.json

configFile=${1:?"Please provide a JSON config file as argument"}
[ -f "$configFile" ] || { echo "Config file $configFile not found"; exit 1; }
command -v jq >/dev/null 2>&1 || { echo "jq not found, needed for JSON parsing"; exit 1; }

genomeID=$(jq -r '.genomeID' "$configFile")
topDir=$(jq -r '.topDir' "$configFile")
singleend=$(jq -r '.singleend // 0' "$configFile")
picardloc=$(jq -r '.picardloc // "/gpfs0/work/suhas/scripts/picard_tools/picard/dist/picard.jar"' "$configFile")
controlonly=$(jq -r '.controlonly // 0' "$configFile")
inputcontrol=$(jq -r '.inputcontrol // 0' "$configFile")
name=${topDir##*/}

echo $genomeID 

## Read arguments
usageHelp="Usage: ${0##*/} -g genomeID [-d topDir] [-p picard-tools location] [-s] [-h] [-i input_control] [-c]"
genomeHelp="   genomeID must be one of \"mm9\" (mouse), \"hg18\" \"hg19\" (human), \"sCerS288c\" (yeast), \"GMh1\" (human diploid maternal), \"GMh2\" (human diploid paternal), \"Pf3D7\" (Plasmodium falciparum aka malaria), \"canFam3 (dog)\", or \"dMel\" (fly)"
dirHelp="   [topDir] is the top level directory (default \"$topDir\")\n     [topDir]/fastq must contain the fastq files\n     [topDir]/splits will be created to contain the temporary split files\n     [topDir]/aligned will be created for the final alignment"
picardHelp="	[picard-tools location] is the path for the picard tools jar file. Default is /aidenlab/work/suhas/scripts/picard_tools/picard/dist/picard.jar"
inputcontrolHelp="	[input_control} is a dedup-ed bam file for the input control for the experiment"
inputHelp="	-c: this option indicates that you are only aligning and dedup-ing an input control experiment"
readtypeHelp="   -s: single-end reads (default paired-end)"
helpHelp="   -h: print this help and exit"

printHelpAndExit() {
    echo "$usageHelp"
    echo "$genomeHelp"
    echo -e "$dirHelp"
    echo "$picardHelp"
    echo "$readtypeHelp"
    echo "$inputcontrolHelp"
    echo "$inputHelp"
    echo "$helpHelp"
    exit $1
}

while getopts "d:g:p:i:shc" opt; do
    case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
	d) topDir=$OPTARG ;;
	s) singleend=1 ;;
	p) picardloc=$OPTARG ;;
	i) inputcontrol=$OPTARG ;;
	c) controlonly=1  ;;
	[?]) printHelpAndExit 1;;
    esac
done


## Set reference sequence based on genome ID
case $genomeID in
    hg19) refSeq="/gpfs0/juicer2/references/hg19/Homo_sapiens_assembly19.fasta"
	  chromosomelengths="/gpfs0/work/suhas/scripts/chipseq/hg19.len";;
    mm9) refSeq="/gpfs0/juicer2/references/Mus_musculus_assembly9_norandom.fasta"
	 chromosomelengths="/gpfs0/work/suhas/scripts/chipseq/mm9.len";;
    GRCh38) refSeq="/gpfs0/juicer2/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
	    chromosomelengths="/gpfs0/work/suhas/scripts/chipseq/GRCh38.len";;
	GMh1) refSeq="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs/hap1_corrected.fa"
         chromosomelengths="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs/hap1_corrected.len";;
    GMh2) refSeq="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs/hap2_corrected.fa"
         chromosomelengths="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs/hap2_corrected.len";;
    GMdip) refSeq="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs/diploid.fa"
         chromosomelengths="/mnt/altnas/work/ishawnia/empvsag/process-chipseq/refs/diploid.len";;
    *)  echo "$usageHelp"
        echo "$genomeHelp"
        exit 1;;
esac


## Check that refSeq exists 
if [ ! -e "$refSeq" ]; then
    echo "Reference sequence $refSeq does not exist";
    exit 1;
fi

## Hard-coded directory here
if [ $singleend == 0 ]; then
	fastqdir=$topDir"/fastq/*_R*.fastq"
else
	fastqdir=$topDir"/fastq/*.fastq"
fi
outputdir=$topDir"/aligned"
tmpdir=$topDir"/tmp"

############### ISHAWNIA's VARIABLES ###################
bedfile1="anchors_h1.bed"
bedfile2="anchors_h2.bed"

#checking libraries
# . ./check_libraries.sh

echo "Arguments:
  Genome ID: $genomeID
  Top level directory: $topDir
  Picard tools location: $picardloc
  Single-end reads: $singleend
  Control-only experiment: $controlonly
  Input control file: $inputcontrol
  Fastq directory: $fastqdir
  Output directory: $outputdir
  Tmp directory: $tmpdir
  Reference sequence: $refSeq
  Chromosome lengths file: $chromosomelengths
  BED file for counting: $bedfile1
  BED file for counting: $bedfile2
"

## Check that fastq directory exists and has proper fastq files
if [ ! -d "$topDir/fastq" ]; then
		echo "Directory \"$topDir/fastq\" does not exist."
		echo "Create \"$topDir/$fastq\" and put fastq files to be aligned there."
		echo "Type \"align.sh -h \" for help"
		exit 1
else 
		if stat -t ${fastqdir} >/dev/null 2>&1
				then
				echo "Looking for fastq files...fastq files exist"
		else
			echo "Failed to find any files matching ${fastqdir}"
			echo "Type \"align.sh -h \" for help"
			exit 1			
		fi
fi

## Create output and tmp directories
# if [ -d "$outputdir" ]; then
#     echo "Move or remove directory \"$outputdir\" before proceeding."
# 		echo "Type \"align.sh -h \" for help"
# 		exit 1			
# fi

mkdir -p $outputdir
mkdir -p $tmpdir

echo -e "Aligning files matching $fastqdir to genome $genomeID"


if [ $singleend == 0 ]; then

	Combining all R1s and R2s into one file
	echo "Combining R1s and R2s into one file"
	cat $topDir"/fastq/"*"R1"*".fastq" > $tmpdir/.tmp_R1.fastq
	cat $topDir"/fastq/"*"R2"*".fastq" > $tmpdir/.tmp_R2.fastq

	Align with bwa-mem
	-M flag used to mark smaller split alignments as secondary for Picard compatibility
	bwa mem -t 20 -SP5M $refSeq \
		$tmpdir/.tmp_R1.fastq \
		$tmpdir/.tmp_R2.fastq \
		| samtools sort -o $tmpdir/output.sorted.bam
	
	echo "BWA alignment and  sorting completed. Now trying to mark duplicates with Picard."
	# java -jar $picardloc MarkDuplicates \
	# 	I=$tmpdir/output.sorted.bam \
	# 	O=$tmpdir/output.nodups.bam \
	# 	METRICS_FILE=$tmpdir/dedup_metrics.txt \
	# 	REMOVE_DUPLICATES=TRUE
	# samtools markdup -r $tmpdir/output.sorted.bam $tmpdir/output.nodups.bam

	echo "Deduplication completed. Trying to index/resort"
	# samtools sort -@ 12 -o $topDir/aligned/aligned1.sorted.bam $tmpdir/output.nodups.bam
	samtools sort -o hap1_renamed.sorted.bam hap1_renamed.bam
	samtools index hap1_renamed.sorted.bam

	echo "Indexing final bam file"
	samtools index $topDir/aligned/aligned1.sorted.bam
	
	echo "Getting counts"
	bedtools multicov -bams $topDir/aligned/aligned1.sorted.bam -bed $bedfile1 > ${topDir}/1counts1.tsv
	echo "Counts written to ${topDir}/1counts1.tsv"
	bedtools multicov -bams $topDir/aligned/aligned1.sorted.bam -bed $bedfile2 > ${topDir}/1counts2.tsv
	echo "Counts written to ${topDir}/1counts2.tsv"


#### not gonna use bc it's paired end
else

	## Combining all R1s into one file
	cat $topDir"/fastq/"*".fastq" > $tmpdir/.tmp_R1.fastq

	## Align with bwa-mem
	bwa mem -t 20 -SP5M $refSeq \
		$tmpdir/.tmp_R1.fastq \
		| samtools view -bS - > $tmpdir/output.bam

fi