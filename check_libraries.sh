# Tools
command -v bwa >/dev/null 2>&1 || { echo "bwa not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not found"; exit 1; }
command -v java >/dev/null 2>&1 || { echo "java not found"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found"; exit 1; }

# Picard jar
[ -f "$picardloc" ] || { echo "Picard jar not found at $picardloc"; exit 1; }

# Reference genome
[ -f "$refSeq" ] || { echo "Reference $refSeq not found"; exit 1; }

# BWA index files
for ext in amb ann bwt pac sa; do
    [ -f "${refSeq}.${ext}" ] || { echo "BWA index missing: ${refSeq}.${ext}"; exit 1; }
done

# Chromosome lengths file
[ -f "$chromosomelengths" ] || { echo "Chromosome lengths file not found: $chromosomelengths"; exit 1; }

# Fastq directory and files
[ -d "$topDir/fastq" ] || { echo "Fastq directory not found"; exit 1; }

# BED file for counting
[ -f "$bedfile" ] || { echo "BED file $bedfile not found"; exit 1; }

# Output directory doesn't already exist
[ ! -d "$outputdir" ] || { echo "Output directory $outputdir already exists"; exit 1; }

# No leftover tmp files
[ ! -f "$topDir/fastq/tmp_R1.fastq" ] || { echo "Leftover tmp_R1.fastq found, remove it first"; exit 1; }
[ ! -f "$topDir/fastq/tmp_R2.fastq" ] || { echo "Leftover tmp_R2.fastq found, remove it first"; exit 1; }

echo "All checks passed."