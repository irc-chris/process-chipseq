#!/bin/bash
# ChIP-seq pipeline orchestrator.
# Run directly (do NOT sbatch this script) — it submits each step as its own
# SBATCH job with afterok dependencies so each step has its own log/err/out.
#
# Usage:
#   bash diploid_chipseq_pipeline.sh config.json [--from-step N]
#
# Required config keys:
#   genomeID   - reference genome (hg19, mm9, GRCh38, GMh1, GMh2, GMdip)
#   topDir     - top-level directory; FASTQs must be in topDir/fastq/
#   mode       - "diploid" or "haploid"
#   singleend  - 0 for paired-end (default), 1 for single-end
#   vcf        - path to phased VCF (required for haploid mode)
#   bedfile1   - BED file for hap1 counts
#   bedfile2   - BED file for hap2 counts
#
# Steps:
#   1  align               (bwa mem -> output.sorted.bam)
#   2  split/haplotag      (diploid: split by chrom; haploid: whatshap + split by HP tag)
#   3  sort_index          (hap1.bam + hap2.bam -> hap1.sorted.bam + hap2.sorted.bam)
#   4  filter_mapq         (output.sorted.bam -> output_mapq{30,10,0}.bam) [parallel with step 3]
#   5  counts              (bedtools multicov on hap1.sorted.bam and hap2.sorted.bam)

set -euo pipefail

# ---- Argument parsing ----
configFile=""
from_step=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from-step) from_step="$2"; shift 2 ;;
        --from-step=*) from_step="${1#*=}"; shift ;;
        *) configFile="$1"; shift ;;
    esac
done

configFile=$(realpath "${configFile:?"Usage: $0 config.json [--from-step N]"}")
[ -f "$configFile" ] || { echo "Config file not found: $configFile"; exit 1; }

[[ "$from_step" =~ ^[1-5]$ ]] || { echo "--from-step must be 1-5"; exit 1; }

command -v jq     >/dev/null 2>&1 || { echo "jq is required"; exit 1; }
command -v sbatch >/dev/null 2>&1 || { echo "sbatch not found — are you on the cluster?"; exit 1; }

SCRIPT_DIR=$(dirname "$(realpath "$0")")

topDir=$(jq -r '.topDir'          "$configFile")
mode=$(jq -r   '.mode // "haploid"' "$configFile")
name=${topDir##*/}
tmpdir="$topDir/tmp"
outputdir="$topDir/aligned"

[[ "$mode" == "diploid" || "$mode" == "haploid" ]] \
    || { echo "mode must be \"diploid\" or \"haploid\""; exit 1; }

# ---- File presence checks per step ----
# Validates that required inputs exist before submitting from that step.

check_files() {
    local step=$1; shift
    local all_ok=true
    for f in "$@"; do
        # Support glob patterns (e.g. for fastq dir)
        if compgen -G "$f" > /dev/null 2>&1; then
            echo "  ✅  $f"
        else
            echo "  ❌  $f  (not found)"
            all_ok=false
        fi
    done
    $all_ok || { echo "Missing required inputs for step $step. Aborting."; exit 1; }
}

echo "=== ChIP-seq Pipeline ==="
echo "  config:     $configFile"
echo "  mode:       $mode"
echo "  name:       $name"
echo "  from step:  $from_step"
echo ""

if [ "$from_step" -gt 1 ]; then
    echo "Checking required files for step $from_step..."
    case "$from_step" in
        2)
            check_files 2 "$tmpdir/output.sorted.bam"
            ;;
        3)
            check_files 3 "$outputdir/hap1.bam" "$outputdir/hap2.bam"
            ;;
        4)
            check_files 4 "$tmpdir/output.sorted.bam"
            ;;
        5)
            check_files 5 "$outputdir/hap1.sorted.bam" "$outputdir/hap2.sorted.bam"
            ;;
    esac
    echo ""
fi

# ---- Job submission ----
# PREV_JOB tracks the dependency for step 3 (set by step 1 or 2).
# PREV_JOB4 tracks the dependency for step 4 (same as PREV_JOB since step 4
# runs in parallel with step 3, both depending on step 2).
PREV_JOB=""

submit() {
    # submit <job-name-suffix> <dep_job_id_or_empty> <script> [args...]
    local suffix=$1; local dep=$2; local script=$3; shift 3
    local dep_flag=""
    [ -n "$dep" ] && dep_flag="--dependency=afterok:$dep"
    sbatch --parsable \
        --job-name="${name}_${suffix}" \
        $dep_flag \
        "$SCRIPT_DIR/$script" "$@"
}

# Step 1
if [ "$from_step" -le 1 ]; then
    JOB1=$(submit align "" step1_align.sh "$configFile")
    echo "Submitted STEP 1 (align):         job $JOB1  ->  ${name}_align_${JOB1}.{out,err}"
    PREV_JOB=$JOB1
fi

# Step 2 — dispatch by mode
if [ "$from_step" -le 2 ]; then
    if [ "$mode" = "diploid" ]; then
        JOB2=$(submit split "$PREV_JOB" step2_split_haplotype.sh "$configFile")
        echo "Submitted STEP 2 (split haplos):  job $JOB2  ->  ${name}_split_${JOB2}.{out,err}"
    else
        JOB2=$(submit whatshap "$PREV_JOB" step2_whatshap.sh "$configFile")
        echo "Submitted STEP 2 (whatshap):      job $JOB2  ->  ${name}_whatshap_${JOB2}.{out,err}"
    fi
    PREV_JOB=$JOB2
fi

# Steps 3 and 4 run in parallel (both depend on step 2 / PREV_JOB)
if [ "$from_step" -le 3 ]; then
    JOB3=$(submit sort "$PREV_JOB" step3_sort_index.sh "$configFile")
    echo "Submitted STEP 3 (sort/index):    job $JOB3  ->  ${name}_sort_${JOB3}.{out,err}"
fi

if [ "$from_step" -le 4 ]; then
    JOB4=$(submit filter "$PREV_JOB" step4_filter_mapq.sh "$configFile")
    echo "Submitted STEP 4 (MAPQ filter):   job $JOB4  ->  ${name}_filter_${JOB4}.{out,err}"
fi

# Step 5 depends only on step 3 (needs sorted hap BAMs)
if [ "$from_step" -le 5 ]; then
    JOB5=$(submit counts "$JOB3" step5_counts.sh "$configFile")
    echo "Submitted STEP 5 (counts):        job $JOB5  ->  ${name}_counts_${JOB5}.{out,err}"
fi

echo ""
echo "Monitor with: squeue -u $USER"
