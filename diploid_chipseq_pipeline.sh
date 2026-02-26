#!/bin/bash
# ChIP-seq pipeline orchestrator.
# Run directly (do NOT sbatch this script) — it submits each step as its own
# SBATCH job with afterok dependencies so each step has its own log/err/out.
#
# Usage: bash diploid_chipseq_pipeline.sh config.json
#
# Required config keys:
#   genomeID   - reference genome (hg19, mm9, GRCh38, GMh1, GMh2, GMdip)
#   topDir     - top-level directory; FASTQs must be in topDir/fastq/
#   mode       - "diploid" (align + split by haplotype) or "haploid" (align only)
#   singleend  - 0 for paired-end (default), 1 for single-end

set -euo pipefail

configFile=$(realpath "${1:?"Usage: $0 config.json"}")
[ -f "$configFile" ] || { echo "Config file not found: $configFile"; exit 1; }
command -v jq     >/dev/null 2>&1 || { echo "jq is required"; exit 1; }
command -v sbatch >/dev/null 2>&1 || { echo "sbatch not found — are you on the cluster?"; exit 1; }

SCRIPT_DIR=$(dirname "$(realpath "$0")")

topDir=$(jq -r '.topDir' "$configFile")
mode=$(jq -r '.mode // "haploid"' "$configFile")
name=${topDir##*/}

[[ "$mode" == "diploid" || "$mode" == "haploid" ]] \
    || { echo "mode must be \"diploid\" or \"haploid\""; exit 1; }

echo "=== ChIP-seq Pipeline Submission ==="
echo "  config:  $configFile"
echo "  topDir:  $topDir"
echo "  mode:    $mode"
echo "  name:    $name"
echo ""

# Step 1: Align
JOB1=$(sbatch --parsable \
    --job-name="${name}_align" \
    "$SCRIPT_DIR/step1_align.sh" "$configFile")
echo "Submitted STEP 1 (align):         job $JOB1  ->  ${name}_align_${JOB1}.{out,err}"

PREV_JOB=$JOB1

if [ "$mode" = "diploid" ]; then
    # Step 2: Split by haplotype and strip tags (diploid only)
    JOB2=$(sbatch --parsable \
        --job-name="${name}_split" \
        --dependency=afterok:$JOB1 \
        "$SCRIPT_DIR/step2_split_haplotype.sh" "$configFile")
    echo "Submitted STEP 2 (split haplos):  job $JOB2  ->  ${name}_split_${JOB2}.{out,err}"
    PREV_JOB=$JOB2
fi

# Step 3: Sort and index (1 BAM for haploid, 2 for diploid)
JOB3=$(sbatch --parsable \
    --job-name="${name}_sort" \
    --dependency=afterok:$PREV_JOB \
    "$SCRIPT_DIR/step3_sort_index.sh" "$configFile")
echo "Submitted STEP 3 (sort/index):    job $JOB3  ->  ${name}_sort_${JOB3}.{out,err}"

# Step 4: Filter by MAPQ (1 BAM for haploid, 2 for diploid)
JOB4=$(sbatch --parsable \
    --job-name="${name}_filter" \
    --dependency=afterok:$JOB3 \
    "$SCRIPT_DIR/step4_filter_mapq.sh" "$configFile")
echo "Submitted STEP 4 (MAPQ filter):   job $JOB4  ->  ${name}_filter_${JOB4}.{out,err}"

echo ""
echo "Monitor with: squeue -u $USER"
