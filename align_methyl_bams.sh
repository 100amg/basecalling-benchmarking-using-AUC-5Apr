#!/usr/bin/env bash
# Aligns all unaligned methylation BAMs to a reference, producing new
# sorted+indexed BAMs in a separate output directory.
# Original BAMs are NEVER modified.
# Basecall folders are skipped entirely.
#
# Supports both Dorado and Guppy comparison folders:
#   dorado_fast_comparison/dorado_fast_methyl_q*/
#   dorado_hac_comparison/dorado_hac_methyl_q*/
#   dorado_sup_comparison/dorado_sup_methyl_q*/
#   guppy_comparison/guppy_hac_methyl_q*/
#
# Pipeline per BAM:
#   samtools fastq -T MM,ML  ->  minimap2 -y  ->  samtools sort  ->  .bam + .bai
#
# Usage:
#   bash align_methyl_bams.sh <bam_root_dir> <reference_fasta> <aligned_output_dir> [threads]
#
# Example:
#   bash align_methyl_bams.sh \
#       /Volumes/Amishi_SSD/bio_data/5Apr/new_EpiC_runs \
#       /Volumes/Amishi_SSD/bio_data/5Apr/new_EpiC_runs/reference.fasta \
#       /Volumes/Amishi_SSD/bio_data/5Apr/aligned_methyl \
#       8

set -euo pipefail

# ------------------------------
# Arguments
# ------------------------------
if [ -z "${1:-}" ] || [ -z "${2:-}" ] || [ -z "${3:-}" ]; then
    echo "Usage: $0 <bam_root_dir> <reference_fasta> <aligned_output_dir> [threads]"
    exit 1
fi

BAM_ROOT="$(realpath "$1")"
REF="$(realpath "$2")"
ALIGNED_ROOT="$3"
THREADS="${4:-4}"

[ -d "$BAM_ROOT" ] || { echo "ERROR: BAM root not found: $BAM_ROOT"; exit 1; }
[ -f "$REF"      ] || { echo "ERROR: Reference FASTA not found: $REF"; exit 1; }

for tool in minimap2 samtools; do
    command -v "$tool" &>/dev/null || { echo "ERROR: $tool not found in PATH"; exit 1; }
done

mkdir -p "$ALIGNED_ROOT"
ALIGNED_ROOT="$(realpath "$ALIGNED_ROOT")"

echo "BAM root      : $BAM_ROOT"
echo "Reference     : $REF"
echo "Aligned output: $ALIGNED_ROOT"
echo "Threads       : $THREADS"
echo ""

TOTAL=0
ALIGNED_COUNT=0
SKIPPED=0
FAILED=0

# ------------------------------
# Walk: *_comparison/ -> *_methyl_q*/ -> *.bam
# Matches both dorado_*_comparison and guppy_comparison folders.
# Matches both dorado_*_methyl_q* and guppy_*_methyl_q* subfolders.
# ------------------------------
for COMP_DIR in "$BAM_ROOT"/*_comparison/; do
    [ -d "$COMP_DIR" ] || continue
    COMP_NAME="$(basename "$COMP_DIR")"

    echo "========================================"
    echo "Model folder: $COMP_NAME"
    echo "========================================"

    # Match any *_methyl_q* subfolder (covers dorado_fast_methyl_q0,
    # dorado_hac_methyl_q5, guppy_hac_methyl_q10, etc.)
    # basecall folders are implicitly skipped since they don't match *_methyl_q*
    for QDIR in "$COMP_DIR"*_methyl_q*/; do
        [ -d "$QDIR" ] || continue
        QDIR_NAME="$(basename "$QDIR")"

        OUT_DIR="$ALIGNED_ROOT/$COMP_NAME/$QDIR_NAME"
        mkdir -p "$OUT_DIR"

        echo ""
        echo "  Folder : $QDIR_NAME"
        echo "  Out dir: $OUT_DIR"

        for BAM in "$QDIR"*.bam; do
            [[ "$BAM" == *.bai ]] && continue
            [ -f "$BAM" ]        || continue

            SAMPLE="$(basename "$BAM" .bam)"
            TOTAL=$((TOTAL + 1))

            SORTED_BAM="$OUT_DIR/${SAMPLE}_sorted.bam"

            # Skip if already successfully aligned
            if [ -f "$SORTED_BAM" ] && [ -f "${SORTED_BAM}.bai" ]; then
                SQ=$(samtools view -H "$SORTED_BAM" 2>/dev/null | grep -c "^@SQ" || true)
                if [ "$SQ" -gt 0 ]; then
                    echo "    [SKIP - exists] $SAMPLE"
                    SKIPPED=$((SKIPPED + 1))
                    continue
                fi
            fi

            echo "    [ALIGNING] $SAMPLE"

            if samtools fastq \
                    -T MM,ML \
                    "$BAM" \
                2>/dev/null \
                | minimap2 \
                    -a \
                    -x map-ont \
                    -y \
                    --secondary=no \
                    -t "$THREADS" \
                    "$REF" \
                    - \
                2>/dev/null \
                | samtools sort \
                    -@ "$THREADS" \
                    -o "$SORTED_BAM" \
                    - \
                2>/dev/null; then

                samtools index -@ "$THREADS" "$SORTED_BAM"
                echo "    [DONE]     $SAMPLE  ->  $(basename "$SORTED_BAM")"
                ALIGNED_COUNT=$((ALIGNED_COUNT + 1))

            else
                echo "    [FAILED]   $SAMPLE"
                rm -f "$SORTED_BAM" "${SORTED_BAM}.bai"
                FAILED=$((FAILED + 1))
            fi
        done
    done
done

# ------------------------------
# Summary
# ------------------------------
echo ""
echo "========================================"
echo "  SUMMARY"
echo "========================================"
echo "  Total BAMs found      : $TOTAL"
echo "  Successfully aligned  : $ALIGNED_COUNT"
echo "  Skipped (already done): $SKIPPED"
echo "  Failed                : $FAILED"
echo ""

if [ "$FAILED" -gt 0 ]; then
    echo "  WARNING: $FAILED BAMs failed. Check output above."
    exit 1
fi

echo "All methyl BAMs aligned. Sorted BAMs are in:"
echo "  $ALIGNED_ROOT"
echo ""
echo "Now run the AUC analysis with:"
echo "  python auc_model_qscore_analysis.py $ALIGNED_ROOT <reference_fasta>"