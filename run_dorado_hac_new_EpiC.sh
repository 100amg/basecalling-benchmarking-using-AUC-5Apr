#!/bin/bash
set -euo pipefail

# ─────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────

POD5_DIR="/Volumes/Amishi_SSD/bio_data/5Apr/new_EpiC_runs/pod5_data"
BASE_OUTPUT_DIR="/Volumes/Amishi_SSD/bio_data/5Apr/new_EpiC_runs/dorado_hac_comparison"

DORADO_BIN="/Volumes/Amishi_SSD/bio_data/6Feb/models_and_libraries/dorado-1.3.1-osx-arm64/bin/dorado"

# Models
BASE_MODEL="/Volumes/Amishi_SSD/bio_data/6Feb/models_and_libraries/dorado-1.3.1-osx-arm64/dorado_models/dna_r10.4.1_e8.2_400bps_hac@v5.2.0"
METHYL_MODEL="/Volumes/Amishi_SSD/bio_data/6Feb/models_and_libraries/dorado-1.3.1-osx-arm64/dorado_models/dna_r10.4.1_e8.2_400bps_hac@v5.2.0_5mCG_5hmCG@v1"

Q_SCORES=(5 7 11 15 18)

# ─────────────────────────────────────────────
# SETUP
# ─────────────────────────────────────────────

mkdir -p "$BASE_OUTPUT_DIR"
LOG="$BASE_OUTPUT_DIR/run_log.txt"
echo "Run started: $(date)" | tee "$LOG"

# ─────────────────────────────────────────────
# STEP 1: BASECALLING ONLY
# ─────────────────────────────────────────────

echo "" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"
echo "  STEP 1: BASECALLING ONLY" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"

for Q in "${Q_SCORES[@]}"; do

    FOLDER_NAME="dorado_hac_basecall_q${Q}"
    OUTPUT_DIR="$BASE_OUTPUT_DIR/$FOLDER_NAME"
    mkdir -p "$OUTPUT_DIR"

    echo "" | tee -a "$LOG"
    echo "-> Basecalling | Q${Q}" | tee -a "$LOG"

    for pod5 in "$POD5_DIR"/*.pod5; do
        [ -e "$pod5" ] || continue
        name=$(basename "$pod5" .pod5)
        FINAL_BAM="$OUTPUT_DIR/${name}_hac_basecall_q${Q}.bam"
        TEMP_DIR="$OUTPUT_DIR/${name}_tempdir"

        if [ -f "$FINAL_BAM" ]; then
            echo "  [SKIP] $name already processed" | tee -a "$LOG"
            continue
        fi

        echo "  Processing: $name" | tee -a "$LOG"

        mkdir -p "$TEMP_DIR"

        "$DORADO_BIN" basecaller \
            "$BASE_MODEL" \
            "$pod5" \
            --min-qscore "$Q" \
            --emit-moves \
            --output-dir "$TEMP_DIR" \
            2>> "$LOG"

        samtools merge -f "${FINAL_BAM%.bam}_unsorted.bam" "$TEMP_DIR"/*.bam
        samtools sort -@ 4 -o "$FINAL_BAM" "${FINAL_BAM%.bam}_unsorted.bam"
        rm -f "${FINAL_BAM%.bam}_unsorted.bam"
        rm -rf "$TEMP_DIR"

        samtools index "$FINAL_BAM"

        READ_COUNT=$(samtools view -c "$FINAL_BAM")
        echo "  $READ_COUNT reads -> $(basename $FINAL_BAM)" | tee -a "$LOG"

    done
done

# ─────────────────────────────────────────────
# STEP 2: METHYLATION CALLING
# ─────────────────────────────────────────────

echo "" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"
echo "  STEP 2: METHYLATION CALLING" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"

for Q in "${Q_SCORES[@]}"; do

    FOLDER_NAME="dorado_hac_methyl_q${Q}"
    OUTPUT_DIR="$BASE_OUTPUT_DIR/$FOLDER_NAME"
    mkdir -p "$OUTPUT_DIR"

    echo "" | tee -a "$LOG"
    echo "-> Methylation calling | Q${Q}" | tee -a "$LOG"

    for pod5 in "$POD5_DIR"/*.pod5; do
        [ -e "$pod5" ] || continue
        name=$(basename "$pod5" .pod5)
        FINAL_BAM="$OUTPUT_DIR/${name}_hac_methyl_q${Q}.bam"
        TEMP_DIR="$OUTPUT_DIR/${name}_tempdir"

        if [ -f "$FINAL_BAM" ]; then
            echo "  [SKIP] $name already processed" | tee -a "$LOG"
            continue
        fi

        echo "  Processing: $name" | tee -a "$LOG"

        mkdir -p "$TEMP_DIR"

        "$DORADO_BIN" basecaller \
            "$BASE_MODEL" \
            "$pod5" \
            --modified-bases-models "$METHYL_MODEL" \
            --min-qscore "$Q" \
            --emit-moves \
            --output-dir "$TEMP_DIR" \
            2>> "$LOG"

        samtools merge -f "${FINAL_BAM%.bam}_unsorted.bam" "$TEMP_DIR"/*.bam
        samtools sort -@ 4 -o "$FINAL_BAM" "${FINAL_BAM%.bam}_unsorted.bam"
        rm -f "${FINAL_BAM%.bam}_unsorted.bam"
        rm -rf "$TEMP_DIR"

        samtools index "$FINAL_BAM"

        READ_COUNT=$(samtools view -c "$FINAL_BAM")
        echo "  $READ_COUNT reads -> $(basename $FINAL_BAM)" | tee -a "$LOG"

        MM_CHECK=$(samtools view "$FINAL_BAM" | head -5 | grep -c "MM:Z:" || true)
        if [ "$MM_CHECK" -gt 0 ]; then
            echo "  MM/ML methylation tags confirmed present" | tee -a "$LOG"
        else
            echo "  WARNING: MM/ML tags not found in $name" | tee -a "$LOG"
        fi

    done
done

# ─────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────

echo "" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"
echo "  FINAL SUMMARY" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "--- Basecall outputs ---" | tee -a "$LOG"
for bam in "$BASE_OUTPUT_DIR"/dorado_hac_basecall_q*/*.bam; do
    [ -f "$bam" ] || continue
    count=$(samtools view -c "$bam")
    printf "  %-10s reads -> %s\n" "$count" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "--- Methylation outputs ---" | tee -a "$LOG"
for bam in "$BASE_OUTPUT_DIR"/dorado_hac_methyl_q*/*.bam; do
    [ -f "$bam" ] || continue
    count=$(samtools view -c "$bam")
    printf "  %-10s reads -> %s\n" "$count" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "All done: $(date)" | tee -a "$LOG"
