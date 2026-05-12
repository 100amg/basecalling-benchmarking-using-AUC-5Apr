#!/bin/bash
set -euo pipefail

FAST5_DIR="/DATA4/amishi/5Apr/fast5_data"
BASE_OUTPUT_DIR="/DATA4/amishi/5Apr/guppy_comparison"
GUPPY_BIN="/DATA4/amishi/ont-guppy/bin/guppy_basecaller"
BASE_MODEL="/DATA4/amishi/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_hac_prom.cfg"
METHYL_MODEL="/DATA4/amishi/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg"
Q_SCORES=(0 5 10 15 20)

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

    FOLDER_NAME="guppy_hac_basecall_q${Q}"
    OUTPUT_DIR="$BASE_OUTPUT_DIR/$FOLDER_NAME"
    mkdir -p "$OUTPUT_DIR"

    echo "" | tee -a "$LOG"
    echo "-> Basecalling | Q${Q}" | tee -a "$LOG"

    for sample_dir in "$FAST5_DIR"/*/; do
        [ -d "$sample_dir" ] || continue
        sample=$(basename "$sample_dir")
        FINAL_BAM="$OUTPUT_DIR/${sample}_guppy_hac_basecall_q${Q}.bam"

        if [ -f "$FINAL_BAM" ]; then
            echo "  [SKIP] $sample already processed" | tee -a "$LOG"
            continue
        fi

        echo "  Processing: $sample (Q${Q})" | tee -a "$LOG"

        GUPPY_OUT="$OUTPUT_DIR/${sample}_guppy_raw"
        mkdir -p "$GUPPY_OUT"

        "$GUPPY_BIN" \
            -i "$sample_dir" \
            -s "$GUPPY_OUT" \
            -c "$BASE_MODEL" \
            --bam_out \
            --moves_out \
            --min_qscore "$Q" \
            --device cuda:0 \
            2>> "$LOG"

        # Merge pass BAMs
        if ls "$GUPPY_OUT"/pass/*.bam 1>/dev/null 2>&1; then
            samtools merge -f "${FINAL_BAM%.bam}_unsorted.bam" "$GUPPY_OUT"/pass/*.bam
            samtools sort -@ 4 -o "$FINAL_BAM" "${FINAL_BAM%.bam}_unsorted.bam"
            rm -f "${FINAL_BAM%.bam}_unsorted.bam"
            samtools index "$FINAL_BAM"
            READ_COUNT=$(samtools view -c "$FINAL_BAM")
            echo "  $READ_COUNT reads -> $(basename $FINAL_BAM)" | tee -a "$LOG"
        else
            echo "  WARNING: No pass BAMs for $sample at Q${Q}" | tee -a "$LOG"
        fi

        rm -rf "$GUPPY_OUT"

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

    FOLDER_NAME="guppy_hac_methyl_q${Q}"
    OUTPUT_DIR="$BASE_OUTPUT_DIR/$FOLDER_NAME"
    mkdir -p "$OUTPUT_DIR"

    echo "" | tee -a "$LOG"
    echo "-> Methylation calling | Q${Q}" | tee -a "$LOG"

    for sample_dir in "$FAST5_DIR"/*/; do
        [ -d "$sample_dir" ] || continue
        sample=$(basename "$sample_dir")
        FINAL_BAM="$OUTPUT_DIR/${sample}_guppy_hac_methyl_q${Q}.bam"

        if [ -f "$FINAL_BAM" ]; then
            echo "  [SKIP] $sample already processed" | tee -a "$LOG"
            continue
        fi

        echo "  Processing: $sample (Q${Q})" | tee -a "$LOG"

        GUPPY_OUT="$OUTPUT_DIR/${sample}_guppy_raw"
        mkdir -p "$GUPPY_OUT"

        "$GUPPY_BIN" \
            -i "$sample_dir" \
            -s "$GUPPY_OUT" \
            -c "$METHYL_MODEL" \
            --bam_out \
            --moves_out \
            --min_qscore "$Q" \
            --device cuda:0 \
            2>> "$LOG"

        if ls "$GUPPY_OUT"/pass/*.bam 1>/dev/null 2>&1; then
            samtools merge -f "${FINAL_BAM%.bam}_unsorted.bam" "$GUPPY_OUT"/pass/*.bam
            samtools sort -@ 4 -o "$FINAL_BAM" "${FINAL_BAM%.bam}_unsorted.bam"
            rm -f "${FINAL_BAM%.bam}_unsorted.bam"
            samtools index "$FINAL_BAM"
            READ_COUNT=$(samtools view -c "$FINAL_BAM")
            echo "  $READ_COUNT reads -> $(basename $FINAL_BAM)" | tee -a "$LOG"

            MM_CHECK=$(samtools view "$FINAL_BAM" | head -5 | grep -c "MM:Z:" || true)
            if [ "$MM_CHECK" -gt 0 ]; then
                echo "  MM/ML methylation tags confirmed present" | tee -a "$LOG"
            else
                echo "  WARNING: MM/ML tags not found in $sample" | tee -a "$LOG"
            fi
        else
            echo "  WARNING: No pass BAMs for $sample at Q${Q}" | tee -a "$LOG"
        fi

        rm -rf "$GUPPY_OUT"

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
for bam in "$BASE_OUTPUT_DIR"/guppy_hac_basecall_q*/*.bam; do
    [ -f "$bam" ] || continue
    count=$(samtools view -c "$bam")
    printf "  %-10s reads -> %s\n" "$count" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "--- Methylation outputs ---" | tee -a "$LOG"
for bam in "$BASE_OUTPUT_DIR"/guppy_hac_methyl_q*/*.bam; do
    [ -f "$bam" ] || continue
    count=$(samtools view -c "$bam")
    printf "  %-10s reads -> %s\n" "$count" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "All done: $(date)" | tee -a "$LOG"
EOF