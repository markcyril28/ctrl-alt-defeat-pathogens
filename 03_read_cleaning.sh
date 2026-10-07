#!/usr/bin/env bash
# ============================================================================
# Program 03: Read cleaning — fastp and fastplong
# ============================================================================
# Trims adapters and low-quality ends, and drops reads too short to help an
# assembler. fastp does the paired short reads, fastplong the long ones.
#
# Every setting below throws some data away on purpose. These are the values
# in the trainee guide, and changing one changes your assembly, so they are
# worth understanding rather than copying. Worksheet section 4 asks you to
# explain each one.
#
# The cleaned reads go to a new folder. The raw files are never written to:
# they are the one thing in this workshop you cannot regenerate.
#
# Needs:  conda activate metagenomics_env
# Run from the workshop folder:  bash 03_read_cleaning.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
THREADS=12

# fastp, for the short reads
CUT_RIGHT_WINDOW_SIZE=4          # width of the sliding window, in bases
CUT_RIGHT_MEAN_QUALITY=20        # cut when the window mean Phred drops below this
QUALIFIED_QUALITY_PHRED=20       # a base at or above this counts as "qualified"
UNQUALIFIED_PERCENT_LIMIT=30     # drop the read if more than this % is unqualified
LENGTH_REQUIRED=50               # drop reads shorter than this after trimming

# fastplong, for the long reads. Gentler on purpose: long-read quality is
# lower per base by nature, and hard trimming throws away the length that
# made the read useful in the first place.
LONG_LENGTH_REQUIRED=1000

DATA="Datasets/Metagenome Assembly and Annotation"
OUT="RESULTS/metagenomics/03_Read_Cleaning"

SHORT_R1="$DATA/raw_data/short_reads_R1.fastq"
SHORT_R2="$DATA/raw_data/short_reads_R2.fastq"
LONG_READS="$DATA/raw_data/long_reads.fastq"
# ============================================================

mkdir -p "$OUT"

# --cut_right slides a window along the read and cuts where quality falls off,
# which is why the window size and the threshold travel with it.
fastp \
    -i "$SHORT_R1" \
    -I "$SHORT_R2" \
    -o "$OUT/short_R1.clean.fastq.gz" \
    -O "$OUT/short_R2.clean.fastq.gz" \
    --cut_right \
    --cut_right_window_size "$CUT_RIGHT_WINDOW_SIZE" \
    --cut_right_mean_quality "$CUT_RIGHT_MEAN_QUALITY" \
    --qualified_quality_phred "$QUALIFIED_QUALITY_PHRED" \
    --unqualified_percent_limit "$UNQUALIFIED_PERCENT_LIMIT" \
    --length_required "$LENGTH_REQUIRED" \
    --html "$OUT/fastp.html" \
    --json "$OUT/fastp.json" \
    --thread "$THREADS"

fastplong \
    -i "$LONG_READS" \
    -o "$OUT/long_reads.clean.fastq.gz" \
    --length_required "$LONG_LENGTH_REQUIRED" \
    --html "$OUT/fastplong.html" \
    --json "$OUT/fastplong.json" \
    --thread "$THREADS"
