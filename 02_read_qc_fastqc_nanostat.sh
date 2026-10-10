#!/usr/bin/env bash
# ============================================================================
# Program 02: Read quality control — FastQC and NanoStat
# ============================================================================
# FastQC reports per-base quality, adapter content and duplication for the
# short reads. NanoStat reports the long reads, where the numbers that matter
# are read length N50 and median read quality.
#
# Run this twice: once now on the raw reads, and again after
# 03_read_cleaning_fastp_fastplong.sh with WHICH="cleaned". Cleaning changes your data; it
# does not improve it by definition. Measuring both is the only way to know
# what it did. Worksheet sections 3 and 5.
#
# Needs:  conda activate meta_env
# Run from the workshop folder:  bash 02_read_qc_fastqc_nanostat.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
WHICH="raw"       # raw     = the untouched files in WORKING_FOLDER/INPUT_DATASETS/
                  # cleaned = the output of 03_read_cleaning_fastp_fastplong.sh

THREADS=12

DATA="WORKING_FOLDER/INPUT_DATASETS/Metagenome_Assembly_and_Annotation"
OUT="WORKING_FOLDER/RESULTS/metagenomics/02_Read_QC"
CLEANED="WORKING_FOLDER/RESULTS/metagenomics/03_Read_Cleaning"
# ============================================================

if [[ "$WHICH" == "cleaned" ]]; then
    R1="$CLEANED/short_R1.clean.fastq.gz"
    R2="$CLEANED/short_R2.clean.fastq.gz"
    LONG_READS="$CLEANED/long_reads.clean.fastq.gz"
else
    R1="$DATA/raw_data/short_reads_R1.fastq"
    R2="$DATA/raw_data/short_reads_R2.fastq"
    LONG_READS="$DATA/raw_data/long_reads.fastq"
fi

mkdir -p "$OUT/$WHICH/short_reads"

# Two input files, so two threads. More would sit idle.
fastqc "$R1" "$R2" \
    --outdir "$OUT/$WHICH/short_reads" \
    --threads 4

NanoStat \
    --fastq "$LONG_READS" \
    --threads "$THREADS" \
    > "$OUT/$WHICH/long_reads_nanostat.txt"

cat "$OUT/$WHICH/long_reads_nanostat.txt"
