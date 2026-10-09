#!/usr/bin/env bash
# ============================================================================
# Program 01: Data inventory
# ============================================================================
# Counts the reads and bases in each input file. Fills worksheet sections
# 2.1 to 2.4, and tells you whether every file is complete.
#
# No bioinformatics tool is needed. A FASTQ record is always four lines, so
# line 2 of every 4 holds the bases, and awk can count them.
#
# Run from the workshop folder:  bash 01_data_inventory_awk.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
DATA="Datasets/Metagenome_Assembly_and_Annotation"
OUT="RESULTS/metagenomics/01_Data_Inventory"

SHORT_R1="$DATA/raw_data/short_reads_R1.fastq"
SHORT_R2="$DATA/raw_data/short_reads_R2.fastq"
LONG_READS="$DATA/raw_data/long_reads.fastq"
# ============================================================

mkdir -p "$OUT"

echo -e "file\tbytes\treads\tbases\tmean_length\tmin_length\tmax_length\tintegrity" > "$OUT/inventory.tsv"

for FASTQ in "$SHORT_R1" "$SHORT_R2" "$LONG_READS"; do
    awk -v name="$(basename "$FASTQ")" -v bytes="$(stat -c %s "$FASTQ")" '
        NR % 4 == 2 {
            reads++
            len = length($0)
            bases += len
            if (min == "" || len < min) min = len
            if (len > max) max = len
        }
        END {
            status = "OK"
            if (NR % 4 != 0) status = "TRUNCATED"
            printf "%s\t%s\t%d\t%d\t%.1f\t%d\t%d\t%s\n",
                   name, bytes, reads, bases, bases / reads, min, max, status
        }' "$FASTQ" >> "$OUT/inventory.tsv"
done

column -t -s $'\t' "$OUT/inventory.tsv"

# A TRUNCATED file is not something to work around. Copy it again from the
# original source; do not try to repair reads.
grep -q TRUNCATED "$OUT/inventory.tsv" && exit 1 || exit 0
