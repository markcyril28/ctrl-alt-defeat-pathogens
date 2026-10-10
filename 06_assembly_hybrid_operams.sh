#!/usr/bin/env bash
# ============================================================================
# Program 06: Hybrid assembly — OPERA-MS
# ============================================================================
# Combines both read types: the short reads give accurate bases, the long
# reads give the order and orientation. It starts from the metaSPAdes contigs,
# so run program 04 first. Worksheet section 8.
#
# OPERA-MS is the exception to the easy-install rule in this workshop. It is
# not a conda package; you build it from source, and it brings its own older
# dependencies:
#
#   git clone https://github.com/CSB5/OPERA-MS.git
#   cd OPERA-MS
#   make
#   perl OPERA-MS.pl check-dependency
#
# Then set OPERA_MS_DIR below. If the dependency check fails on your laptop,
# use the workshop server rather than spending the session on it, and record
# in the worksheet that you skipped this step.
#
# Run from the workshop folder:  bash 06_assembly_hybrid_operams.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
OPERA_MS_DIR=""     # the folder containing OPERA-MS.pl
THREADS=12

CLEANED="WORKING_FOLDER/RESULTS/metagenomics/03_Read_Cleaning"
SHORT_ASSEMBLY="WORKING_FOLDER/RESULTS/metagenomics/04_Assembly_Short"
OUT="WORKING_FOLDER/RESULTS/metagenomics/06_Assembly_Hybrid"
ASSEMBLIES="WORKING_FOLDER/RESULTS/metagenomics/Assemblies"
# ============================================================

# OPERA-MS is not installed by default, so this program is skipped until you
# set the path above. Skipping it does not affect the rest of the pipeline.
if [[ -z "$OPERA_MS_DIR" ]]; then
    echo "OPERA_MS_DIR is empty, so the hybrid assembly is skipped."
    exit 0
fi

mkdir -p "$OUT" "$ASSEMBLIES"

# OPERA-MS wants the long reads uncompressed. This copy is large and can be
# deleted once the assembly has finished.
gzip -cd "$CLEANED/long_reads.clean.fastq.gz" > "$OUT/long_reads.clean.fastq"

# --no-ref-clustering skips the reference-guided step, as in the guide example.
perl "$OPERA_MS_DIR/OPERA-MS.pl" \
    --contig-file "$SHORT_ASSEMBLY/contigs.fasta" \
    --short-read1 "$CLEANED/short_R1.clean.fastq.gz" \
    --short-read2 "$CLEANED/short_R2.clean.fastq.gz" \
    --long-read "$OUT/long_reads.clean.fastq" \
    --no-ref-clustering \
    --out-dir "$OUT" \
    --num-processors "$THREADS"

cp "$OUT/contigs.fasta" "$ASSEMBLIES/operams.fasta"
