#!/usr/bin/env bash
# ============================================================================
# Program 05: Long-read assembly — MetaFlye
# ============================================================================
# Assembles the cleaned long reads. A long read can span a repeat that a short
# read cannot, so this assembly is usually far more contiguous than the
# short-read one — but each base is less accurate. Neither is simply better.
# Worksheet section 7.
#
# This is the slowest program in the pipeline. Expect it to run for a while.
#
# READ_TYPE must match how the reads were produced. --nano-raw is for ordinary
# basecalled nanopore reads; use --nano-hq only for Q20+ or SUP basecalling,
# and --pacbio-raw for PacBio. Telling Flye the wrong thing wastes the run.
#
# Needs:  conda activate meta_env
# Run from the workshop folder:  bash 05_assembly_long_metaflye.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
READ_TYPE="--nano-raw"   # --nano-raw | --nano-hq | --pacbio-raw
ITERATIONS=1             # polishing rounds; more is slower, not always better
THREADS=12

CLEANED="WORKING_FOLDER/RESULTS/metagenomics/03_Read_Cleaning"
OUT="WORKING_FOLDER/RESULTS/metagenomics/05_Assembly_Long"
ASSEMBLIES="WORKING_FOLDER/RESULTS/metagenomics/Assemblies"
# ============================================================

mkdir -p "$OUT" "$ASSEMBLIES"

# --meta tells Flye to expect several genomes at different abundances.
flye \
    "$READ_TYPE" "$CLEANED/long_reads.clean.fastq.gz" \
    --meta \
    --out-dir "$OUT" \
    --threads "$THREADS" \
    --iterations "$ITERATIONS"

cp "$OUT/assembly.fasta" "$ASSEMBLIES/metaflye.fasta"
cp "$OUT/assembly_graph.gfa" "$ASSEMBLIES/metaflye.gfa"
