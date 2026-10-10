#!/usr/bin/env bash
# ============================================================================
# Program 04: Short-read assembly — metaSPAdes
# ============================================================================
# Joins the cleaned short reads into contigs. metaSPAdes is SPAdes set up for
# metagenomes: it expects uneven coverage, because in a community the abundant
# organisms are sequenced deeply and the rare ones barely at all.
#
# Expect many contigs. That is not a failure — it is what a mixture of genomes
# looks like when some of them are only partly covered. Worksheet section 6.
#
# If this environment has spades.py but no metaspades.py, use
# "spades.py --meta" below instead. It is the same assembler.
#
# The main output is copied to WORKING_FOLDER/RESULTS/metagenomics/Assemblies/, which is
# where programs 07 to 10 look for assemblies.
#
# Needs:  conda activate meta_env
# Run from the workshop folder:  bash 04_assembly_short_metaspades.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
THREADS=12
MEMORY_GB=20        # SPAdes stops itself rather than exhausting the machine

CLEANED="WORKING_FOLDER/RESULTS/metagenomics/03_Read_Cleaning"
OUT="WORKING_FOLDER/RESULTS/metagenomics/04_Assembly_Short"
ASSEMBLIES="WORKING_FOLDER/RESULTS/metagenomics/Assemblies"
# ============================================================

mkdir -p "$OUT" "$ASSEMBLIES"

metaspades.py \
    -1 "$CLEANED/short_R1.clean.fastq.gz" \
    -2 "$CLEANED/short_R2.clean.fastq.gz" \
    -o "$OUT" \
    --threads "$THREADS" \
    --memory "$MEMORY_GB"

cp "$OUT/contigs.fasta" "$ASSEMBLIES/metaspades.fasta"
cp "$OUT/assembly_graph_with_scaffolds.gfa" "$ASSEMBLIES/metaspades.gfa"
