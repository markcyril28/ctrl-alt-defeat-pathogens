#!/usr/bin/env bash
# ============================================================================
# Program 07: Assembly evaluation — metaQUAST
# ============================================================================
# Measures every assembly in one run, so the report puts them side by side:
# contig count, total length, largest contig, N50, L50 and GC.
# Worksheet section 9.
#
# The dataset ships with the four genomes the mock community was built from.
# Pointing metaQUAST at them adds the genome fraction recovered per organism,
# which is the only way to answer "what did we actually get?" — contiguity
# alone cannot tell you that. Worksheet section 15.
#
# Do not choose an assembly on N50 alone. One long wrong contig raises N50.
# Read total length, genome fraction and misassemblies together.
#
# Needs:  conda activate metagenomics_env
# Run from the workshop folder:  bash 07_assembly_evaluation.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
MIN_CONTIG=500         # ignore contigs shorter than this, in bases
THREADS=12

DATA="Datasets/Metagenome Assembly and Annotation"
ASSEMBLIES="RESULTS/metagenomics/Assemblies"
OUT="RESULTS/metagenomics/07_Assembly_Evaluation"
REFERENCES="$DATA/ref_genomes"      # set to "" to run without references
# ============================================================

mkdir -p "$OUT"

# Collect the assemblies, and their names, so the report columns read
# "metaspades" rather than a long file path.
FASTA_LIST=()
LABEL_LIST=""
for FASTA in "$ASSEMBLIES"/*.fasta; do
    FASTA_LIST+=("$FASTA")
    LABEL_LIST="$LABEL_LIST,$(basename "$FASTA" .fasta)"
done
LABEL_LIST="${LABEL_LIST#,}"    # drop the comma left at the front

# metaQUAST takes a folder of reference genomes here, or a comma-separated
# list of files.
metaquast "${FASTA_LIST[@]}" \
    --labels "$LABEL_LIST" \
    --min-contig "$MIN_CONTIG" \
    --threads "$THREADS" \
    -r "$REFERENCES" \
    -o "$OUT"
