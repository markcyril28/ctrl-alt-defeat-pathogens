#!/usr/bin/env bash
# ============================================================================
# Program 10: AMR and virulence screening — AMRFinderPlus and ABRicate
# ============================================================================
# Searches each assembly for genes associated with antimicrobial resistance
# and with virulence. Two tools, because they ask the question differently:
# AMRFinderPlus applies NCBI's curated rules, ABRicate does a similarity
# search against whichever database you name, VFDB among them.
# Worksheet section 14.
#
# Read the results carefully. A hit means a sequence in your assembly
# resembles a database entry closely enough to pass the thresholds below.
# It is not proof of resistance, of expression, or of a phenotype. And a
# missing hit is not proof of absence: coverage, assembly quality, the
# thresholds and the database's own scope all affect what is found.
#
# No organism is given to AMRFinderPlus on purpose. The --organism option
# applies species-specific rules, and a metagenome is not one species.
#
# Needs:  conda activate meta_env
# Run from the workshop folder:  bash 10_amr_virulence_screening_amrfinderplus_abricate.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
# ABRicate databases. Remove one from this line to skip it.
# See what is installed with:  abricate --list
ABRICATE_DATABASES="vfdb ncbi"

IDENT_MIN=0.9       # AMRFinderPlus: minimum identity, as a fraction
COVERAGE_MIN=0.5    # AMRFinderPlus: minimum coverage of the reference gene
MINID=80            # ABRicate: minimum identity, as a percentage
MINCOV=60           # ABRicate: minimum coverage, as a percentage
THREADS=12

ASSEMBLIES="RESULTS/metagenomics/Assemblies"
OUT="RESULTS/metagenomics/10_AMR_Virulence_Screening"
# ============================================================

mkdir -p "$OUT"

for FASTA in "$ASSEMBLIES"/*.fasta; do
    NAME=$(basename "$FASTA" .fasta)

    # --plus adds the stress, biocide and virulence gene sets.
    amrfinder \
        --nucleotide "$FASTA" \
        --output "$OUT/$NAME.amrfinder.tsv" \
        --ident_min "$IDENT_MIN" \
        --coverage_min "$COVERAGE_MIN" \
        --threads "$THREADS" \
        --plus

    for DB in $ABRICATE_DATABASES; do
        abricate \
            --db "$DB" \
            --minid "$MINID" \
            --mincov "$MINCOV" \
            --threads "$THREADS" \
            "$FASTA" \
            > "$OUT/$NAME.abricate.$DB.tsv"
    done
done
