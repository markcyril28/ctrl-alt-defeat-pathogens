#!/usr/bin/env bash
# ============================================================================
# Program 09: Annotation — Bakta
# ============================================================================
# Finds the genes in each assembly and gives them names: CDS, tRNA, rRNA and
# more, written out as TSV, GFF3, GenBank and protein FASTA.
# Worksheet sections 12 and 13.
#
# Bakta needs its database, which is a separate download with its own version
# number. Record both: the same assembly annotated against a different
# database version can give different answers.
#
#   bash setup_metagenomics_conda_envs.sh --with-db
#   bakta_db list
#
# One thing worth keeping in mind while you read the output. These contigs
# come from a mixture of organisms, so a gene found here belongs to the
# community, not to any named species. Assigning it to one organism needs
# binning first, which this introductory pipeline does not do. Say so in the
# worksheet.
#
# Needs:  conda activate metagenomics_env
# Run from the workshop folder:  bash 09_annotation.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
DB_DIR="$HOME/workshop/databases/db-light"   # where --with-db put it
MIN_CONTIG_LENGTH=500     # below this, a gene call means very little
THREADS=12

ASSEMBLIES="RESULTS/metagenomics/Assemblies"
OUT="RESULTS/metagenomics/09_Annotation"
# ============================================================

if [[ ! -f "$DB_DIR/version.json" ]]; then
    echo "ERROR: no Bakta database at $DB_DIR"
    echo "Download it with:  bash setup_metagenomics_conda_envs.sh --with-db"
    exit 1
fi

mkdir -p "$OUT"

for FASTA in "$ASSEMBLIES"/*.fasta; do
    NAME=$(basename "$FASTA" .fasta)

    # --skip-plot leaves out the circular genome picture: it assumes one
    # genome, and this is a community.
    bakta \
        --db "$DB_DIR" \
        --output "$OUT/$NAME" \
        --prefix "$NAME" \
        --min-contig-length "$MIN_CONTIG_LENGTH" \
        --threads "$THREADS" \
        --skip-plot \
        --force \
        "$FASTA"
done
