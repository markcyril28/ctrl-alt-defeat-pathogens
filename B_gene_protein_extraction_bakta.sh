#!/usr/bin/env bash
# ============================================================================
# Program B: Gene and protein extraction — the metagenome's own copy of each match
# ============================================================================
# Program A ran blastn and kept, for each hypothetical CDS, its best match to
# a reference gene — if the match cleared A's identity and length thresholds.
# This program goes back to the metagenomics annotation for the genes
# themselves: the DNA from Bakta's .ffn, the protein from its .faa, and where
# the gene sits from its .tsv. Each gene is filed under the reference gene it
# matched, one FASTA per reference.
#
# Two settings, both judgements:
#
#   REFERENCES   which reference genes you care about. Empty keeps them all.
#   MIN_AA       a 40-residue protein is not worth modeling, and most of what
#                comes back here is short: a gene Bakta could match end to end
#                would not have been on the hypothetical list.
#
# Within each file the genes are in the order they land on the reference.
# When consecutive locus tags on one contig tile one reference — one covering
# the start, the next the middle — they are pieces of one gene, split by
# frameshifts in the assembly, not several genes. The contig and coordinates
# in targets.tsv are how you tell.
#
# A match is similarity to a reference sequence, nothing more. It is not proof
# of function, of expression, or of a phenotype, and the gene belongs to the
# community rather than to any named species — this pipeline does no binning.
#
# extract_targets.py does the work, and refuses to run on a stale catalog:
# Bakta assigns locus tags per run, so if the annotation changed since program
# A searched it, A's matches now name different genes.
#
# Needs:  conda activate protein_modeling
# Run from the workshop folder:  bash B_gene_protein_extraction_bakta.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
# Reference genes to keep, by gene name ("TNT_CpnT FimH") or by the full name
# program A gave them ("1_Mycobacterium_tuberculosis__TNT_CpnT"). Empty = all.
REFERENCES=""
MIN_AA=50           # drop proteins shorter than this many residues

CATALOG="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/A_Gene_Catalog"
ANNOTATION="WORKING_FOLDER/RESULTS/metagenomics/09_Annotation"
OUT="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/B_Gene_Protein_Extraction"
EXTRACT="modules/protein_modeling/extract_targets.py"
# ============================================================

if [[ ! -s "$CATALOG/matches.tsv" ]]; then
    echo "ERROR: no BLAST matches at $CATALOG/matches.tsv"
    echo "Run:  bash A_gene_catalog_blastn.sh"
    exit 1
fi

mkdir -p "$OUT/genes" "$OUT/proteins"

# Rebuilt each run, so a reference with no targets this time does not keep
# its file from the last run.
rm -f "$OUT"/genes/*.fna "$OUT"/proteins/*.faa

python3 "$EXTRACT" "$CATALOG" "$ANNOTATION" "$OUT" "$REFERENCES" "$MIN_AA"

echo
echo "--- targets ---"
cut -f1,4,5,7,8,10-14 "$OUT/targets.tsv" | column -t -s $'\t'

echo
echo "target table  $OUT/targets.tsv"
echo "genes (DNA)   $OUT/genes/<reference>.fna"
echo "proteins      $OUT/proteins/<reference>.faa"
echo "next          bash C_protein_prep_python.sh"
