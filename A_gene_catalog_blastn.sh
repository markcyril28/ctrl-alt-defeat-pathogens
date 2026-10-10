#!/usr/bin/env bash
# ============================================================================
# Program A: Gene catalog. Do any hypothetical proteins look like known genes?
# ============================================================================
# Bakta labels most metagenome CDS "hypothetical protein" (gene called,
# function unknown) and lists them in <assembly>.hypotheticals.faa. This
# program searches them with blastn against the workshop's reference genes in
# WORKING_FOLDER/INPUT_DATASETS/from_Database/. Pooling and tables are done by
# gene_catalog.py.
#
# The reference genes come from modules/datasets/fetch_gene_cds.py. Run it
# with --offline first to verify them: a window shifted by one base still
# gives a hit, but no longer to the gene its name claims.
#
# Reading the hits:
#   - Expect short hits. Full-length matches were already named by Bakta, so
#     what remains is a shared domain or a fragment.
#   - One gene can come back as several. Long-read indels shift the frame and
#     split a gene into short ORFs; consecutive locus tags hitting consecutive
#     stretches of one reference are one gene. See reference_coverage.tsv.
#   - No binning was done, so a hit shows the sequence is in the community,
#     not which organism carries it.
#
# Needs:  conda activate protein_modeling
# Run from the workshop folder:  bash A_gene_catalog_blastn.sh
# ============================================================================

set -euo pipefail   # exit on error, unset variable, or failed pipe

# ========================= SETTINGS =========================
ANNOTATION="WORKING_FOLDER/RESULTS/metagenomics/09_Annotation"   # Bakta output, one folder per assembly
REFERENCE="WORKING_FOLDER/INPUT_DATASETS/from_Database"                # reference genes, one folder per organism
OUT="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/A_Gene_Catalog"

EVALUE=1e-5         # max e-value reported
MIN_IDENTITY=70     # min % identity for a match
MIN_LENGTH=40       # min aligned bases for a match
THREADS=4

CATALOG="modules/protein_modeling/gene_catalog.py"
# ============================================================

if [[ ! -d "$ANNOTATION" ]]; then
    echo "ERROR: no annotation at $ANNOTATION"
    echo "Run the metagenomics pipeline first:  bash 09_annotation_bakta.sh"
    exit 1
fi

mkdir -p "$OUT/blast"   # -p: make parents, ok if it exists

# Pool the reference genes into reference_genes.fna
python3 "$CATALOG" references "$REFERENCE" "$OUT"

# Build a BLAST database from the reference genes
makeblastdb \
    -in "$OUT/reference_genes.fna" `# input FASTA` \
    -dbtype nucl `# nucleotide` \
    -out "$OUT/blast/reference_genes" `# database name prefix` \
    > "$OUT/blast/makeblastdb.log"   # log the output

# Pool the hypothetical CDS into hypotheticals.fna
python3 "$CATALOG" queries "$ANNOTATION" "$OUT"

# -task blastn, not the default megablast: megablast needs a 28-base exact
# seed and misses the same gene from another strain or species.
# -max_hsps 1 keeps the best alignment per query/reference pair; a query can
# still hit several references (e.g. the paralogs ExoS and ExoT).
blastn \
    -task blastn `# sensitive mode` \
    -query "$OUT/hypotheticals.fna" `# hypothetical CDS` \
    -db "$OUT/blast/reference_genes" `# reference database` \
    -evalue "$EVALUE" `# e-value cutoff` \
    -max_hsps 1 `# best alignment per pair` \
    -num_threads "$THREADS" `# CPU threads` \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore" `# 6 = tab table, these columns` \
    > "$OUT/blast/blastn.out"   # hit table

# Apply the identity/length cutoffs and write the summary tables
python3 "$CATALOG" tables "$OUT" "$MIN_IDENTITY" "$MIN_LENGTH"

# Rerun blastn on the matched loci only, with pairwise output (-outfmt 0), to
# see base by base whether a hit is one clean stretch (a real domain) or
# scattered mismatches (a coincidence). One run per species, since each run
# writes one file; each still searches the full database, so the e-values
# match the tables. Cheap: few loci, few reference genes.

# Write one matched-loci FASTA per species
python3 "$CATALOG" matched "$OUT"

rm -f "$OUT"/alignments_*.txt          # clear old results; -f: ok if none

for QUERY in "$OUT"/blast/matched_*.fna; do
    [[ -e "$QUERY" ]] || continue      # no matches: glob unexpanded, skip

    SPECIES=$(basename "$QUERY" .fna)  # drop folder and .fna
    SPECIES=${SPECIES#matched_}        # drop "matched_" prefix

    blastn \
        -task blastn `# sensitive mode` \
        -query "$QUERY" `# matched loci` \
        -db "$OUT/blast/reference_genes" `# reference database` \
        -evalue "$EVALUE" `# e-value cutoff` \
        -max_hsps 1 `# best alignment per pair` \
        -num_threads "$THREADS" `# CPU threads` \
        -outfmt 0 `# 0 = pairwise alignment` \
        > "$OUT/alignments_$SPECIES.txt"   # one file per species

    echo "alignments       $OUT/alignments_$SPECIES.txt"
done

for TABLE in search_summary reference_coverage; do
    echo
    echo "--- $TABLE ---"
    column -t -s $'\t' "$OUT/$TABLE.tsv"   # -t: align columns; -s: tab-separated
done

echo
echo "every hit        $OUT/hits.tsv"
echo "the matches      $OUT/matches.tsv"
echo "the alignments   $OUT/alignments_<species>.txt   (each match drawn out, base by base)"
echo "reference genes  $OUT/reference_genes.tsv"
echo "all proteins     $OUT/proteins.faa   (to grep a locus tag by hand)"
