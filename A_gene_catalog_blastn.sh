#!/usr/bin/env bash
# ============================================================================
# Program A: Gene catalog — do any hypothetical proteins look like known genes?
# ============================================================================
# The metagenomics pipeline ended with Bakta, and Bakta called most CDS
# "hypothetical protein". That is the normal state of a metagenome, not a
# failure: the gene call stands, but no database match explained what the
# protein does. Bakta lists those genes separately, in
# <assembly>.hypotheticals.faa.
#
# This program asks whether any of them resemble the genes this workshop is
# about — the reference genes in
# WORKING_FOLDER/INPUT_DATASETS/from_Database/ — using blastn. It pools the
# reference genes, pools the hypothetical CDS, searches one against the other,
# and writes the tables. The pooling and the tables are in gene_catalog.py,
# which explains each step at the top of the file.
#
# The reference genes are fetched by modules/datasets/fetch_gene_cds.py, which
# downloads each CDS from NCBI by genome coordinate and checks that it reads as
# a gene. Run it with --offline to confirm the files in
# WORKING_FOLDER/INPUT_DATASETS/from_Database/ are the ones it would write.
# This is worth doing before trusting a result: a reference window shifted by
# one base still blasts and still looks like a hit, it just no longer holds the
# gene its name claims.
#
# Things worth knowing before you read the hits:
#
#   Expect short hits. A gene Bakta could match end to end would not be on
#   the hypothetical list. What is left is a stretch of a gene: a shared
#   domain, or a fragment.
#
#   One gene can come back as several. A long-read assembly keeps some of the
#   sequencing errors, and an insertion or deletion shifts the reading frame,
#   so Bakta calls one gene as two or three short ORFs. Consecutive locus tags
#   hitting consecutive stretches of one reference gene are one gene broken by
#   frameshifts. reference_coverage.tsv shows how much of each reference the
#   pieces cover together.
#
#   A locus tag belongs to the community, not to a species. These contigs come
#   from a mixture of organisms and this pipeline does no binning, so a hit to
#   a Klebsiella gene says the sequence is there, not which organism carries it.
#
# Needs:  conda activate protein_modeling
# Run from the workshop folder:  bash A_gene_catalog_blastn.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
ANNOTATION="WORKING_FOLDER/RESULTS/metagenomics/09_Annotation"   # one Bakta folder per assembly
REFERENCE="WORKING_FOLDER/INPUT_DATASETS/from_Database"                # one folder per organism, genes in GENE_*.gene.fna
OUT="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/A_Gene_Catalog"

EVALUE=1e-5         # blastn reports nothing weaker than this
MIN_IDENTITY=70     # percent identity for a hit to count as a match
MIN_LENGTH=40       # aligned bases for a hit to count as a match
THREADS=4

CATALOG="modules/protein_modeling/gene_catalog.py"
# ============================================================

if [[ ! -d "$ANNOTATION" ]]; then
    echo "ERROR: no annotation at $ANNOTATION"
    echo "Run the metagenomics pipeline first:  bash 09_annotation_bakta.sh"
    exit 1
fi

mkdir -p "$OUT/blast"

python3 "$CATALOG" references "$REFERENCE" "$OUT"

# Making a Database out of the reference genes 
makeblastdb \
    -in "$OUT/reference_genes.fna" \
    -dbtype nucl \
    -out "$OUT/blast/reference_genes" > "$OUT/blast/makeblastdb.log"

python3 "$CATALOG" queries "$ANNOTATION" "$OUT"

# -task blastn rather than the default, megablast. Megablast starts from a
# 28-base exact match and is built for near-identical sequence; the same gene
# from another strain or species is often too different for it to find.
#
# -max_hsps 1 keeps the best alignment of each hypothetical to each reference
# gene. A hypothetical can still hit more than one reference: ExoS and ExoT
# are paralogs, so a stretch of one usually finds the other as well.
blastn \
    -task blastn \
    -query "$OUT/hypotheticals.fna" \
    -db "$OUT/blast/reference_genes" \
    -evalue "$EVALUE" \
    -max_hsps 1 \
    -num_threads "$THREADS" \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore" \
    > "$OUT/blast/blastn.out"

python3 "$CATALOG" tables "$OUT" "$MIN_IDENTITY" "$MIN_LENGTH"

# The same search again, over the matched loci only, asking for blastn's
# pairwise output: the one that prints the two sequences with a line of bars
# between them marking the bases that agree. The tables say a hit is 85%
# identical over 115 bases; this is where you see whether those bases are one
# clean stretch or a mismatch every few bases, which is the difference between
# a real shared domain and a coincidence.
#
# One file per bacterial species, because one blastn run writes one file, so
# the split has to happen on the query side. gene_catalog.py wrote one query
# file per species; each is still searched against the whole reference
# database, so the e-values here are the ones the tables report.
#
# Extra searches rather than a second -outfmt on the first one, because a
# blastn run writes a single format. They cost almost nothing: the matched
# loci are a handful, and the database is the few reference genes.
python3 "$CATALOG" matched "$OUT"

rm -f "$OUT"/alignments_*.txt          # a species may not match this time

for QUERY in "$OUT"/blast/matched_*.fna; do
    [[ -e "$QUERY" ]] || continue      # no matches at all, so no query files

    SPECIES=$(basename "$QUERY" .fna)  # matched_1_Mycobacterium_tuberculosis
    SPECIES=${SPECIES#matched_}        # 1_Mycobacterium_tuberculosis

    blastn \
        -task blastn \
        -query "$QUERY" \
        -db "$OUT/blast/reference_genes" \
        -evalue "$EVALUE" \
        -max_hsps 1 \
        -num_threads "$THREADS" \
        -outfmt 0 \
        > "$OUT/alignments_$SPECIES.txt"

    echo "alignments       $OUT/alignments_$SPECIES.txt"
done

for TABLE in search_summary reference_coverage; do
    echo
    echo "--- $TABLE ---"
    column -t -s $'\t' "$OUT/$TABLE.tsv"
done

echo
echo "every hit        $OUT/hits.tsv"
echo "the matches      $OUT/matches.tsv"
echo "the alignments   $OUT/alignments_<species>.txt   (each match drawn out, base by base)"
echo "reference genes  $OUT/reference_genes.tsv"
echo "all proteins     $OUT/proteins.faa   (to grep a locus tag by hand)"
