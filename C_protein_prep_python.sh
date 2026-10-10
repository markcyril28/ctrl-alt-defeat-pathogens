#!/usr/bin/env bash
# ============================================================================
# Program C: Protein prep — one clean sequence per target, ready to submit
# ============================================================================
# Program B pulled the metagenome's own version of each reference gene out of
# the assembly, as DNA and as protein, grouped by the reference gene it
# matched. Those files are grouped for reading. The modelling services want
# the opposite: one sequence, on its own, with nothing in it they will choke
# on. This program does that split and the cleaning that goes with it —
# removing stop characters, writing unknown residues as X, collapsing
# identical sequences, and dropping anything too short to model. What each of
# those is for, and what it costs, is at the top of prep_proteins.py.
#
# The sequences come out in a folder per reference organism:
#
#   targets/<organism>/<target>.faa  one sequence each, for program D
#   by_species/<organism>.faa        one organism in one file, for program E
#   all_targets.faa                  all of them in one file, to grep by hand
#
# Four reference organisms went in, so a flat folder would mix them, and
# submitting, downloading and comparing is work you do one organism at a time.
# The organism is the folder the reference gene came from, which is already in
# its name — "3_Klebsiella_pneumoniae__FimH". It says which reference gene the
# sequence matched. It does not say the sequence came from that organism: the
# contigs are a community and this pipeline does no binning.
#
# A thing to keep in front of you while reading the table this prints: these
# sequences came from hypothetical-protein calls that matched a reference gene
# over a few dozen bases. A short, partial, low-identity match is a lead, not
# an identification. The verdict column says whether a sequence can be
# submitted, and nothing at all about whether it is the gene you want.
#
# Needs:  conda activate protein_modeling
# Run from the workshop folder:  bash C_protein_prep_python.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
MIN_AA=50           # shorter than this is dropped: too little to model
WARN_AA=2000        # longer than this is kept, with a note to split it
DROP_DUPLICATES=true   # false keeps every copy of an identical sequence

EXTRACTION="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/B_Gene_Protein_Extraction"
OUT="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/C_Protein_Prep"
PREP="modules/protein_modeling/prep_proteins.py"
# ============================================================

if [[ ! -s "$EXTRACTION/targets.tsv" ]]; then
    echo "ERROR: no extracted targets at $EXTRACTION/"
    echo "Run:  bash B_gene_protein_extraction_bakta.sh"
    exit 1
fi

mkdir -p "$OUT"

# Rebuilt each run, folders and all. A target dropped from program B's table
# should stop being submitted, and the programs after this one read the
# folders, not the table — so an organism folder left over from a previous
# reference set has to go too, not just the files in it.
rm -rf "$OUT/targets" "$OUT/by_species"
rm -f "$OUT/all_targets.faa"

# The helper's notes go to a file rather than straight to the screen, so they
# print after the table they are commenting on instead of before it.
python3 "$PREP" "$EXTRACTION" "$OUT" "$MIN_AA" "$WARN_AA" "$DROP_DUPLICATES" \
    2> "$OUT/prep_notes.txt"

echo
echo "--- prepared sequences ---"
cut -f1,2,6,8,9,10,11,12 "$OUT/prep_report.tsv" | column -t -s $'\t'
cat "$OUT/prep_notes.txt"

echo
echo "report      $OUT/prep_report.tsv"
echo "one each    $OUT/targets/<organism>/*.faa"
echo "per species $OUT/by_species/*.faa"
echo "all in one  $OUT/all_targets.faa"
echo
echo "Next, and they are alternatives worth running both of:"
echo "  bash D_swissmodel_inputs.sh    needs a solved template, fails without one"
echo "  bash E_alphafold3_inputs.sh    needs no template, so it always answers"
