#!/usr/bin/env bash
# ============================================================================
# Program D: SWISS-MODEL inputs — one upload per target, and what to record
# ============================================================================
# SWISS-MODEL is a web service, so this program cannot run it. What it can do
# is prepare exactly what the form wants, one file per target, and tell you
# which numbers to bring back. Those numbers are the result; the .pdb is only
# the picture.
#
# SWISS-MODEL builds by homology: it searches for an experimentally solved
# structure whose sequence resembles your target, then threads your sequence
# onto that template's fold. The consequences are worth understanding before
# you submit anything:
#
#   No template, no model. A sequence with nothing solved near it comes back
#   empty. That is not a failure of your target — it is the state of the PDB.
#   Expect this for several of these targets, and expect program E to be the
#   one that answers for them.
#
#   The template is the model's ceiling. Below roughly 30% sequence identity
#   the fold may be right while the details are not, and the details are what
#   you dock into.
#
#   The model is usually apo. Templates' metals, cofactors and ligands are
#   generally not carried over. A metalloenzyme model that arrives without its
#   metal is not the enzyme — program G puts a metal back if you ask it to,
#   and you must be able to say where it came from.
#
# The uploads are grouped the way program C grouped the sequences: one folder
# per reference organism, with a matching folder under models/ for what comes
# back. Fifteen uploads done by hand, then fifteen downloads, then the numbers
# copied off each results page — doing that an organism at a time is the only
# way it stays straight, and submission_manifest.tsv names the exact folder
# each model goes back into.
#
# Where to go:  https://swissmodel.expasy.org/interactive
# Download each model into models/<organism>/, the folder this program creates
# for it, keeping whatever name SWISS-MODEL gave the file. Program F collects
# from there, and models/ is what tells it the model came from SWISS-MODEL, so
# nothing needs renaming.
#
# Needs:  conda activate protein_modeling
# Run from the workshop folder:  bash D_swissmodel_inputs.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
WARN_LENGTH=2000       # longer than this: model it domain by domain instead

PREP="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/C_Protein_Prep"
OUT="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/D_SWISS_MODEL_Inputs"
DROP="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/D_SWISS_MODEL_Inputs/models"
UPLOAD="modules/protein_modeling/swissmodel_upload.py"
# ============================================================

mkdir -p "$OUT" "$DROP"

# Rebuilt each run, so the folders are only ever the current target set: an
# upload left over from a previous selection is the easiest way to submit the
# wrong sequence. Both patterns, because the uploads used to sit loose in
# $OUT/ before they were grouped by organism. Only .fasta is touched, so a
# model already downloaded into models/ is left where it is.
rm -f "$OUT"/*.fasta "$OUT"/*/*.fasta

python3 "$UPLOAD" "$PREP" "$OUT" "$DROP" "$WARN_LENGTH"

echo
echo "--- upload these ---"
cut -f1-6 "$OUT/submission_manifest.tsv" | column -t -s $'\t'

echo
echo "uploads    $OUT/<organism>/*.fasta"
echo "manifest   $OUT/submission_manifest.tsv   (names the folder for each model)"
echo "submit at  https://swissmodel.expasy.org/interactive"
echo "return to  $DROP/<organism>/   (any file name; unzip the .zip first)"
echo
echo "Record per target, in the worksheet, before you close the browser tab:"
echo "  template PDB + chain, sequence identity, coverage, QMEANDisCo global,"
echo "  oligomeric state, and whether any ligand or metal was transferred."
echo
echo "A target that comes back with no model is a result, not a gap: write down"
echo "that the PDB has nothing near it, and let program E answer for that one."
