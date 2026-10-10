#!/usr/bin/env bash
# ============================================================================
# Program F: Model QC — collect what came back, then judge it
# ============================================================================
# Programs D and E sent sequences away. This is where what came back gets
# collected and checked. Every model is a prediction, including the confident
# ones, and the only question that matters is whether it is good enough for
# what you intend to do with it. For docking, that means the pocket, not the
# protein.
#
# Where to put your downloads — next to the files you uploaded, in the folder
# for that organism:
#
#   SWISS-MODEL   WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/D_SWISS_MODEL_Inputs/models/<organism>/
#   AlphaFold3    WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/E_AlphaFold3_Inputs/alphafoldserver/models/<organism>/
#
# Keep whatever name the service gave the file. The drop folder is what tells
# this program which service produced a model, and everything collected is
# copied into WORKING_FOLDER/RESULTS/protein_modeling/Models/ under one naming scheme — the
# single folder programs G and I read. That one stays flat: the organism is
# already in every report, and docking cares which model it has rather than
# which folder it arrived in. It is also derived — this program rebuilds what
# it put there, and a model you place in it by hand is left alone.
#
# Reading the table:
#
#   A high mean hides a bad loop. A chain at 90 with the active-site loop at
#   35 is a good model of the wrong thing. Open it and colour by B-factor.
#
#   The scale is not the same between services. SWISS-MODEL writes QMEANDisCo
#   from 0 to 1; AlphaFold writes pLDDT from 0 to 100; a crystal structure
#   writes a real B-factor, where high is bad instead of good. The helper
#   detects which and says so in the conf_scale column.
#
#   Confidence is not accuracy. It is the model's opinion of itself.
#
#   Two models of the same target are the useful case, not a duplicate. Where
#   the homology model and the prediction agree on the fold you have something
#   to work with; where they disagree, the disagreement is the finding, and
#   neither one settles it.
#
# This program measures, converts and copies. It rejects nothing: a
# low-confidence model is still worth docking into if you report the
# confidence alongside the score.
#
# Needs:  conda activate protein_modeling
# Run from the workshop folder:  bash F_model_qc_pymol.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
CONF_MIN=70         # confidence below this counts as low, on the 0-100 scale
CONVERT="--convert" # write a .pdb beside every .cif, for programs G to I.
                    # Set to "" to leave the .cif files alone.

SWISSMODEL_DROP="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/D_SWISS_MODEL_Inputs/models"
AF3_DROP="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/E_AlphaFold3_Inputs/alphafoldserver/models"
MODELS="WORKING_FOLDER/RESULTS/protein_modeling/Models"      # derived pool, read by G and I
TARGETS="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/C_Protein_Prep/targets"   # one folder per organism
OUT="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/F_Model_QC"
COLLECT="modules/protein_modeling/collect_models.py"
MEASURE="modules/protein_modeling/model_qc.py"
# ============================================================

mkdir -p "$SWISSMODEL_DROP" "$AF3_DROP" "$MODELS" "$OUT"

# Nothing to judge yet is a normal state, not an error: the models arrive by
# download, on whatever day the services get to them.
#
# Two levels deep, which is as far as a download counts: one folder per
# organism inside each drop folder, or left loose beside those folders.
# Deeper is an unzipped job folder, holding five models of one target and the
# templates it was built from, and collect_models.py leaves those alone too.
DOWNLOADED=$(find "$SWISSMODEL_DROP" "$AF3_DROP" "$MODELS" -maxdepth 2 \
    -name "*.pdb" -o -name "*.cif" | wc -l || true)

if [[ "$DOWNLOADED" -eq 0 ]]; then
    echo "No models downloaded yet, so program F is skipped."
    echo "Drop them next to the files you uploaded, in the folder for that"
    echo "organism, keeping whatever name the service gave the file — the drop"
    echo "folder is what tells this program which service produced it, and the"
    echo "two use different confidence scales:"
    echo "  SWISS-MODEL  ->  $SWISSMODEL_DROP/<organism>/"
    echo "  AlphaFold3   ->  $AF3_DROP/<organism>/"
    exit 0
fi

echo "--- collecting downloads ---"
python3 "$COLLECT" "$MODELS" "$TARGETS" "$SWISSMODEL_DROP" "$AF3_DROP"

# The notes go to a file rather than straight to the screen, so they print
# after the table they are commenting on instead of before it.
python3 "$MEASURE" "$MODELS" --conf-min "$CONF_MIN" $CONVERT \
    > "$OUT/model_qc.tsv" 2> "$OUT/qc_notes.txt"

echo
echo "--- model QC ---"
column -t -s $'\t' "$OUT/model_qc.tsv"
cat "$OUT/qc_notes.txt"

echo
echo "ready to dock  $(ls "$MODELS" | grep -c '\.pdb$' || true) .pdb in $MODELS/"
echo "next           bash G_receptor_prep_meeko.sh"
