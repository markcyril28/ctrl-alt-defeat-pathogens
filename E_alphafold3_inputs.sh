#!/usr/bin/env bash
# ============================================================================
# Program E: AlphaFold3 inputs — the same targets, asked a different way
# ============================================================================
# SWISS-MODEL needs a template. AlphaFold3 does not: it predicts from the
# sequence and the alignment of its relatives, so it returns a model for
# targets that have no solved homologue at all. Running both on the same
# sequence is the point of doing both — where they agree you have a fold you
# can work with, and where they disagree you have a question.
#
# That cuts both ways, and it is the trap of this program. AlphaFold always
# answers. A sequence with no relatives and no template still comes back as a
# structure, and what comes back is confident-looking ribbon with nothing
# behind it. The pLDDT column is where that shows up, which is why program F
# measures it before anything gets docked.
#
# This writes job files in both dialects, because they are not the same
# format and sending one to the other fails:
#
#   alphafoldserver   the web form at alphafoldserver.com. A file is an array
#                     of jobs, so one upload can carry every target.
#   alphafold3        a local AlphaFold3 install, which also needs the genetic
#                     databases and the model weights — the weights are
#                     requested from Google DeepMind, not downloaded. The JSON
#                     is written here either way, so it is ready when they are.
#
# Reading the output when it comes back:
#
#   pLDDT lives in the B-factor column. It is per-residue confidence on a
#   0-100 scale, not a crystallographic B-factor, and it is not accuracy. Judge
#   it over the site you intend to dock into, not over the whole chain: a 95
#   mean with a 40 loop across the active site is useless for docking.
#
#   PAE says whether two parts are placed correctly relative to each other.
#   Two confident domains can be confidently positioned wrong.
#
#   A seed is a starting point, not a replicate. Several seeds show you how
#   stable the prediction is; one seed shows you nothing about that.
#
#   Asking for a zinc or a cofactor is a hypothesis, not an observation.
#   AlphaFold3 will place what you name whether or not the protein binds it.
#
# The job files are grouped the way program C grouped the sequences: one
# folder per reference organism, and one batch file per organism rather than
# one for everything. That is the unit the rest of the work happens in — you
# upload an organism, wait, download it, and read it against the other models
# of the same organism. A batch that fails or comes back short then belongs to
# an organism instead of to a pile of fifteen.
#
# Download the models into models/<organism>/, the folder this program creates
# for each one, keeping whatever name the server gave each file. Program F
# collects from there, and models/ is what tells it the model came from
# AlphaFold, so nothing needs renaming.
#
# Upload the batch files, not the one-per-target files — one upload per
# organism carries all of that organism's targets. Write down the date you
# submitted and the model version the site reports: neither is in the file
# that comes back, and program J leaves a blank for both.
#
# A local run also needs the databases and the weights:
#   alphafold3 --json_path=<job>.json --db_dir=<databases> --model_dir=<weights>
# The weights are requested from Google DeepMind and are not part of this
# workshop.
#
# Needs:  conda activate protein_modeling
# Run from the workshop folder:  bash E_alphafold3_inputs.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
SERVER_SEEDS=0      # 0 = let the server choose; any number = ask for that many
LOCAL_SEEDS=1       # a local run states its own seeds; more = slower, steadier
BATCH_SIZE=20       # jobs per server upload file
PREFIX="AF3"        # stem for the batch upload file

# Ions and cofactors to include in every job, space separated. Empty by
# default on purpose: add one only for a target where you can say why it
# should be there. They go into every job in the batch, so a cofactor that
# suits one target is being asserted about all of them.
#   IONS="ZN"        e.g. a zinc metalloprotease
#   LIGANDS="NAD"    e.g. an NAD-dependent enzyme (CCD codes)
IONS=""
LIGANDS=""

PREP="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/C_Protein_Prep"
OUT="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/E_AlphaFold3_Inputs"
DROP="WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/E_AlphaFold3_Inputs/alphafoldserver/models"
HELPER="modules/protein_modeling/make_af3_json.py"
# ============================================================

if ! ls "$PREP"/by_species/*.faa > /dev/null 2>&1; then
    echo "ERROR: no prepared sequences at $PREP/by_species/"
    echo "Run:  bash C_protein_prep_python.sh"
    exit 1
fi

mkdir -p "$OUT/alphafoldserver" "$OUT/alphafold3_local" "$DROP"

# Rebuilt each run. A stale job file would be uploaded alongside the current
# ones and come back as a model of a target you are no longer working on.
# Also the flat layout these files used to be written in, before they were
# grouped by organism.
rm -f "$OUT"/alphafoldserver/*.json "$OUT"/alphafold3_local/*.json

# One folder per organism, on both sides: the job files to upload, and the
# empty folder the models come back into. The organisms are whichever ones
# program C ended up with, so this list is never out of step with it.
for SPECIES_FASTA in "$PREP"/by_species/*.faa; do
    SPECIES=$(basename "$SPECIES_FASTA" .faa)
    mkdir -p "$OUT/alphafoldserver/$SPECIES" "$OUT/alphafold3_local/$SPECIES" \
             "$DROP/$SPECIES"
    rm -f "$OUT/alphafoldserver/$SPECIES"/*.json \
          "$OUT/alphafold3_local/$SPECIES"/*.json
done

# The helper takes the whole by_species folder and writes one organism per
# subfolder, so both manifests come out of a single call and stay one table.
python3 "$HELPER" "$PREP/by_species" "$OUT/alphafoldserver" alphafoldserver \
    --seeds "$SERVER_SEEDS" --ions "$IONS" --ligands "$LIGANDS" \
    --batch-size "$BATCH_SIZE" --batch-name "${PREFIX}_server_batch" \
    > "$OUT/server_manifest.tsv"

python3 "$HELPER" "$PREP/by_species" "$OUT/alphafold3_local" alphafold3 \
    --seeds "$LOCAL_SEEDS" --ions "$IONS" --ligands "$LIGANDS" \
    > "$OUT/local_manifest.tsv"

echo
echo "--- server jobs ---"
column -t -s $'\t' "$OUT/server_manifest.tsv"

echo
echo "upload     $OUT/alphafoldserver/<organism>/${PREFIX}_server_batch_*.json"
echo "           one file per organism, to https://alphafoldserver.com"
echo "local jobs $OUT/alphafold3_local/<organism>/"
echo "return to  $DROP/<organism>/   (any file name; unzip each job's .zip first)"
