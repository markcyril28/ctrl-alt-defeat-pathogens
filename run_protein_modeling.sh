#!/usr/bin/env bash
# run_protein_modeling.sh — scripted Part I/II of the Protein Modeling Workshop
#
# Runs the reproducible (non-GUI) half of docs/Drafts/Protein_Modeling_Workshop_Manual.md:
# load the structure dataset, export clean receptors, obtain models for the targets
# that have only a sequence, and judge how much each model can be trusted.
#
# Driven entirely by the CONTROL PANEL variables below — this script takes no
# command-line arguments. Set STEP, then run:  bash run_protein_modeling.sh
#
# Steps (set STEP to one of these):
#   session      Load every PDB, colour by chain, save a combined PyMOL session   (Step 4)
#   clean        Export docking-ready receptors, one per target recipe            (Step 6)
#   check        Report which targets have an AlphaFold DB model                  (Step 8a)
#   fetch        Download the available AlphaFold models into Datasets/MODELS     (Step 8b)
#   plddt        Mean pLDDT per model and per region of interest                  (Step 8c)
#   validate     Superpose each model on its crystal structure (super + align)    (Step 8d)
#   identity     Sequence identity of BoNT/C1 to its BoNT/A and /B templates      (Step 8e)
#   trim         Cut low-confidence tails and linkers off the models              (Step 8f)
#   zn-transfer  Give the BoNT/E and /F models the catalytic Zn they lack         (Step 10)
#   all          session, clean, check, fetch, plddt, validate, identity, trim, zn-transfer
#
# Models — fetched, trimmed and Zn-loaded — always use the MODEL_ prefix, never PDB_:
# predicted and experimental coordinates carry different kinds of error.
#
# Needs: pymol (open-source), python 3 with biopython. Docking tools are not used here —
# see run_docking.sh for those.

set -euo pipefail

# ────────────────────────────────────────────────────────────────────
# CONTROL PANEL — set the run here; the script takes no arguments.
# ────────────────────────────────────────────────────────────────────

# ── What to run ─────────────────────────────────────────────────────
# session | clean | check | fetch | plddt | validate | identity | trim | zn-transfer | all
STEP=""

# ── Which dataset ──────────────────────────────────────────────────
# Point to one of the run_protein_modeling_D*.toml files to run only that species.
# Leave empty to run all four species (original behaviour).
CONFIG=""               # e.g. "run_protein_modeling_D1_Mycobacterium_tuberculosis.toml"

# ── Environment and paths ───────────────────────────────────────────
ENV_WANTED="${OLOGIST_ENV:-protein_modeling}"   # conda env name; OLOGIST_ENV wins
DATASETS=""             # dataset folder; empty = <script dir>/Datasets
OUT_DIR=""              # receptors/sessions output folder; empty = <script dir>/RESULTS/2_Protein_Model_Visualization

# ── Behavior toggles ────────────────────────────────────────────────
DRY_RUN=false           # true = print what would run, touch nothing
# ────────────────────────────────────────────────────────────────────

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -z "$DATASETS" ]] && DATASETS="$ROOT/Datasets"
[[ -d "$DATASETS" ]] && DATASETS="$(cd "$DATASETS" && pwd)"
[[ -z "$OUT_DIR"  ]] && OUT_DIR="$ROOT/RESULTS/2_Protein_Model_Visualization"

# ── Dataset config ──────────────────────────────────────────────────
OW_CONFIG_JSON=""
if [[ -n "$CONFIG" ]]; then
    _cfg_path="$CONFIG"
    [[ "$_cfg_path" = /* ]] || _cfg_path="$ROOT/$_cfg_path"
    [[ -f "$_cfg_path" ]] || { printf '\033[0;31m[FAIL]\033[0m config not found: %s\n' "$_cfg_path" >&2; exit 1; }
fi
export OW_CONFIG_JSON   # set to a temp JSON path in preflight, after TMP is ready

# ── Output helpers (same vocabulary as setup_conda_envs.sh) ─────────
GREEN='\033[0;32m'  YELLOW='\033[0;33m'  RED='\033[0;31m'  BLUE='\033[0;34m'  NC='\033[0m'
info()  { printf "${GREEN}[OK]${NC}   %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC} %s\n" "$*"; }
fail()  { printf "${RED}[FAIL]${NC} %s\n" "$*" >&2; }
step()  { printf "\n${BLUE}━━━ %s ━━━${NC}\n" "$*"; }

# ── Validate the control-panel settings ─────────────────────────────
case "$STEP" in
    session|clean|check|fetch|plddt|validate|identity|trim|zn-transfer|all) ;;
    "") fail "STEP is not set. Edit the CONTROL PANEL at the top of this script."
        fail "Choose one of: session clean check fetch plddt validate identity trim zn-transfer all"
        exit 1 ;;
    *)  fail "unknown STEP: '$STEP'"
        fail "Choose one of: session clean check fetch plddt validate identity trim zn-transfer all"
        exit 1 ;;
esac

RECEPTORS="$OUT_DIR/receptors"
TMP=""
cleanup() { [[ -n "$TMP" && -d "$TMP" ]] && rm -rf "$TMP"; }
trap cleanup EXIT

# ── Environment resolution ──────────────────────────────────────────
tools_present() {
    local dir="$1" tool
    for tool in pymol python; do
        if [[ -n "$dir" ]]; then [[ -x "$dir/$tool" ]] || return 1
        else command -v "$tool" &>/dev/null || return 1; fi
    done
}

resolve_env() {
    if [[ -z "$ENV_WANTED" ]] && tools_present ""; then
        info "using pymol/python already on PATH ($(command -v pymol))"
        return
    fi
    local candidates=() name prefix
    [[ -n "$ENV_WANTED" ]] && candidates+=("$ENV_WANTED")
    candidates+=(protein_modeling)
    for name in "${candidates[@]}"; do
        prefix="$(conda info --envs 2>/dev/null | awk -v n="$name" '$1==n {print $NF}')" || true
        [[ -z "$prefix" ]] && continue
        if tools_present "$prefix/bin"; then
            export PATH="$prefix/bin:$PATH"
            info "using conda environment '$name' ($prefix)"
            return
        fi
        warn "environment '$name' exists but lacks pymol or python"
    done
    fail "no environment with pymol + python found. Run: bash setup_conda_envs.sh"
    fail "then set ENV_WANTED in the CONTROL PANEL if your environment is named differently."
    exit 1
}

run_pymol() {
    local script="$1"
    if $DRY_RUN; then printf '  [dry-run] pymol -cq %s\n' "$script"; return 0; fi
    ( cd "$DATASETS" && pymol -cq "$script" )
}

run_python() {
    local script="$1"
    if $DRY_RUN; then printf '  [dry-run] python %s\n' "$script"; return 0; fi
    ( cd "$DATASETS" && python "$script" )
}

preflight() {
    resolve_env
    TMP="$(mktemp -d)"

    if [[ -n "$CONFIG" ]]; then
        OW_CONFIG_JSON="$TMP/config.json"
        python "$ROOT/modules/parse_config.py" to-json "$_cfg_path" "$OW_CONFIG_JSON" \
            || { fail "could not parse $CONFIG"; exit 1; }
        export OW_CONFIG_JSON
        local species_name species_dir
        species_name="$(python "$ROOT/modules/parse_config.py" species-name "$OW_CONFIG_JSON")"
        species_dir="$(python "$ROOT/modules/parse_config.py" species-dir "$OW_CONFIG_JSON")"
        info "config: $CONFIG — $species_name ($species_dir)"
        [[ -d "$DATASETS/$species_dir" ]] || { fail "species directory not found: $DATASETS/$species_dir"; exit 1; }
    else
        local species_dirs=()
        shopt -s nullglob; species_dirs=("$DATASETS"/[1-9]_*/); shopt -u nullglob
        [[ ${#species_dirs[@]} -gt 0 ]] || { fail "no species directories in $DATASETS"; exit 1; }
    fi

    $DRY_RUN || mkdir -p "$RECEPTORS"
    export OW_DATASETS="$DATASETS" OW_RECEPTORS="$RECEPTORS"
}

# ── Step 4 — combined session ───────────────────────────────────────
do_session() {
    step "Step 4 — load the whole dataset into one session"
    run_pymol "$ROOT/modules/protein_modeling/build_session.py"
    info "session written to $RECEPTORS/toxin_dataset.pse"
}

# ── Step 6 — clean receptors ────────────────────────────────────────
do_clean() {
    step "Step 6 — export docking-ready receptors"
    run_pymol "$ROOT/modules/protein_modeling/clean_receptors.py"
    info "receptors written to $RECEPTORS/"
}

# ── Step 8a/8b — AlphaFold models ───────────────────────────────────
_fetch_accessions() {
    [[ -z "$OW_CONFIG_JSON" ]] && return
    python "$ROOT/modules/parse_config.py" accessions "$OW_CONFIG_JSON"
}

do_check() {
    step "Step 8a — which targets need a model"
    if $DRY_RUN; then printf '  [dry-run] python fetch_models.py --check\n'; return; fi
    local acc; acc="$(_fetch_accessions)"
    ( cd "$DATASETS" && python fetch_models.py --check $acc )
}

do_fetch() {
    step "Step 8b — download the available AlphaFold models"
    if $DRY_RUN; then printf '  [dry-run] python fetch_models.py\n'; return; fi
    local acc; acc="$(_fetch_accessions)"
    ( cd "$DATASETS" && python fetch_models.py $acc )
    local search_dir="$DATASETS"
    [[ -n "$OW_CONFIG_JSON" ]] && {
        local sd; sd="$(python "$ROOT/modules/parse_config.py" species-dir "$OW_CONFIG_JSON")"
        search_dir="$DATASETS/$sd"
    }
    local model_count; model_count="$(find "$search_dir" -name 'MODEL_*.pdb' | wc -l)"
    info "models in $search_dir ($model_count files)"
    [[ -z "$OW_CONFIG_JSON" ]] && \
        warn "BoNT/C1, /D and /G have no AlphaFold entry — build them by homology (Step 8e)."
}

# ── Step 8c — confidence ────────────────────────────────────────────
do_plddt() {
    step "Step 8c — mean pLDDT, whole chain and per region"
    run_python "$ROOT/modules/protein_modeling/plddt_analysis.py"
}

# ── Step 8d — model vs crystal ──────────────────────────────────────
do_validate() {
    step "Step 8d — superpose each model on its crystal structure"
    run_pymol "$ROOT/modules/protein_modeling/validate_models.py"
}

# ── Step 8e — template identity ─────────────────────────────────────
do_identity() {
    step "Step 8e — identity of BoNT/C1 to its candidate templates"
    run_python "$ROOT/modules/protein_modeling/sequence_identity.py"
}

# ── Step 8f — trim models ───────────────────────────────────────────
do_trim() {
    step "Step 8f — trim low-confidence regions off the models"
    run_pymol "$ROOT/modules/protein_modeling/trim_models.py"
}

# ── Step 10 note — Zn for the BoNT/E and /F models ──────────────────
do_zn_transfer() {
    step "Step 10 — transfer the catalytic Zn into the BoNT/E and /F models"
    run_pymol "$ROOT/modules/protein_modeling/zn_transfer.py"
    warn "A borrowed metal is a modelling assumption. Say so in any result that uses it."
}

# ── Dispatch ────────────────────────────────────────────────────────
preflight

case "$STEP" in
    session)     do_session ;;
    clean)       do_clean ;;
    check)       do_check ;;
    fetch)       do_fetch ;;
    plddt)       do_plddt ;;
    validate)    do_validate ;;
    identity)    do_identity ;;
    trim)        do_trim ;;
    zn-transfer) do_zn_transfer ;;
    all)
        do_session; do_clean; do_check; do_fetch
        do_plddt; do_validate; do_identity; do_trim; do_zn_transfer ;;
esac

step "Done"
echo "results folder       : $OUT_DIR"
echo "receptors + sessions : $RECEPTORS"
echo "models               : $DATASETS/*/MODEL_*.pdb"
echo "next                 : set COMMAND=control in run_docking.sh and run it"
