#!/usr/bin/env bash
# setup_protein_modeling_conda_envs.sh — Create the Ologist Workshop protein modeling conda environment
# Run from: Ubuntu/WSL terminal  or  Mac Terminal
# Usage:    bash setup_protein_modeling_conda_envs.sh
#           bash setup_protein_modeling_conda_envs.sh --dry-run     # show what would be created
#           bash setup_protein_modeling_conda_envs.sh --no-mamba    # install with conda only
#
# Environment created (name set at the top of the script):
#   protein_modeling — python 3.11, the lettered pipeline's tools (BLAST, PyMOL, Meeko, RDKit, Vina, PDBFixer) and
#                      the manuals' session tools (Biopython, pandas, FastQC, SeqKit, BWA, samtools, bcftools)
#
# Idempotent: an environment that already exists is kept, and any package missing from it is added.
# Channel priority: conda-forge first, bioconda second (required by bioconda).
# bioconda packages (samtools, bwa, etc.) require Linux or macOS.
# Every run is also written to logs/setup_protein_modeling_conda_envs_<date>_<time>.log next to this script.
#
# This is the environment the lettered pipeline runs in. Programs A to J take the genes the metagenomics
# pipeline found all the way to a docking score, and between them they need:
#
#   A gene catalog  BLAST    D SWISS-MODEL   python      G receptor prep  PyMOL + Meeko + PDBFixer
#   B extraction    python   E AlphaFold3    python      H ligand prep    RDKit + Meeko
#   C protein prep  python   F model QC      PyMOL       I docking        Vina
#                                                        J report         python
#
# F and G run PyMOL inside their helpers, as a Python library (from pymol import cmd), not as the viewer.
# G and H call Meeko's mk_prepare_receptor.py and mk_prepare_ligand.py, and I calls the vina program.
# G also rebuilds the side-chain atoms a crystal left out, with PDBFixer. Without it Meeko deletes those
# residues without saying so, even when one of them lines the pocket.
#
# Each lettered program is a short shell script: the settings you change, a loop, and a call to the tool
# or to a helper in modules/. The helpers are where the table and sequence handling lives, and each one
# starts with a comment explaining what it does and why that step matters.
#
# Programs A and B read RESULTS/metagenomics/, so run the metagenomics pipeline first — its tools live in
# a separate environment on purpose (the assemblers would otherwise force conda to downgrade these):
#   bash setup_metagenomics_conda_envs.sh --with-db
#

set -euo pipefail

# ── Conda environment ───────────────────────────────────────────────
ENV_NAME="protein_modeling"        # the env students activate
# ────────────────────────────────────────────────────────────────────

# ── Configuration ──────────────────────────────────────
# mamba installs faster and is used when present; if it cannot build the env, conda retries it
USE_MAMBA=true
PYTHON_VERSION="3.11"
CHANNEL_MAIN="conda-forge"
CHANNEL_BIO="bioconda"
# Programs A to J: BLAST (A), PyMOL (F, G), Meeko (G, H), PDBFixer (G), RDKit (H), Vina (I)
PIPELINE_PACKAGES="blast pymol-open-source meeko pdbfixer rdkit vina"
# The command-line sessions in the workshop manuals; programs A to J run without these
SESSION_PACKAGES="biopython pandas fastqc seqkit bwa samtools bcftools"
# ────────────────────────────────────────────────────────────────────

DRY_RUN=false

for arg in "$@"; do
    case "$arg" in
        --dry-run)   DRY_RUN=true ;;
        --no-mamba)  USE_MAMBA=false ;;
        -h|--help)
            sed -n '2,15p' "$0"
            exit 0 ;;
        *)
            echo "Unknown option: $arg"; exit 1 ;;
    esac
done

# ── Logging ─────────────────────────────────────────────────────────
# The script reruns itself through tee: the terminal shows the usual output, and logs/ keeps a copy
# without colour codes — including conda's full solver error, which scrolls out of view.
if [[ -z "${SETUP_LOG:-}" ]]; then
    LOG_DIR="$(cd "$(dirname "$0")" && pwd)/logs"
    mkdir -p "$LOG_DIR"
    export SETUP_LOG="$LOG_DIR/$(basename "$0" .sh)_$(date +%Y%m%d_%H%M%S).log"
    echo "# $(date)  bash $0 $*" > "$SETUP_LOG"
    # awk flushes every line, so a run stopped with Ctrl-C still leaves a complete log. A progress meter
    # redraws its line with carriage returns; the log keeps only the last state of it.
    bash "$0" "$@" 2>&1 | tee >(awk '{ gsub(/\033\[[0-9;]*m/, "")
                                      if (match($0, /\r[^\r]+$/)) $0 = substr($0, RSTART + 1)
                                      print; fflush() }' >> "$SETUP_LOG")
    exit "${PIPESTATUS[0]}"
fi

# ── Helpers ─────────────────────────────────────────────────────────
GREEN='\033[0;32m'  YELLOW='\033[0;33m'  RED='\033[0;31m'  NC='\033[0m'

info()  { printf "${GREEN}[OK]${NC}  %s\n" "$*"; }
warn()  { printf "${YELLOW}[SKIP]${NC} %s\n" "$*"; }
fail()  { printf "${RED}[FAIL]${NC} %s\n" "$*"; }

# grep reads to the end, not -q: under pipefail, grep -q quitting early can SIGPIPE conda and hide the env
env_exists() { conda info --envs 2>/dev/null | grep -E "^$1 |/$1\$" >/dev/null; }

env_prefix() { conda info --envs 2>/dev/null | awk -v name="$1" '$1==name {print $NF}'; }

# Run a tool from the env at prefix $1 without activating it — avoids "conda run" activation bugs in conda ≥26.
# The env's bin goes first on PATH, so FastQC finds the env's Java, not the stub /usr/bin/java on a Mac.
# PyMOL is a Qt GUI program; the offscreen platform lets it answer --version on headless WSL.
in_env() { local prefix="$1"; shift; PATH="$prefix/bin:$PATH" CONDA_PREFIX="$prefix" QT_QPA_PLATFORM=offscreen "$@"; }

# Version of Python module $1, imported by the python on PATH — a broken library fails here, not in a session
py_module() { python -c 'import importlib, sys; print(importlib.import_module(sys.argv[1]).__version__)' "$1"; }

# bwa has no version option: bare "bwa" prints its usage, which names the version, and exits 1
bwa_version() { { bwa || true; } 2>&1 | awk '/^Version:/ {print "bwa", $2; v = 1} END {exit !v}'; }

# Meeko ships receptor and ligand preparation as command-line scripts, and programs G and H call them by
# name rather than through the library — so an install with the module but not the scripts must fail here.
# --help is the only flag they answer without an input file. It also lists the receptor script's options,
# so a Meeko that lacks one program G passes fails here rather than halfway through G.
meeko_cli() {
    local receptor_help flag
    receptor_help="$(mk_prepare_receptor.py --help)" || return 1
    for flag in --read_pdb --write_pdbqt --write_vina_box --default_altloc --box_center --box_size; do
        if ! grep -q -- "$flag" <<< "$receptor_help"; then
            echo "mk_prepare_receptor.py has no $flag option"
            return 1
        fi
    done
    mk_prepare_ligand.py --help >/dev/null || return 1
    echo "mk_prepare_receptor.py, mk_prepare_ligand.py"
}

# Programs F and G use PyMOL as a Python library: F reads the downloaded .cif models and writes a .pdb of
# each, G cleans those. A residue taken through the same .cif-to-.pdb path checks that without a display.
pymol_api() {
    python - <<'PYTHON'
import os, tempfile
from pymol import cmd
cmd.fragment('ala')
with tempfile.TemporaryDirectory() as folder:
    cmd.save(os.path.join(folder, 'model.cif'), 'ala')
    cmd.delete('all')
    cmd.load(os.path.join(folder, 'model.cif'), 'check')
    cmd.save(os.path.join(folder, 'model.pdb'), 'check')
    assert os.path.getsize(os.path.join(folder, 'model.pdb')) > 0, 'no .pdb written'
assert cmd.count_atoms('check') > 0, 'no atoms read from the .cif'
print(f'pymol {cmd.get_version()[0]}, .cif to .pdb ok')
PYTHON
}

# Program G rebuilds missing side-chain atoms with PDBFixer, which imports OpenMM. PDBFixer has no
# __version__ of its own, so the check reports OpenMM's. An env without either would not fail G: G would
# print a note that PDBFixer is not installed and carry on, leaving Meeko to delete the incomplete residues.
pdbfixer_api() {
    python -c 'import openmm, pdbfixer; from openmm.app import PDBFile; print("pdbfixer, openmm", openmm.__version__)'
}

# Program H turns a SMILES into a 3-D conformer. An RDKit that imports but cannot embed one fails here,
# in setup, instead of in the middle of a session.
rdkit_embed() {
    python - <<'PYTHON'
from rdkit import Chem
from rdkit.Chem import AllChem
molecule = Chem.AddHs(Chem.MolFromSmiles('OCC(O)CO'))
parameters = AllChem.ETKDGv3()
parameters.randomSeed = 42
assert AllChem.EmbedMolecule(molecule, parameters) == 0, 'no conformer embedded'
print('ETKDG conformer ok')
PYTHON
}

# --override-channels keeps a "defaults" entry in ~/.condarc out of the solve (conda and mamba both read it).
CHANNELS=(--override-channels -c "$CHANNEL_MAIN" -c "$CHANNEL_BIO")

# Create env $1 from the package specs that follow. Returns 1 if it cannot be built.
# An env that already exists is kept and topped up instead: the specs go to "install", which adds a package
# the env lacks (one added to the lists since it was built) and leaves the rest as it is. python= is left out
# so the env's Python is not touched. If the top-up fails the env is not removed — the verify step below
# then names the tool that is missing.
# mamba goes first when $INSTALLER is mamba. Both solve with libmamba, but mamba 2.8's sharded index once
# reported conda-forge packages as missing, so a mamba failure is retried with conda rather than reported.
make_env() {
    local name="$1"
    shift
    echo ""
    echo "━━━ $name ━━━"
    if env_exists "$name"; then
        local spec
        local -a add=()
        for spec in "$@"; do
            [[ $spec == python=* ]] || add+=("$spec")
        done
        warn "$name already exists — adding any package it lacks"
        if $DRY_RUN; then
            echo "  [dry-run] $INSTALLER install -y -n $name --strict-channel-priority ${CHANNELS[*]} ${add[*]}"
            return 0
        fi
        if [[ $INSTALLER == mamba ]]; then
            if MAMBA_ROOT_PREFIX="$CONDA_BASE" mamba install -y -n "$name" --strict-channel-priority \
                    "${CHANNELS[@]}" "${add[@]}"; then
                info "$name is up to date"
                return 0
            fi
            warn "mamba could not update $name — retrying with conda"
        fi
        if conda install -y -n "$name" --strict-channel-priority "${CHANNELS[@]}" "${add[@]}"; then
            info "$name is up to date"
        else
            fail "Could not add the missing packages to $name — conda's reason is above, and in the log"
        fi
        return 0
    fi
    if $DRY_RUN; then
        echo "  [dry-run] $INSTALLER create -y -n $name --strict-channel-priority ${CHANNELS[*]} $*"
        if [[ $INSTALLER == mamba ]]; then
            echo "  [dry-run] if mamba fails: conda create with the same arguments"
        fi
        return 0
    fi
    echo "Creating $name with $INSTALLER (several minutes of downloading is normal) ..."
    if [[ $INSTALLER == mamba ]]; then
        # MAMBA_ROOT_PREFIX puts the env in conda's envs folder, where conda (and env_exists) look for it
        if MAMBA_ROOT_PREFIX="$CONDA_BASE" mamba create -y -n "$name" --strict-channel-priority \
                "${CHANNELS[@]}" "$@" && env_exists "$name"; then
            info "$name created"
            return 0
        fi
        warn "mamba could not create $name — retrying with conda"
        conda env remove -y -n "$name" &>/dev/null || true   # clear any half-built env
    fi
    if conda create -y -n "$name" --strict-channel-priority "${CHANNELS[@]}" "$@"; then
        info "$name created"
        return 0
    fi
    fail "Could not create $name — conda's reason is above, and in the log"
    conda env remove -y -n "$name" &>/dev/null || true   # clear any half-built env
    return 1
}

FAILED=0
ENV_PREFIX=""

# ── Preflight ───────────────────────────────────────────────────────
if ! command -v conda &>/dev/null; then
    fail "conda not found. Install Miniforge first (see docs/Drafts/Mine/Bioinformatics_Workshop_Setup_Instructions.md)."
    exit 1
fi
info "conda $(conda --version 2>&1 | awk '{print $2}') found"

PLATFORM="$(uname -s)-$(uname -m)"
info "Platform: $PLATFORM"
info "Log: $SETUP_LOG"

case "$(uname -s)" in
    Linux|Darwin) ;;
    *)
        fail "bioconda packages require Linux or macOS (detected: $(uname -s))."
        fail "Run this script inside WSL Ubuntu, not Windows Miniforge Prompt."
        exit 1 ;;
esac

INSTALLER=conda
if $USE_MAMBA && command -v mamba &>/dev/null; then
    INSTALLER=mamba
    CONDA_BASE="$(conda info --base)"
    info "mamba $(mamba --version 2>&1 | awk 'NR==1 {print $NF}') found — installs use mamba, with conda as the fallback"
elif $USE_MAMBA; then
    warn "mamba not found — installing with conda (Miniforge includes mamba)"
fi

# ── Create environment ──────────────────────────────────────────────
if ! make_env "$ENV_NAME" python="$PYTHON_VERSION" $PIPELINE_PACKAGES $SESSION_PACKAGES; then
    fail "Every session needs $ENV_NAME — save the log for the trainer: $SETUP_LOG"
    exit 1
fi
ENV_PREFIX="$(env_prefix "$ENV_NAME")"

# ── Verify ──────────────────────────────────────────────────────────
# Everything is checked from the env's bin, the way students run it after conda activate.
# label|command[|fallback command]
TOOLS=(
    # Programs A to J
    "BLAST|blastn -version"                 # A
    "makeblastdb|makeblastdb -version"      # A
    "PyMOL library|pymol_api"               # F, G
    "Meeko|py_module meeko"                 # G, H
    "Meeko CLI|meeko_cli"                   # G, H
    "PDBFixer|pdbfixer_api"                 # G
    "RDKit|py_module rdkit"                 # H
    "RDKit 3-D|rdkit_embed"                 # H
    "Vina|vina --version"                   # I
    # The manuals' sessions
    "PyMOL|pymol --version"
    "Biopython|py_module Bio"
    "pandas|py_module pandas"
    "FastQC|fastqc --version"
    "SeqKit|seqkit version"
    "BWA|bwa_version"
    "samtools|samtools --version"
    "bcftools|bcftools --version"
)

if ! $DRY_RUN; then
    echo ""
    echo "Verifying $ENV_NAME ..."
    if [[ -z "$ENV_PREFIX" ]]; then
        fail "Could not locate $ENV_NAME prefix"
        FAILED=$((FAILED + 1))
    else
        for entry in "${TOOLS[@]}"; do
            IFS='|' read -ra parts <<< "$entry"
            ok=false
            for cmd in "${parts[@]:1}"; do
                read -ra argv <<< "$cmd"
                if out="$(in_env "$ENV_PREFIX" "${argv[@]}" 2>&1)"; then
                    printf "  %-14s %s\n" "${parts[0]}" "${out%%$'\n'*}"
                    ok=true
                    break
                fi
            done
            if ! $ok; then
                fail "${parts[0]} did not answer (${parts[1]})"
                [[ -n "$out" ]] && echo "       ${out##*$'\n'}"   # last line: the error, or the check's reason
                FAILED=$((FAILED + 1))
            fi
        done
        if [[ $FAILED -eq 0 ]]; then
            info "$ENV_NAME verified"
        else
            fail "$FAILED tool(s) failed verification — activate and debug manually"
        fi
    fi
fi

# ── Summary ─────────────────────────────────────────────────────────
echo ""
echo "━━━ Summary ━━━"
conda info --envs 2>/dev/null | grep -E "^$ENV_NAME " || true
echo ""
if [[ $FAILED -gt 0 ]]; then
    fail "$FAILED check(s) failed — see [FAIL] above. Rebuild the env with: conda env remove -n $ENV_NAME, then rerun"
fi
echo "Activate with:   conda activate $ENV_NAME"
echo "Record versions: conda list -n $ENV_NAME > ${ENV_NAME}_versions.txt"
echo "Workshop data:   cd $(cd "$(dirname "$0")" && pwd)/WORKING_FOLDER/INPUT_DATASETS"
echo "Log:             $SETUP_LOG"
echo ""
echo "The lettered pipeline, once the metagenomics pipeline has been run:"
echo "  bash run_protein_modeling.sh      # or one program at a time, A to J"
echo ""
echo "It runs in two sittings. Programs A to E end with files to upload to SWISS-MODEL and AlphaFold3;"
echo "download the models into the models/ folder beside each set of uploads — into the folder for that"
echo "organism, keeping whatever name the service gave them —"
echo "  RESULTS/protein_modeling/for_metagenomics_dataset/D_SWISS_MODEL_Inputs/models/<organism>/"
echo "  RESULTS/protein_modeling/for_metagenomics_dataset/E_AlphaFold3_Inputs/alphafoldserver/models/<organism>/"
echo "then run it again: F to J skip themselves until those files are there."

exit $(( FAILED > 0 ))
