#!/usr/bin/env bash
# setup_metagenomics_conda_envs.sh — Create the Ologist Workshop metagenomics conda environment
# Run from: Ubuntu/WSL terminal  or  Mac Terminal
# Usage:    bash setup_metagenomics_conda_envs.sh
#           bash setup_metagenomics_conda_envs.sh --dry-run        # show what would be created
#           bash setup_metagenomics_conda_envs.sh --with-db        # also download the annotation databases
#           bash setup_metagenomics_conda_envs.sh --with-db --db-dir=/path/to/databases
#
# Environment created:
#   metagenomics_env — python 3.10, FastQC, NanoStat, fastp, fastplong, SPAdes (metaSPAdes),
#                      Flye (MetaFlye), QUAST (metaQUAST), Bandage, Bakta, AMRFinderPlus, ABRicate
#
# Not installed: OPERA-MS — not a conda package; build it from source only if the trainer asks
#   (Module 7 Step 5 of docs/Drafts/Mine/Bioinformatics_Workshop_Setup_Instructions.md).
#
# Idempotent: skips the environment if it already exists, and the Bakta database if present.
# Kept separate from protein_modeling: the assemblers' dependencies would downgrade its packages.
# bioconda packages require Linux or macOS.

set -euo pipefail

# ── Configuration ──────────────────────────────────────
ENV_NAME="metagenomics_env"
PYTHON_VERSION="3.10"
CHANNEL_MAIN="conda-forge"
CHANNEL_BIO="bioconda"
# QC → cleaning → assembly → evaluation → graph. metaQUAST ships inside the quast package.
CORE_PACKAGES="fastqc nanostat fastp fastplong spades flye quast bandage"
# Annotation and AMR/virulence screening. Dropped (with a warning) if they cannot be solved.
ANNOT_PACKAGES="bakta ncbi-amrfinderplus abricate"
DB_DIR="$HOME/workshop/databases"
# ────────────────────────────────────────────────────────────────────

DRY_RUN=false
WITH_DB=false

for arg in "$@"; do
    case "$arg" in
        --dry-run)   DRY_RUN=true ;;
        --with-db)   WITH_DB=true ;;
        --db-dir=*)  DB_DIR="${arg#--db-dir=}" ;;
        -h|--help)
            sed -n '2,18p' "$0"
            exit 0 ;;
        *)
            echo "Unknown option: $arg"; exit 1 ;;
    esac
done
DB_DIR="${DB_DIR/#\~/$HOME}"   # --db-dir=~/x is not tilde-expanded by the shell

# ── Helpers ─────────────────────────────────────────────────────────
GREEN='\033[0;32m'  YELLOW='\033[0;33m'  RED='\033[0;31m'  NC='\033[0m'

info()  { printf "${GREEN}[OK]${NC}  %s\n" "$*"; }
warn()  { printf "${YELLOW}[SKIP]${NC} %s\n" "$*"; }
fail()  { printf "${RED}[FAIL]${NC} %s\n" "$*"; }

env_exists() {
    conda info --envs 2>/dev/null | grep -qE "^$1 " || \
    conda info --envs 2>/dev/null | grep -qE "/$1\$"
}

# Run a tool from the env without activating it — avoids "conda run" activation bugs in conda ≥26.
# Bandage is a Qt GUI program; the offscreen platform lets it answer --version on headless WSL.
ENV_PREFIX=""
in_env() { PATH="$ENV_PREFIX/bin:$PATH" CONDA_PREFIX="$ENV_PREFIX" QT_QPA_PLATFORM=offscreen "$@"; }

FAILED=0

# ── Preflight ───────────────────────────────────────────────────────
if ! command -v conda &>/dev/null; then
    fail "conda not found. Install Miniforge first (see docs/Drafts/Mine/Bioinformatics_Workshop_Setup_Instructions.md)."
    exit 1
fi
info "conda $(conda --version 2>&1 | awk '{print $2}') found"

PLATFORM="$(uname -s)-$(uname -m)"
info "Platform: $PLATFORM"

case "$(uname -s)" in
    Linux|Darwin) ;;
    *)
        fail "bioconda packages require Linux or macOS (detected: $(uname -s))."
        fail "Run this script inside WSL Ubuntu, not Windows Miniforge Prompt."
        exit 1 ;;
esac

# ── Create environment ──────────────────────────────────────────────
# conda, not mamba: conda ≥23.10 already solves with libmamba, and mamba 2.8's sharded index
# wrongly reported conda-forge packages (isa-l, libdeflate) as missing for this package set.
# --override-channels keeps a "defaults" entry in ~/.condarc out of the solve.
CONDA_CREATE=(conda create -y -n "$ENV_NAME" --override-channels --strict-channel-priority
              -c "$CHANNEL_MAIN" -c "$CHANNEL_BIO" python="$PYTHON_VERSION")
ANNOT_INSTALLED=true

echo ""
echo "━━━ $ENV_NAME ━━━"

if env_exists "$ENV_NAME"; then
    warn "$ENV_NAME already exists — skipping creation"
elif $DRY_RUN; then
    echo "  [dry-run] ${CONDA_CREATE[*]} $CORE_PACKAGES $ANNOT_PACKAGES"
    echo "  [dry-run] if that fails: ${CONDA_CREATE[*]} $CORE_PACKAGES"
else
    echo "Creating $ENV_NAME (several minutes of downloading is normal) ..."
    if "${CONDA_CREATE[@]}" $CORE_PACKAGES $ANNOT_PACKAGES; then
        info "$ENV_NAME created"
    else
        fail "Full install failed — retrying with the core tools only"
        conda env remove -y -n "$ENV_NAME" &>/dev/null || true   # clear any half-built env
        if ! "${CONDA_CREATE[@]}" $CORE_PACKAGES; then
            fail "Core install failed too — save the error above for the trainer, or use the workshop server"
            exit 1
        fi
        ANNOT_INSTALLED=false
        warn "$ENV_NAME created WITHOUT $ANNOT_PACKAGES — use the workshop server for annotation and AMR screening"
    fi
fi

# ── Verify ──────────────────────────────────────────────────────────
# label|command[|fallback command]
TOOLS=(
    "FastQC|fastqc --version"
    "NanoStat|NanoStat --version"
    "fastp|fastp --version"
    "fastplong|fastplong --version"
    "metaSPAdes|metaspades.py --version|spades.py --version"
    "MetaFlye|flye --version"
    "metaQUAST|metaquast --version|metaquast.py --version"
    "Bandage|Bandage --version"
    "Bakta|bakta --version"
    "AMRFinderPlus|amrfinder --version"
    "ABRicate|abricate --version"
)

if ! $DRY_RUN; then
    echo "Verifying $ENV_NAME ..."
    ENV_PREFIX="$(conda info --envs 2>/dev/null | awk -v name="$ENV_NAME" '$1==name {print $NF}')"
    if [[ -z "$ENV_PREFIX" ]]; then
        fail "Could not locate $ENV_NAME prefix"
        FAILED=$((FAILED + 1))
    else
        for entry in "${TOOLS[@]}"; do
            IFS='|' read -ra parts <<< "$entry"
            ok=false
            for cmd in "${parts[@]:1}"; do
                read -ra argv <<< "$cmd"
                if out="$(in_env "${argv[@]}" 2>&1)"; then
                    printf "  %-14s %s\n" "${parts[0]}" "${out%%$'\n'*}"
                    ok=true
                    break
                fi
            done
            if ! $ok; then
                fail "${parts[0]} did not answer (${parts[1]})"
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

# ── Databases (--with-db) ───────────────────────────────────────────
# The tools and their databases are separate downloads; a tool with no database fails only when run.
if $WITH_DB; then
    echo ""
    echo "━━━ Databases → $DB_DIR ━━━"
    if $DRY_RUN; then
        echo "  [dry-run] bakta_db download --output $DB_DIR --type light"
        echo "  [dry-run] amrfinder -u"
        echo "  [dry-run] abricate --setupdb"
    elif [[ -z "$ENV_PREFIX" ]] || ! $ANNOT_INSTALLED; then
        warn "annotation tools not installed — skipping databases"
    else
        mkdir -p "$DB_DIR"
        # light, not full: the full Bakta database runs to tens of gigabytes
        if [[ -d "$DB_DIR/db-light" ]]; then
            warn "Bakta database already at $DB_DIR/db-light — skipping download"
        elif in_env bakta_db download --output "$DB_DIR" --type light; then
            info "Bakta database: $DB_DIR/db-light"
        else
            fail "bakta_db download failed"
            FAILED=$((FAILED + 1))
        fi

        if in_env amrfinder -u; then
            info "AMRFinderPlus database up to date"
        else
            fail "amrfinder -u failed"
            FAILED=$((FAILED + 1))
        fi

        # VFDB ships with ABRicate; --setupdb indexes it alongside the bundled AMR databases
        if in_env abricate --setupdb && [[ "$(in_env abricate --list)" == *vfdb* ]]; then
            info "ABRicate databases indexed (VFDB included)"
        else
            fail "abricate --setupdb failed or VFDB missing"
            FAILED=$((FAILED + 1))
        fi
    fi
fi

# ── Summary ─────────────────────────────────────────────────────────
echo ""
echo "━━━ Summary ━━━"
conda info --envs 2>/dev/null | grep -E "$ENV_NAME" || true
echo ""
if [[ $FAILED -gt 0 ]]; then
    fail "$FAILED check(s) failed. To rebuild: conda env remove -n $ENV_NAME, then rerun this script"
fi
echo "Activate with:   conda activate $ENV_NAME"
if $WITH_DB; then
    echo "Bakta --db:      $DB_DIR/db-light"
else
    echo "Databases:       rerun with --with-db (Bakta light, AMRFinderPlus, ABRicate/VFDB)"
fi
echo "Record versions: conda list -n $ENV_NAME > ${ENV_NAME}_versions.txt"
echo "OPERA-MS:        not installed — build from source only if the trainer asks"
echo "Workshop data:   cd $(cd "$(dirname "$0")" && pwd)/Datasets"

exit $(( FAILED > 0 ))
