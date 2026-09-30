#!/usr/bin/env bash
# setup_conda_envs.sh — Create the Ologist Workshop conda environment
# Run from: Ubuntu/WSL terminal  or  Mac Terminal
# Usage:    bash setup_conda_envs.sh
#           bash setup_conda_envs.sh --dry-run     # show what would be created
#
# Environment created:
#   protein_model — python 3.11, biopython, pandas, pymol, NGS tools, docking (Vina)
#
# Idempotent: skips the environment if it already exists.
# Channel priority: conda-forge first, bioconda second (required by bioconda).
# bioconda packages (samtools, bwa, etc.) require Linux or macOS.

set -euo pipefail

# ── Configuration ──────────────────────────────────────
ENV_NAME="protein_modeling"
PYTHON_VERSION="3.11"
CHANNEL_MAIN="conda-forge"
CHANNEL_BIO="bioconda"
PACKAGES="biopython pandas pymol-open-source fastqc seqkit bwa samtools bcftools blast vina meeko rdkit"
# ────────────────────────────────────────────────────────────────────

DRY_RUN=false

for arg in "$@"; do
    case "$arg" in
        --dry-run)   DRY_RUN=true ;;
        -h|--help)
            sed -n '2,12p' "$0"
            exit 0 ;;
        *)
            echo "Unknown option: $arg"; exit 1 ;;
    esac
done

# ── Helpers ─────────────────────────────────────────────────────────
GREEN='\033[0;32m'  YELLOW='\033[0;33m'  RED='\033[0;31m'  NC='\033[0m'

info()  { printf "${GREEN}[OK]${NC}  %s\n" "$*"; }
warn()  { printf "${YELLOW}[SKIP]${NC} %s\n" "$*"; }
fail()  { printf "${RED}[FAIL]${NC} %s\n" "$*"; }

env_exists() {
    conda info --envs 2>/dev/null | grep -qE "^$1 " || \
    conda info --envs 2>/dev/null | grep -qE "/$1\$"
}

# ── Preflight ───────────────────────────────────────────────────────
if ! command -v conda &>/dev/null; then
    fail "conda not found. Install Miniforge first (see docs/Bioinformatics_Workshop_Manual.md)."
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
echo ""
echo "━━━ $ENV_NAME ━━━"

if env_exists "$ENV_NAME"; then
    warn "$ENV_NAME already exists — skipping creation"
else
    if $DRY_RUN; then
        echo "  [dry-run] conda create -n $ENV_NAME -c $CHANNEL_MAIN -c $CHANNEL_BIO python=$PYTHON_VERSION $PACKAGES"
    else
        echo "Creating $ENV_NAME ..."
        conda create -y -n "$ENV_NAME" \
            -c "$CHANNEL_MAIN" -c "$CHANNEL_BIO" \
            python="$PYTHON_VERSION" $PACKAGES
        info "$ENV_NAME created"
    fi
fi

# ── Verify ──────────────────────────────────────────────────────────
if ! $DRY_RUN; then
    echo "Verifying $ENV_NAME ..."
    # Resolve the env prefix — avoids "conda run" activation bugs in conda ≥26
    ENV_PREFIX="$(conda info --envs 2>/dev/null | awk -v name="$ENV_NAME" '$1==name {print $NF}')"
    if [[ -z "$ENV_PREFIX" ]]; then
        fail "Could not locate $ENV_NAME prefix"
    elif "$ENV_PREFIX/bin/python" -c \
            "import Bio, pandas, pymol, meeko; print('biopython', Bio.__version__, 'pandas', pandas.__version__, 'pymol', pymol.cmd.get_version()[0])" \
         && "$ENV_PREFIX/bin/samtools" --version | head -1 \
         && "$ENV_PREFIX/bin/vina" --version \
         && "$ENV_PREFIX/bin/seqkit" version; then
        info "$ENV_NAME verified"
    else
        fail "$ENV_NAME verification failed — activate and debug manually"
    fi
fi

# ── Summary ─────────────────────────────────────────────────────────
echo ""
echo "━━━ Summary ━━━"
conda info --envs 2>/dev/null | grep -E "$ENV_NAME" || true
echo ""
echo "Activate with:  conda activate $ENV_NAME"
echo "Workshop data:  cd $(cd "$(dirname "$0")" && pwd)/Datasets"
