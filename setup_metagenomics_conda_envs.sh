#!/usr/bin/env bash
# setup_metagenomics_conda_envs.sh — Create the Ologist Workshop metagenomics conda environments
# Run from: Ubuntu/WSL terminal  or  Mac Terminal
# Usage:    bash setup_metagenomics_conda_envs.sh
#           bash setup_metagenomics_conda_envs.sh --dry-run        # show what would be created
#           bash setup_metagenomics_conda_envs.sh --no-db          # skip the database downloads
#           bash setup_metagenomics_conda_envs.sh --db-dir=/path/to/databases
#           bash setup_metagenomics_conda_envs.sh --no-mamba       # install with conda only
#
# Environments created (names set at the top of the script):
#   meta_env          — python 3.10, FastQC, NanoStat, fastp, fastplong, SPAdes (metaSPAdes),
#                       Flye (MetaFlye), QUAST (metaQUAST), Bandage
#   meta_annot_env    — python 3.10, Bakta, AMRFinderPlus
#   meta_abricate_env — ABRicate (VFDB and the other screening databases)
# Only meta_env is activated: its bakta, amrfinder and abricate commands are wrappers
# that run the programs in the other two environments.
#
# Databases, downloaded by default into ~/workshop/databases: AMRFinderPlus, ABRicate's VFDB, and
# Bakta light (1.3 GB from Zenodo — an interrupted download resumes where it stopped on the next run).
#
# Not installed: OPERA-MS — not a conda package; build it from source only if the trainer asks.
#
# Idempotent: skips an environment that already exists, and the Bakta database if present.
# Exceptions: annotation tools that older runs put in the core env are removed from it,
# and an older-schema Bakta database is replaced.
# Kept separate from protein_modeling: the assemblers' dependencies would downgrade its packages.
# bioconda packages require Linux or macOS.
# Every run is also written to logs/setup_metagenomics_conda_envs_<date>_<time>.log next to this script.

#set -euo pipefail

# ── Conda environments ──────────────────────────────────────────────
CORE_ENV="meta_env"                # QC → assembly → evaluation; the env students activate
ANNOT_ENV="meta_annot_env"         # Bakta and AMRFinderPlus, reached through wrappers in CORE_ENV
ABRICATE_ENV="meta_abricate_env"   # ABRicate, reached the same way
# ────────────────────────────────────────────────────────────────────

# ── Configuration ──────────────────────────────────────
# mamba installs faster and is used when present; any env it cannot build is retried with conda
USE_MAMBA=true
PYTHON_VERSION="3.10"
CHANNEL_MAIN="conda-forge"
CHANNEL_BIO="bioconda"
# QC → cleaning → assembly → evaluation → graph. metaQUAST ships inside the quast package.
CORE_PACKAGES="fastqc nanostat fastp fastplong spades flye quast bandage"
# Why three envs — each pair clashes on a shared library, so no single env can hold them:
#   QUAST:    every build pins BLAST below 2.17 (or zlib below 1.3); Bakta 1.12 needs BLAST 2.17 and DIAMOND 2.2.
#   ABRicate: its BioPerl pulls in perl-bio-samtools, which needs zlib below 1.3 or samtools 0.1.19
#             (OpenSSL 1.1); AMRFinderPlus 4.2 and Bakta 1.12 need zlib 1.3.2 and OpenSSL 3.
# The floors keep the solver off Bakta 1.9, which caps AMRFinderPlus at 3.12: its database froze at
# 2024-07-22.1, and "amrfinder -u" only warns that newer databases need AMRFinderPlus 4.2.
AMR_MIN="4.2"
ANNOT_PACKAGES="bakta>=1.12 ncbi-amrfinderplus>=$AMR_MIN"
# ABRicate is Perl: no python spec, so the solver is free to pick what its BioPerl needs
ABRICATE_PACKAGES="abricate"
DB_DIR="$HOME/workshop/databases"
# ────────────────────────────────────────────────────────────────────

DRY_RUN=false
WITH_DB=true

for arg in "$@"; do
    case "$arg" in
        --dry-run)   DRY_RUN=true ;;
        --with-db)   WITH_DB=true ;;    # the default; still accepted from older instructions
        --no-db)     WITH_DB=false ;;
        --db-dir=*)  DB_DIR="${arg#--db-dir=}" ;;
        --no-mamba)  USE_MAMBA=false ;;
        -h|--help)
            sed -n '2,28p' "$0"
            exit 0 ;;
        *)
            echo "Unknown option: $arg"; exit 1 ;;
    esac
done
DB_DIR="${DB_DIR/#\~/$HOME}"   # --db-dir=~/x is not tilde-expanded by the shell

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
# AMRFinderPlus looks for its database under $CONDA_PREFIX, so that follows the env too.
# Bandage is a Qt GUI program; the offscreen platform lets it answer --version on headless WSL.
in_env() { local prefix="$1"; shift; PATH="$prefix/bin:$PATH" CONDA_PREFIX="$prefix" QT_QPA_PLATFORM=offscreen "$@"; }

# True if the AMRFinderPlus reached from env prefix $1 is older than $AMR_MIN (false if it is missing — verify reports that)
amr_outdated() {
    local have have_major have_minor need_major need_minor
    have="$(in_env "$1" amrfinder --version 2>/dev/null)" || return 1
    IFS=. read -r have_major have_minor _ <<< "$have"
    IFS=. read -r need_major need_minor <<< "$AMR_MIN"
    (( have_major < need_major || (have_major == need_major && have_minor < need_minor) ))
}

# True if the Bakta database in $1 has the schema the installed Bakta reads (bakta_db list shows it)
bakta_db_current() {
    in_env "$ANNOT_PREFIX" python - "$1/version.json" 2>/dev/null <<'EOF'
import bakta, json, sys
sys.exit(json.load(open(sys.argv[1]))["major"] != bakta.__db_schema_version__)
EOF
}

# Prints "<url> <md5>" of the newest light database this Bakta reads — the lookup bakta_db download makes
bakta_db_source() {
    in_env "$ANNOT_PREFIX" python - <<'EOF'
import bakta
from bakta.db import fetch_db_versions
v = max((v for v in fetch_db_versions() if v["major"] == bakta.__db_schema_version__), key=lambda v: v["minor"])
print(f"https://zenodo.org/records/{v['record']}/files/db-light.tar.xz", v["md5-light"])
EOF
}

# MD5 of file $1, computed the way bakta_db does
file_md5() {
    in_env "$ANNOT_PREFIX" python -c 'import sys, pathlib, bakta.db; print(bakta.db.calc_md5_sum(pathlib.Path(sys.argv[1])))' "$1"
}

size_of() { if [[ -f "$1" ]]; then wc -c < "$1" | tr -d ' '; else echo 0; fi; }

# curl URL $1 into file $2, continuing from whatever part of $2 is already there, and retrying a dropped
# or stalled connection. Gives up after 5 tries in a row that add nothing.
download_resumable() {
    local url="$1" out="$2" stalls=0 before code rc
    while true; do
        before="$(size_of "$out")"
        # -f keeps an HTTP error page out of the file; the speed limit turns a stalled connection into a retry
        code="$(curl -fL -C - --connect-timeout 30 --speed-limit 1024 --speed-time 120 \
                     -o "$out" -w '%{http_code}' "$url")" && return 0
        rc=$?
        # HTTP 416 or curl exit 33: the server will not continue this partial file — start it over
        if [[ $code == 416 || $rc == 33 ]]; then
            rm -f "$out"
        fi
        if (( $(size_of "$out") > before )); then stalls=0; else stalls=$((stalls + 1)); fi
        (( stalls < 5 )) || return 1
        echo "Download interrupted at $(( $(size_of "$out") / 1048576 )) MB (curl exit $rc, HTTP $code) — resuming in 30 s"
        sleep 30
    done
}

# Version of the AMRFinderPlus database in folder $1 — the folder its "latest" link points to (empty if none)
amr_db_version() { local link; link="$(readlink "$1/latest" 2>/dev/null)" && basename "$link"; return 0; }

# Install the Bakta light database at $DB_DIR/db-light (light, not full: the full one runs to tens of gigabytes).
# bakta_db download cannot resume, and Zenodo serves the 1.3 GB tarball slowly and drops long transfers — one run
# lost the connection 50 MB in. So curl fetches it into a staging folder that a rerun resumes from, and
# bakta_db install unpacks it there, checks it and updates the AMRFinderPlus database bundled inside.
# The old database is swapped out only after that succeeds.
install_bakta_db() {
    local stage="$DB_DIR/.bakta_db_download" src url md5
    local tarball="$stage/db-light.tar.xz"
    if ! src="$(bakta_db_source)"; then
        fail "Could not look up the current Bakta database — check the internet connection, then rerun"
        return 1
    fi
    read -r url md5 <<< "$src"
    mkdir -p "$stage"
    if [[ -f "$tarball" && "$(file_md5 "$tarball")" == "$md5" ]]; then
        echo "Bakta database already downloaded: $tarball"
    else
        echo "Downloading the Bakta light database (1.3 GB) from Zenodo — over an hour on a slow connection ..."
        if [[ -f "$tarball" ]]; then
            echo "Resuming from $(( $(size_of "$tarball") / 1048576 )) MB"
        fi
        if ! download_resumable "$url" "$tarball"; then
            fail "Bakta database download keeps failing — rerun later; it resumes from $(( $(size_of "$tarball") / 1048576 )) MB"
            return 1
        fi
        if [[ "$(file_md5 "$tarball")" != "$md5" ]]; then
            rm -f "$tarball"
            fail "Bakta database download is corrupt (MD5 mismatch) — deleted; rerun to download it again"
            return 1
        fi
    fi
    rm -rf "${stage:?}/db-light"   # half-unpacked by an interrupted run
    if ! in_env "$ANNOT_PREFIX" bakta_db install --db-file "$tarball" --output "$stage"; then
        fail "bakta_db install failed — the download is kept, so a rerun only repeats the install"
        return 1
    fi
    rm -rf "${DB_DIR:?}/db-light"
    mv "$stage/db-light" "$DB_DIR/db-light"
    rm -rf "${stage:?}"
}

# --override-channels keeps a "defaults" entry in ~/.condarc out of the solve (conda and mamba both read it).
CHANNELS=(--override-channels -c "$CHANNEL_MAIN" -c "$CHANNEL_BIO")

# Create env $1 from the package specs that follow, unless it already exists. Returns 1 if it cannot be built.
# mamba goes first when $INSTALLER is mamba. Both solve with libmamba, but mamba 2.8's sharded index once
# reported conda-forge packages (isa-l, libdeflate) as missing for this package set, so a mamba failure
# is retried with conda rather than reported.
make_env() {
    local name="$1"
    shift
    echo ""
    echo "━━━ $name ━━━"
    if env_exists "$name"; then
        warn "$name already exists — skipping creation"
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

# Put command $1 from env $2 (at prefix $3) in $CORE_ENV/bin as a wrapper that runs it on that env's PATH,
# so each tool finds its own BLAST, Perl and libraries rather than the core env's. Never overwrites a conda file.
WRAPPER_TAG="wrapper written by setup_metagenomics_conda_envs.sh"
write_wrapper() {
    local target="$CORE_PREFIX/bin/$1"
    if [[ -e "$target" ]] && ! grep -qF "$WRAPPER_TAG" "$target"; then
        fail "$target was installed by conda, not this script — left alone"
        FAILED=$((FAILED + 1))
        return 1
    fi
    printf '#!/usr/bin/env bash\n# %s — runs %s from %s\nexport PATH=%q:"$PATH" CONDA_PREFIX=%q\nexec %q "$@"\n' \
        "$WRAPPER_TAG" "$1" "$2" "$3/bin" "$3" "$3/bin/$1" > "$target"
    chmod +x "$target"
}

# Wrap every command in env $1 (at prefix $2) whose name starts with one of the words that follow
wrap_tools() {
    local name="$1" prefix="$2" stem exe wrapped=""
    shift 2
    for stem in "$@"; do
        for exe in "$prefix/bin/$stem"*; do
            [[ -f "$exe" && -x "$exe" ]] || continue
            if write_wrapper "$(basename "$exe")" "$name" "$prefix"; then
                wrapped+=" $(basename "$exe")"
            fi
        done
    done
    info "$CORE_ENV runs these from $name:$wrapped"
}

FAILED=0
BAKTA_DB_READY=false
CORE_PREFIX=""
ANNOT_PREFIX=""
ABRICATE_PREFIX=""

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

# ── Core environment ────────────────────────────────────────────────
if ! make_env "$CORE_ENV" python="$PYTHON_VERSION" $CORE_PACKAGES; then
    fail "Every program needs $CORE_ENV — save the error for the trainer, or use the workshop server"
    exit 1
fi
CORE_PREFIX="$(env_prefix "$CORE_ENV")"

# Older runs installed the annotation tools here, where they held Bakta at 1.9 and AMRFinderPlus at 3.12
OLD_ANNOT=""
if [[ -n "$CORE_PREFIX" ]]; then
    for spec in $ANNOT_PACKAGES $ABRICATE_PACKAGES; do
        pkg="${spec%%[<>=]*}"
        if compgen -G "$CORE_PREFIX/conda-meta/$pkg-[0-9]*.json" >/dev/null; then
            OLD_ANNOT+=" $pkg"
        fi
    done
fi
if [[ -z "$OLD_ANNOT" ]]; then
    :
elif $DRY_RUN; then
    echo "  [dry-run] conda remove -y -n $CORE_ENV ${CHANNELS[*]}$OLD_ANNOT"
else
    echo "$CORE_ENV still holds$OLD_ANNOT from an older run — removing them; $ANNOT_ENV and $ABRICATE_ENV replace them"
    if conda remove -y -n "$CORE_ENV" "${CHANNELS[@]}" $OLD_ANNOT; then
        info "Old annotation tools removed from $CORE_ENV"
    else
        fail "Could not remove$OLD_ANNOT from $CORE_ENV — rebuild it: conda env remove -n $CORE_ENV, then rerun"
        FAILED=$((FAILED + 1))
    fi
fi

# ── Annotation environments ─────────────────────────────────────────
if make_env "$ANNOT_ENV" python="$PYTHON_VERSION" $ANNOT_PACKAGES; then
    ANNOT_PREFIX="$(env_prefix "$ANNOT_ENV")"
else
    warn "No Bakta or AMRFinderPlus — use the workshop server for annotation and AMR screening"
    FAILED=$((FAILED + 1))
fi

if make_env "$ABRICATE_ENV" $ABRICATE_PACKAGES; then
    ABRICATE_PREFIX="$(env_prefix "$ABRICATE_ENV")"
else
    warn "No ABRicate — use the workshop server for virulence screening"
    FAILED=$((FAILED + 1))
fi

# ── Wrappers ────────────────────────────────────────────────────────
# Rewritten on every run, so they follow a rebuilt env. Removing an env leaves its wrappers
# pointing at nothing; verify then reports those tools as not answering.
echo ""
if $DRY_RUN; then
    echo "  [dry-run] wrappers in $CORE_ENV/bin for every bakta* and amrfinder* command in $ANNOT_ENV"
    echo "  [dry-run] wrappers in $CORE_ENV/bin for every abricate* command in $ABRICATE_ENV"
elif [[ -n "$CORE_PREFIX" ]]; then
    if [[ -n "$ANNOT_PREFIX" ]]; then
        wrap_tools "$ANNOT_ENV" "$ANNOT_PREFIX" bakta amrfinder
    fi
    if [[ -n "$ABRICATE_PREFIX" ]]; then
        wrap_tools "$ABRICATE_ENV" "$ABRICATE_PREFIX" abricate
    fi
fi

# ── Verify ──────────────────────────────────────────────────────────
# Everything is checked from $CORE_ENV, the way students run it — the annotation tools through their wrappers.
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
    echo "Verifying $CORE_ENV ..."
    VERIFY_FAILED=0
    if [[ -z "$CORE_PREFIX" ]]; then
        fail "Could not locate $CORE_ENV prefix"
        VERIFY_FAILED=$((VERIFY_FAILED + 1))
    else
        for entry in "${TOOLS[@]}"; do
            IFS='|' read -ra parts <<< "$entry"
            ok=false
            for cmd in "${parts[@]:1}"; do
                read -ra argv <<< "$cmd"
                if out="$(in_env "$CORE_PREFIX" "${argv[@]}" 2>&1)"; then
                    printf "  %-14s %s\n" "${parts[0]}" "${out%%$'\n'*}"
                    ok=true
                    break
                fi
            done
            if ! $ok; then
                fail "${parts[0]} did not answer (${parts[1]})"
                VERIFY_FAILED=$((VERIFY_FAILED + 1))
            fi
        done
        if amr_outdated "$CORE_PREFIX"; then
            fail "AMRFinderPlus $(in_env "$CORE_PREFIX" amrfinder --version) is older than $AMR_MIN — its database cannot update past 2024-07-22.1"
            VERIFY_FAILED=$((VERIFY_FAILED + 1))
        fi
        if [[ $VERIFY_FAILED -eq 0 ]]; then
            info "$CORE_ENV verified"
        else
            fail "$VERIFY_FAILED tool(s) failed verification — activate and debug manually"
        fi
    fi
    FAILED=$((FAILED + VERIFY_FAILED))
fi

# ── Databases (on by default; --no-db skips them) ───────────────────
# The tools and their databases are separate downloads; a tool with no database fails only when run.
if $WITH_DB; then
    echo ""
    echo "━━━ Databases → $DB_DIR ━━━"
    if $DRY_RUN; then
        echo "  [dry-run] amrfinder -u"
        echo "  [dry-run] curl the Bakta light database from Zenodo (resumable), then bakta_db install it as $DB_DIR/db-light"
        echo "  [dry-run] abricate --setupdb"
    else
        if [[ -z "$ANNOT_PREFIX" ]]; then
            warn "$ANNOT_ENV not installed — skipping the Bakta and AMRFinderPlus databases"
        else
            mkdir -p "$DB_DIR"
            if in_env "$ANNOT_PREFIX" amrfinder -u; then
                info "AMRFinderPlus database up to date"
            else
                fail "amrfinder -u failed"
                FAILED=$((FAILED + 1))
            fi

            if bakta_db_current "$DB_DIR/db-light"; then
                warn "Bakta database already at $DB_DIR/db-light — skipping download"
                BAKTA_DB_READY=true
                # Bakta runs AMRFinderPlus on its own copy of the AMRFinderPlus database, which bakta_db install
                # updated; it is updated again (Bakta's documented command) once NCBI has published a newer one.
                BUNDLED_AMR="$DB_DIR/db-light/amrfinderplus-db"
                BUNDLED_VERSION="$(amr_db_version "$BUNDLED_AMR")"
                if [[ -n "$BUNDLED_VERSION" && "$BUNDLED_VERSION" == "$(amr_db_version "$ANNOT_PREFIX/share/amrfinderplus/data")" ]]; then
                    info "AMRFinderPlus database inside the Bakta database already at $BUNDLED_VERSION"
                elif in_env "$ANNOT_PREFIX" amrfinder_update --force_update --database "$BUNDLED_AMR"; then
                    info "AMRFinderPlus database inside the Bakta database updated"
                else
                    fail "amrfinder_update failed for $BUNDLED_AMR"
                    FAILED=$((FAILED + 1))
                fi
            else
                if [[ -d "$DB_DIR/db-light" ]]; then
                    echo "Bakta database at $DB_DIR/db-light is incomplete or the wrong schema for this Bakta — replacing it"
                fi
                if install_bakta_db; then
                    info "Bakta database: $DB_DIR/db-light"
                    BAKTA_DB_READY=true
                else
                    FAILED=$((FAILED + 1))
                fi
            fi
        fi

        # VFDB ships with ABRicate; --setupdb indexes it alongside the bundled AMR databases
        if [[ -z "$ABRICATE_PREFIX" ]]; then
            warn "$ABRICATE_ENV not installed — skipping the ABRicate databases"
        elif in_env "$ABRICATE_PREFIX" abricate --setupdb && [[ "$(in_env "$ABRICATE_PREFIX" abricate --list)" == *vfdb* ]]; then
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
conda info --envs 2>/dev/null | grep -E "^($CORE_ENV|$ANNOT_ENV|$ABRICATE_ENV) " || true
echo ""
if [[ $FAILED -gt 0 ]]; then
    fail "$FAILED check(s) failed — see [FAIL] above. Rerun to retry (downloads resume); rebuild a broken env with: conda env remove -n <name>"
fi
echo "Activate with:   conda activate $CORE_ENV"
echo "                 (bakta and amrfinder run from $ANNOT_ENV, abricate from $ABRICATE_ENV)"
if ! $WITH_DB; then
    echo "Databases:       skipped (--no-db) — rerun without it for Bakta light, AMRFinderPlus, ABRicate/VFDB"
elif $BAKTA_DB_READY || $DRY_RUN; then
    echo "Bakta --db:      $DB_DIR/db-light"
else
    echo "Bakta --db:      not ready — rerun this script; the download resumes where it stopped"
fi
echo "Record versions: conda list -n $CORE_ENV > ${CORE_ENV}_versions.txt"
echo "                 conda list -n $ANNOT_ENV > ${ANNOT_ENV}_versions.txt"
echo "                 conda list -n $ABRICATE_ENV > ${ABRICATE_ENV}_versions.txt"
echo "OPERA-MS:        not installed — build from source only if the trainer asks"
echo "Workshop data:   cd $(cd "$(dirname "$0")" && pwd)/Datasets"
echo "Log:             $SETUP_LOG"

exit $(( FAILED > 0 ))
