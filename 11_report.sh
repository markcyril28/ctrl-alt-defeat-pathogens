#!/usr/bin/env bash
# ============================================================================
# Program 11: Report — the worksheet's tables, and every version
# ============================================================================
# Each tool in this pipeline wrote its own report in its own format. This
# gathers the numbers the worksheet asks for into three tables, and records
# the version of every tool and database you used.
#
# The versions are not an afterthought. Software and databases change, and the
# same assembly screened a year later can give a different answer. A result
# without its versions cannot be reproduced.
#
# Nothing here recalculates anything, so it is always safe to run again.
#
# Needs:  conda activate metagenomics_env
# Run from the workshop folder:  bash 11_report.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
RESULTS="RESULTS/metagenomics"
OUT="RESULTS/metagenomics/11_Report"
HELPER="modules/metagenomics/collect_stats.py"
# ============================================================

mkdir -p "$OUT"

# reads      -> worksheet 2.4 and 5   read counts, lengths, quality, raw vs cleaned
# annotation -> worksheet 13          Bakta feature counts per assembly
# screening  -> worksheet 14          AMR and virulence hits per assembly
for TABLE in reads annotation screening; do
    python3 "$HELPER" "$TABLE" "$RESULTS" > "$OUT/$TABLE.tsv"
    echo "--- $TABLE ---"
    column -t -s $'\t' "$OUT/$TABLE.tsv"
    echo
done

echo "recorded $(date)" > "$OUT/versions.txt"
echo "machine  $(uname -s) $(uname -m)" >> "$OUT/versions.txt"

# Each tool's version, and its database version where it has one. A tool that
# is not installed is recorded as missing rather than stopping the report.
for TOOL in "conda --version" \
            "fastqc --version" \
            "NanoStat --version" \
            "fastp --version" \
            "fastplong --version" \
            "spades.py --version" \
            "flye --version" \
            "metaquast --version" \
            "Bandage --version" \
            "bakta --version" \
            "bakta_db list" \
            "amrfinder --version" \
            "amrfinder --database_version" \
            "abricate --version" \
            "abricate --list"; do
    echo "===== $TOOL" >> "$OUT/versions.txt"
    $TOOL >> "$OUT/versions.txt" 2>&1 || echo "(not available)" >> "$OUT/versions.txt"
done
