#!/usr/bin/env bash
# ============================================================================
# Run the whole metagenomics pipeline, in order
# ============================================================================
# Each program is also a script you can run on its own, and on the day you
# should: read the output of one before starting the next. This exists for the
# unattended case — you have been through it once and want the whole thing to
# run while you do something else.
#
# It stops at the first program that fails.
#
# Program 02 runs on the raw reads. For the second pass, set WHICH="cleaned"
# in 02_read_qc.sh and run that one again afterwards.
#
# Needs:  conda activate metagenomics_env
# Run from the workshop folder:  bash run_metagenomics.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
# Comment out a line to skip that program.
PROGRAMS="
01_data_inventory
02_read_qc
03_read_cleaning
04_assembly_short
05_assembly_long
06_assembly_hybrid
07_assembly_evaluation
08_assembly_graph
09_annotation
10_amr_virulence_screening
11_report
"
# ============================================================

for PROGRAM in $PROGRAMS; do
    echo
    echo "=============================================================="
    echo "  $PROGRAM"
    echo "=============================================================="
    bash "$PROGRAM.sh"
done
