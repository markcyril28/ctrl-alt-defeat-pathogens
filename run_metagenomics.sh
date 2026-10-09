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
# in 02_read_qc_fastqc_nanostat.sh and run that one again afterwards.
#
# Needs:  conda activate metagenomics_env
# Run from the workshop folder:  bash run_metagenomics.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
# Comment out a line to skip that program.
PROGRAMS="
01_data_inventory_awk
02_read_qc_fastqc_nanostat
03_read_cleaning_fastp_fastplong
04_assembly_short_metaspades
05_assembly_long_metaflye
06_assembly_hybrid_operams
07_assembly_evaluation_metaquast
08_assembly_graph_bandage
09_annotation_bakta
10_amr_virulence_screening_amrfinderplus_abricate
11_report_python
"
# ============================================================

for PROGRAM in $PROGRAMS; do
    echo
    echo "=============================================================="
    echo "  $PROGRAM"
    echo "=============================================================="
    bash "$PROGRAM.sh"
done
