#!/usr/bin/env bash
# ============================================================================
# Run the whole protein modeling pipeline, in order
# ============================================================================
# This picks up where the metagenomics pipeline stopped: programs 01 to 11
# assembled and annotated the reads, and programs A to J match the annotation
# against the reference genes, pull out the metagenome's own copy of each one,
# and take it to a structure and a docking score.
#
# Each program is also a script you can run on its own, and on the day you
# should: read the output of one before starting the next. This exists for the
# unattended case — you have been through it once and want the whole thing to
# run while you do something else.
#
# It stops at the first program that fails.
#
# The pipeline has a gap in the middle, and it is not a bug. Programs D and E
# prepare submissions for SWISS-MODEL and AlphaFold3, which are web services:
# nothing comes back until you upload the files and download the results into
# the models/ folder beside each set of uploads, in the folder for that
# organism:
#
#   WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/D_SWISS_MODEL_Inputs/models/<organism>/
#   WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/E_AlphaFold3_Inputs/alphafoldserver/models/<organism>/
#
# Programs F to I skip themselves, without failing, until models are there. So
# a first run goes A to E, then stops being useful; run it again once the
# models have arrived and F to J pick up. Program F collects the downloads
# into WORKING_FOLDER/RESULTS/protein_modeling/Models/, which is derived and yours to delete.
#
# Needs:  conda activate protein_modeling
# Run from the workshop folder:  bash run_protein_modeling.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
# Comment out a line to skip that program.
PROGRAMS="
A_gene_catalog_blastn
B_gene_protein_extraction_bakta
C_protein_prep_python
D_swissmodel_inputs
E_alphafold3_inputs
F_model_qc_pymol
G_receptor_prep_meeko
H_ligand_prep_rdkit_meeko
I_docking_autodock_vina
J_report_python
"
# ============================================================

for PROGRAM in $PROGRAMS; do
    echo
    echo "=============================================================="
    echo "  $PROGRAM"
    echo "=============================================================="
    bash "$PROGRAM.sh"
done
