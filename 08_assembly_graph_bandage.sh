#!/usr/bin/env bash
# ============================================================================
# Program 08: Assembly graph — Bandage
# ============================================================================
# An assembly graph shows what a contig list cannot: where the assembler had a
# choice. Branches, dead ends, loops and disconnected pieces are all visible
# in the graph and invisible in N50. Worksheet section 10.
#
# This renders each graph to a PNG so you have a record of it. The PNG is not
# the exercise. Open the graph in Bandage yourself and look at it:
#
#   Bandage load WORKING_FOLDER/RESULTS/metagenomics/Assemblies/metaspades.gfa
#
# OPERA-MS writes no graph file, so only the two single-technology assemblies
# appear here.
#
# Needs:  conda activate meta_env
# Run from the workshop folder:  bash 08_assembly_graph_bandage.sh
# ============================================================================

set -euo pipefail   # stop on an error, on an unset variable, and on a failed pipe

# ========================= SETTINGS =========================
HEIGHT=1200        # height of each PNG, in pixels

ASSEMBLIES="WORKING_FOLDER/RESULTS/metagenomics/Assemblies"
OUT="WORKING_FOLDER/RESULTS/metagenomics/08_Assembly_Graph"
# ============================================================

mkdir -p "$OUT"

for GFA in "$ASSEMBLIES"/*.gfa; do
    NAME=$(basename "$GFA" .gfa)

    # Bandage is a window program. On a machine with no screen it still works
    # for this, because QT_QPA_PLATFORM tells Qt to draw off-screen.
    QT_QPA_PLATFORM=offscreen Bandage image "$GFA" "$OUT/$NAME.png" --height "$HEIGHT"
done
