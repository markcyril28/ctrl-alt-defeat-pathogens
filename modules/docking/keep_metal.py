#!/usr/bin/env python3
"""Rebuild receptor keeping the catalytic metal in place.

A zinc metalloprotease without its zinc is not the enzyme you meant to dock
into. This script strips solvent and hydrogens but preserves the named metal.

Env vars:
    OW_IN    Input PDB path (relative to WORKING_FOLDER/INPUT_DATASETS/from_Database/)
    OW_OUT   Output PDB path (absolute)
    OW_METAL Residue name of the metal to keep (e.g. ZN)
"""
import os
from pymol import cmd

metal = os.environ["OW_METAL"]
cmd.load(os.environ["OW_IN"], "rec")
cmd.remove("rec and solvent")
cmd.remove("rec and hydro")
count = cmd.count_atoms("rec and resn {}".format(metal))
cmd.save(os.environ["OW_OUT"], "rec")
print("receptor keeps its cofactor: {} x {} ({} heavy atoms total)".format(
    count, metal, cmd.count_atoms("rec")))
if count == 0:
    raise SystemExit("ERROR: no {} left in the receptor - check the pre-selection".format(metal))
