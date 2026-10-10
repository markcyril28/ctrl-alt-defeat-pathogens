#!/usr/bin/env python3
"""Strip solvent and hydrogens from a receptor PDB.

Env vars:
    OW_IN  Input PDB path (relative to WORKING_FOLDER/INPUT_DATASETS/from_Database/)
    OW_OUT Output PDB path (absolute)
"""
import os
from pymol import cmd

cmd.load(os.environ["OW_IN"], "rec")
cmd.remove("rec and solvent")
cmd.remove("rec and hydro")
cmd.save(os.environ["OW_OUT"], "rec")
print("receptor: {} heavy atoms".format(cmd.count_atoms("rec")))
