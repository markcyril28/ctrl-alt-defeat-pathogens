#!/usr/bin/env python3
"""Restrict receptor to chosen chains/selection before splitting.

Env vars:
    OW_IN   Input PDB path (relative to WORKING_FOLDER/INPUT_DATASETS/from_Database/)
    OW_OUT  Output PDB path (absolute)
    OW_KEEP PyMOL selection expression for atoms to keep
"""
import os
from pymol import cmd

cmd.load(os.environ["OW_IN"], "src")
before = cmd.count_atoms("src")
cmd.create("sub", "src and ({})".format(os.environ["OW_KEEP"]))
cmd.remove("sub and solvent")
cmd.save(os.environ["OW_OUT"], "sub")
print("subset: {} -> {} atoms kept by `{}`".format(
    before, cmd.count_atoms("sub"), os.environ["OW_KEEP"]))
