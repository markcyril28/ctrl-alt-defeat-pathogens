#!/usr/bin/env python3
"""Compute the centroid of the chain A / chain B contact surface.

The contact surface between the effector and its chaperone, not either
chain's own centre of mass — a box on the whole complex would search
mostly solvent.

Env vars:
    OW_IN  Input PDB path (relative to WORKING_FOLDER/INPUT_DATASETS/from_Database/)

Prints a CENTRE line that the calling bash script parses:
    CENTRE <x> <y> <z>
"""
import os
from pymol import cmd

cmd.load(os.environ["OW_IN"], "src")
cmd.select("contacts", "(src and chain A and polymer and not hydro) within 4.5 of "
                       "(src and chain B and polymer and not hydro)")
n = cmd.count_atoms("contacts")
if n == 0:
    raise SystemExit("ERROR: no A/B interface contacts found")
coords = cmd.get_coords("contacts")
x, y, z = (float(v) for v in coords.mean(axis=0))
print("interface contact atoms: {}".format(n))
print("CENTRE {:.2f} {:.2f} {:.2f}".format(x, y, z))
