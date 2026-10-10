#!/usr/bin/env python3
"""Step 4 — Load the dataset into one combined PyMOL session.

Env vars:
    OW_RECEPTORS   Output directory for the session file
    OW_CONFIG_JSON (optional) JSON config for single-species mode
"""
import os
import glob
import inspect
import json
from pymol import cmd

# toxin_load.py sits next to this file. Under `pymol script.py` __file__ is
# PyMOL's own, so take this script's path from its code object instead.
_here = os.path.abspath(inspect.getfile(inspect.currentframe()))
TOXIN_LOAD = os.path.join(os.path.dirname(_here), "toxin_load.py")

out =os.path.join(os.environ["OW_RECEPTORS"], "toxin_dataset.pse")
_cfg_json = os.environ.get("OW_CONFIG_JSON", "")
if _cfg_json:
    _cfg = json.load(open(_cfg_json))
    _species_dir = _cfg["species"]["dir"]
    for pdb in sorted(glob.glob(os.path.join(_species_dir, "*.pdb"))):
        name = os.path.splitext(os.path.basename(pdb))[0]
        cmd.load(pdb, name)
else:
    cmd.do("run " + TOXIN_LOAD)
    cmd.do("load_all_pdbs()")

cmd.do("run " + TOXIN_LOAD)
cmd.do("color_by_chain_all()")
loaded = cmd.get_object_list()
cmd.save(out)
print("loaded {} structures: {}".format(len(loaded), ", ".join(sorted(loaded))))
print("session saved to", out)
