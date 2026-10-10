#!/usr/bin/env python3
"""Step 8d — Superpose each model on its crystal structure (super + align).

Env vars:
    OW_CONFIG_JSON (optional) JSON config for single-species mode
"""
import os
import json
from pymol import cmd

_ALL_PAIRS = [
    # label,        model path,                                                      model selection,  crystal path,                                                crystal selection
    ("TNT domain",  "1_Mycobacterium_tuberculosis/MODEL_O05442_CpnT_full_AF.pdb",   "resi 651-846",   "1_Mycobacterium_tuberculosis/PDB_4QLP_TNT_immunity.pdb",     "chain B and polymer"),
    ("FimH lectin", "3_Klebsiella_pneumoniae/MODEL_A0A0H3H2I8_FimH_full_AF.pdb",    "polymer",        "3_Klebsiella_pneumoniae/PDB_9AT9_FimH_lectin_mannose.pdb",  "polymer"),
]

_cfg_json = os.environ.get("OW_CONFIG_JSON", "")
if _cfg_json:
    _cfg = json.load(open(_cfg_json))
    _sd = _cfg["species"]["dir"]
    PAIRS = [
        (v["label"], "{}/{}".format(_sd, v["model"]), v["model_sel"],
         "{}/{}".format(_sd, v["crystal"]), v["crystal_sel"])
        for v in _cfg.get("validate", [])
    ]
else:
    PAIRS = _ALL_PAIRS

print("{:14s} {:>12s} {:>7s} {:>12s} {:>7s}".format(
    "comparison", "super RMSD", "atoms", "align RMSD", "atoms"))
print("-" * 58)
for label, model_path, model_sel, xtal_file, xtal_sel in PAIRS:
    if not (os.path.exists(model_path) and os.path.exists(xtal_file)):
        print("SKIP {}: run the fetch step first".format(label))
        continue
    cmd.delete("all")
    cmd.load(model_path, "pred_src")
    cmd.load(xtal_file, "xtal_src")
    cmd.create("pred1", "pred_src and ({})".format(model_sel))
    cmd.create("xtal1", "xtal_src and ({})".format(xtal_sel))
    sup_rms, sup_atoms = cmd.super("pred1", "xtal1")[:2]
    cmd.create("pred2", "pred_src and ({})".format(model_sel))
    cmd.create("xtal2", "xtal_src and ({})".format(xtal_sel))
    ali_rms, ali_atoms = cmd.align("pred2", "xtal2")[:2]
    print("{:14s} {:12.2f} {:7d} {:12.2f} {:7d}".format(
        label, sup_rms, sup_atoms, ali_rms, ali_atoms))

print("\nWhere a crystal structure of your site exists, dock into the crystal.")
print("Use models only for the sites where nothing else exists.")
