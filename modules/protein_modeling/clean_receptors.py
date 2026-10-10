#!/usr/bin/env python3
"""Step 6 — Export docking-ready receptors, one per target recipe.

Env vars:
    OW_RECEPTORS   Output directory for cleaned receptor PDBs
    OW_CONFIG_JSON (optional) JSON config for single-species mode
"""
import os
import json
from pymol import cmd

OUT = os.environ["OW_RECEPTORS"]

_ALL_RECIPES = [
    # name,                pdb path (relative to WORKING_FOLDER/INPUT_DATASETS/from_Database/),  keep selection,                              drop resn,         Zn, outfile
    ("TNT",                "1_Mycobacterium_tuberculosis/PDB_4QLP_TNT_immunity.pdb",             "chain B and polymer",                       "",                0, "TNT_receptor_clean.pdb"),
    ("BoNT/A holotoxin",   "2_Clostridium_botulinum/PDB_3BTA_BoNT_A_holotoxin.pdb",              "(chain A and polymer) or resn ZN",          "",                1, "BONT_A_holo_clean.pdb"),
    ("BoNT/A LC",          "2_Clostridium_botulinum/PDB_1XTG_BoNT_A_LC_SNAP25.pdb",              "(chain A and polymer) or resn ZN",          "CL",              1, "BONT_A_LC_clean.pdb"),
    ("BoNT/B catalytic",   "2_Clostridium_botulinum/PDB_1EPW_BoNT_B.pdb",                        "polymer or resn ZN",                        "SO4",             1, "BONT_B_catalytic_clean.pdb"),
    ("BoNT/B receptor-bd", "2_Clostridium_botulinum/PDB_1I1E_BoNT_B_doxorubicin.pdb",            "polymer",                                   "DM2+SO4",         0, "BONT_B_HC_clean.pdb"),
    ("FimH lectin",        "3_Klebsiella_pneumoniae/PDB_9AT9_FimH_lectin_mannose.pdb",           "polymer",                                   "MAN",             0, "FimH_receptor_clean.pdb"),
    ("ExoA",               "4_Pseudomonas_aeruginosa/PDB_1AER_ExoA.pdb",                        "chain A and polymer",                       "TAD+TIA+AMP",     0, "ExoA_receptor_clean.pdb"),
    ("ExoT-SpcS interface","4_Pseudomonas_aeruginosa/PDB_6JNP_ExoT_SpcS_complex.pdb",           "(chain A or chain B) and polymer",          "GOL",             0, "ExoT_SpcS_clean.pdb"),
]

_cfg_json = os.environ.get("OW_CONFIG_JSON", "")
if _cfg_json:
    _cfg = json.load(open(_cfg_json))
    _sd = _cfg["species"]["dir"]
    RECIPES = [
        (r["name"], "{}/{}".format(_sd, r["pdb"]), r["keep"], r["drop"], r["zn"], r["outfile"])
        for r in _cfg.get("clean", [])
    ]
else:
    RECIPES = _ALL_RECIPES

for name, pdb, keep, drop, want_zn, outfile in RECIPES:
    path = pdb
    if not os.path.exists(path):
        print("SKIP {}: {} not found".format(name, path))
        continue
    cmd.delete("all")
    cmd.load(path, "rec")
    before = cmd.count_atoms("rec")
    cmd.remove("rec and not ({})".format(keep))
    if drop:
        cmd.remove("rec and resn {}".format(drop))
    cmd.remove("rec and solvent")
    cmd.remove("rec and hydro")
    zinc = cmd.count_atoms("rec and resn ZN")
    target = os.path.join(OUT, outfile)
    cmd.save(target, "rec")
    print("{:22s} {:5d} -> {:5d} heavy atoms, Zn {}  ->  {}".format(
        name, before, cmd.count_atoms("rec"), zinc, outfile))
    if zinc != want_zn:
        print("   WARNING: expected {} Zn atom(s), kept {}".format(want_zn, zinc))

print("\nNote: BoNT/B receptor-binding keeps protein only — its target site is the")
print("doxorubicin pocket, 81 A from the catalytic Zn, so the metal is not needed there.")
print("For any catalytic-site receptor the Zn must survive; the counts above are the check.")
