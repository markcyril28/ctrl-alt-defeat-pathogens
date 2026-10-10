#!/usr/bin/env python3
"""Step 8f — Trim low-confidence regions off the models.

Env vars:
    OW_CONFIG_JSON (optional) JSON config for single-species mode
"""
import os
import json
from pymol import cmd

_ALL_TRIMS = [
    ("1_Mycobacterium_tuberculosis/MODEL_O05442_CpnT_full_AF.pdb",   "resi 651-846", "CpnT TNT domain",  "1_Mycobacterium_tuberculosis", "MODEL_CpnT_TNT_domain_trimmed.pdb"),
    ("3_Klebsiella_pneumoniae/MODEL_P21648_MrkD_AF.pdb",             "not resi 1-7", "MrkD",             "3_Klebsiella_pneumoniae",      "MODEL_MrkD_trimmed.pdb"),
    ("4_Pseudomonas_aeruginosa/MODEL_G3XDA1_ExoS_AF.pdb",           "not resi 1-96","ExoS GAP+ADPRT",   "4_Pseudomonas_aeruginosa",     "MODEL_ExoS_trimmed.pdb"),
    ("4_Pseudomonas_aeruginosa/MODEL_Q9I788_ExoT_AF.pdb",           "resi 236-457", "ExoT ADPRT",       "4_Pseudomonas_aeruginosa",     "MODEL_ExoT_ADPRT_trimmed.pdb"),
    ("4_Pseudomonas_aeruginosa/MODEL_O34208_ExoU_AF.pdb",           "resi 107-357", "ExoU PLA2 region", "4_Pseudomonas_aeruginosa",     "MODEL_ExoU_PLA2_trimmed.pdb"),
]

_cfg_json = os.environ.get("OW_CONFIG_JSON", "")
if _cfg_json:
    _cfg = json.load(open(_cfg_json))
    _sd = _cfg["species"]["dir"]
    TRIMS = [
        ("{}/{}".format(_sd, t["source"]), t["keep"], t["label"], _sd, t["outfile"])
        for t in _cfg.get("trim", [])
    ]
else:
    TRIMS = _ALL_TRIMS

for srcpath, keep, label, outdir, outfile in TRIMS:
    if not os.path.exists(srcpath):
        print("SKIP {}: {} not downloaded".format(label, srcpath))
        continue
    cmd.delete("all")
    cmd.load(srcpath, "pred")
    before = cmd.count_atoms("pred and name CA")
    cmd.create("trimmed", "pred and ({})".format(keep))
    cmd.remove("trimmed and b < 50")
    kept = cmd.count_atoms("trimmed and name CA")
    first = cmd.get_model("trimmed and name CA").atom[0].resi if kept else "-"
    last = cmd.get_model("trimmed and name CA").atom[-1].resi if kept else "-"
    cmd.save(os.path.join(outdir, outfile), "trimmed")
    print("{:18s} {:4d} -> {:4d} residues (now {}-{})  ->  {}/{}".format(
        label, before, kept, first, last, outdir, outfile))

print("\nWrite down what you removed. A pose reported against residue numbering you")
print("silently altered is not reproducible.")
print("ExoU is the weakest model in the set (PLA2 region ~72) — interpret it cautiously.")
