#!/usr/bin/env python3
"""Step 8c — Mean pLDDT per model and per region of interest.

pLDDT lives in the B-factor column of an AlphaFold file. Judge it over the
site you intend to dock into, not over the whole chain.

Env vars:
    OW_DATASETS    Path to WORKING_FOLDER/INPUT_DATASETS/from_Database
    OW_CONFIG_JSON (optional) JSON config for single-species mode
"""
import os
import statistics
import json

DATASETS = os.environ["OW_DATASETS"]

_ALL_REGIONS = [
    # path (relative to WORKING_FOLDER/INPUT_DATASETS/from_Database/),     label,                      first, last
    ("1_Mycobacterium_tuberculosis/MODEL_O05442_CpnT_full_AF.pdb",         "whole chain",              None, None),
    ("1_Mycobacterium_tuberculosis/MODEL_O05442_CpnT_full_AF.pdb",         "TNT domain 651-846",        651,  846),
    ("1_Mycobacterium_tuberculosis/MODEL_O05442_CpnT_full_AF.pdb",         "linker 397-634",            397,  634),
    ("4_Pseudomonas_aeruginosa/MODEL_Q9I788_ExoT_AF.pdb",                  "whole chain",              None, None),
    ("4_Pseudomonas_aeruginosa/MODEL_Q9I788_ExoT_AF.pdb",                  "ADPRT domain 236-457",      236,  457),
    ("4_Pseudomonas_aeruginosa/MODEL_Q9I788_ExoT_AF.pdb",                  "GAP domain 78-235",          78,  235),
    ("4_Pseudomonas_aeruginosa/MODEL_Q9I788_ExoT_AF.pdb",                  "chaperone-binding 23-79",    23,   79),
    ("4_Pseudomonas_aeruginosa/MODEL_G3XDA1_ExoS_AF.pdb",                  "whole chain",              None, None),
    ("4_Pseudomonas_aeruginosa/MODEL_G3XDA1_ExoS_AF.pdb",                  "ADPRT domain 232-453",      232,  453),
    ("3_Klebsiella_pneumoniae/MODEL_P21648_MrkD_AF.pdb",                    "whole chain",              None, None),
    ("3_Klebsiella_pneumoniae/MODEL_P21648_MrkD_AF.pdb",                    "minus signal region 8-321",   8,  321),
    ("3_Klebsiella_pneumoniae/MODEL_A0A0H3H2I8_FimH_full_AF.pdb",          "whole chain",              None, None),
    ("3_Klebsiella_pneumoniae/MODEL_A0A0H3H2I8_FimH_full_AF.pdb",          "lectin domain 1-160",         1,  160),
    ("2_Clostridium_botulinum/MODEL_Q00496_BoNT_E_AF.pdb",                 "whole chain",              None, None),
    ("2_Clostridium_botulinum/MODEL_P30996_BoNT_F_AF.pdb",                 "whole chain",              None, None),
    ("4_Pseudomonas_aeruginosa/MODEL_O34208_ExoU_AF.pdb",                  "whole chain",              None, None),
    ("4_Pseudomonas_aeruginosa/MODEL_O34208_ExoU_AF.pdb",                  "PLA2 region 107-357",       107,  357),
]

_cfg_json = os.environ.get("OW_CONFIG_JSON", "")
if _cfg_json:
    _cfg = json.load(open(_cfg_json))
    _sd = _cfg["species"]["dir"]
    REGIONS = [
        ("{}/{}".format(_sd, r["model"]), r["label"], r.get("first"), r.get("last"))
        for r in _cfg.get("plddt", [])
    ]
else:
    REGIONS = _ALL_REGIONS


def per_residue_plddt(path):
    scores = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("ATOM"):
                scores[int(line[22:26])] = float(line[60:66])
    return scores


def verdict(mean):
    if mean >= 90: return "very high — backbone and side chains reliable"
    if mean >= 70: return "confident backbone"
    if mean >= 50: return "low — treat with care"
    return "no structure predicted here"

print("{:55s} {:26s} {:>6s} {:>6s}  reading".format("model", "region", "pLDDT", "res"))
print("-" * 116)
missing = []
for relpath, label, first, last in REGIONS:
    path = relpath
    if not os.path.exists(path):
        if relpath not in missing:
            missing.append(relpath)
        continue
    scores = per_residue_plddt(path)
    if first is not None:
        scores = {r: b for r, b in scores.items() if first <= r <= last}
    if not scores:
        print("{:55s} {:26s} {:>6s}".format(os.path.basename(relpath), label, "n/a"))
        continue
    mean = statistics.mean(scores.values())
    print("{:55s} {:26s} {:6.1f} {:6d}  {}".format(
        os.path.basename(relpath), label, mean, len(scores), verdict(mean)))

for relpath in missing:
    print("SKIP (not downloaded yet):", relpath)
print("\nA low whole-chain mean does not condemn a model: CpnT averages ~65 because of a")
print("disordered linker you will delete anyway, while its TNT domain scores ~90.")
