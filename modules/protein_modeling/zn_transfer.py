#!/usr/bin/env python3
"""Step 10 — Transfer the catalytic Zn into models that lack it.

Env vars:
    OW_CONFIG_JSON (optional) JSON config for single-species mode
"""
import os
import json
from pymol import cmd

_cfg_json = os.environ.get("OW_CONFIG_JSON", "")
if _cfg_json:
    _cfg = json.load(open(_cfg_json))
    _sd = _cfg["species"]["dir"]
    _zn_pairs = _cfg.get("zn_transfer", [])
    _zn_tmpl = _cfg.get("zn_template", {})
    if not _zn_pairs or not _zn_tmpl:
        print("no Zn transfer defined for this dataset — skipping")
        raise SystemExit(0)
    TEMPLATE = "{}/{}".format(_sd, _zn_tmpl["path"])
    BONT_MODELS = _sd
    PAIRS = [(p["label"], p["model"], p["outfile"]) for p in _zn_pairs]
else:
    TEMPLATE = "2_Clostridium_botulinum/PDB_1XTG_BoNT_A_LC_SNAP25.pdb"
    BONT_MODELS = "2_Clostridium_botulinum"
    PAIRS = [
        ("BoNT/E", "MODEL_Q00496_BoNT_E_AF.pdb", "MODEL_BoNT_E_with_Zn.pdb"),
        ("BoNT/F", "MODEL_P30996_BoNT_F_AF.pdb", "MODEL_BoNT_F_with_Zn.pdb"),
    ]

if not os.path.exists(TEMPLATE):
    raise SystemExit("template {} not found".format(TEMPLATE))

for label, filename, outfile in PAIRS:
    path = os.path.join(BONT_MODELS, filename)
    if not os.path.exists(path):
        print("SKIP {}: {} not downloaded".format(label, filename))
        continue
    cmd.delete("all")
    cmd.load(path, "pred")
    cmd.load(TEMPLATE, "template")
    rms, atoms = cmd.super("template and polymer", "pred")[:2]
    cmd.create("with_zn", "pred or (template and resn ZN)")
    out = os.path.join(BONT_MODELS, outfile)
    cmd.save(out, "with_zn")
    zinc = cmd.count_atoms("with_zn and resn ZN")
    print("{}: template superposed at {:.2f} A over {} atoms, Zn atoms now {}  ->  {}/{}".format(
        label, rms, atoms, zinc, BONT_MODELS, outfile))
    cmd.select("zn_shell", "byres (with_zn and polymer and not hydro within 2.6 of "
                           "(with_zn and resn ZN))")
    shell = []
    cmd.iterate("zn_shell and name CA", "shell.append((resn, resi))", space={"shell": shell})
    print("   Zn coordination shell: " + (", ".join(r + i for r, i in shell) or "none within 2.6 A"))
    if len(shell) < 3:
        print("   WARNING: fewer than three coordinating residues — check the placement")
        print("   in PyMOL before docking into this site.")
