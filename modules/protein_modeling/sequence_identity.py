#!/usr/bin/env python3
"""Step 8e — Template identity via pairwise global alignment.

Env vars:
    OW_CONFIG_JSON (optional) JSON config for single-species mode
"""
import json
import os
from pathlib import Path

try:
    from Bio import Align
except ImportError:
    raise SystemExit("biopython is missing. conda install -c conda-forge biopython")


def sequence(path):
    return "".join(line.strip() for line in Path(path).read_text().splitlines()
                   if not line.startswith(">"))

_cfg_json = os.environ.get("OW_CONFIG_JSON", "")
if _cfg_json:
    _cfg = json.load(open(_cfg_json))
    _sd = _cfg["species"]["dir"]
    _q = _cfg.get("identity_query", [])
    _t = _cfg.get("identity_template", [])
    if not _q or not _t:
        print("no identity comparisons defined for this dataset — skipping")
        raise SystemExit(0)
    queries   = [(q["name"], "{}/{}".format(_sd, q["fasta"])) for q in _q]
    templates = [(t["name"], "{}/{}".format(_sd, t["fasta"])) for t in _t]
else:
    queries = [
        ("BoNT/C1", "2_Clostridium_botulinum/FASTA_P18640_BoNT_C1.fasta"),
        ("BoNT/D",  "2_Clostridium_botulinum/FASTA_P19321_BoNT_D.fasta"),
        ("BoNT/G",  "2_Clostridium_botulinum/FASTA_Q60393_BoNT_G.fasta"),
    ]
    templates = [
        ("BoNT/B (1EPW)", "2_Clostridium_botulinum/FASTA_P10844_BoNT_B.fasta"),
        ("BoNT/A1 (3BTA)", "2_Clostridium_botulinum/FASTA_P0DPI1_BoNT_A1.fasta"),
    ]

aligner = Align.PairwiseAligner(scoring="blastp", mode="global")

for qname, qpath in queries:
    if not Path(qpath).exists():
        print("SKIP {}: {} not found".format(qname, qpath))
        continue
    target = sequence(qpath)
    for tname, tpath in templates:
        alignment = aligner.align(target, sequence(tpath))[0]
        top, bottom = alignment[0], alignment[1]
        pairs = [(x, y) for x, y in zip(top, bottom) if x != "-" and y != "-"]
        identity = 100 * sum(x == y for x, y in pairs) / len(pairs)
        band = ("good enough to locate a binding site" if identity >= 50 else
                "locates domains; check every conclusion" if identity >= 30 else
                "fold hypothesis only")
        print("{:8s} vs {:16s} {:3.0f}% over {} aligned positions — {}".format(
            qname, tname, identity, len(pairs), band))

print("\nThese three serotypes have no AlphaFold DB entry. Build them at")
print("https://swissmodel.expasy.org/interactive using the templates above, and save")
print("the result as <species>/<ACCESSION>_<name>_SWISSMODEL.pdb with its QMEANDisCo score.")
print("Always align before estimating identity: comparing these sequences position by")
print("position without alignment reports ~9%, which is chance level and simply wrong.")
