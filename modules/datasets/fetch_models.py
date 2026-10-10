#!/usr/bin/env python
"""fetch_models.py — download AlphaFold DB models for the toxin targets.

Lives in modules/datasets/. Run it from WORKING_FOLDER/INPUT_DATASETS/from_Database/, because the species folders are
relative to it. Needs only Python 3 and internet access:

    python ../../../modules/datasets/fetch_models.py             # download every available model
    python ../../../modules/datasets/fetch_models.py --check     # report availability, download nothing
    python ../../../modules/datasets/fetch_models.py P21648      # one accession only

Models land in <species>/ as MODEL_<ACCESSION>_<Name>_AF.pdb. Experimental
structures (PDB_*.pdb) in the same folder are never touched.

Predicted models are NOT experimental data. The B-factor column of an AlphaFold
file holds the per-residue pLDDT confidence (0-100), not a crystallographic
B-factor. Check it before docking — see Step 8 of the PyMOL workshop manual.
"""

import argparse
import json
import os
import sys
import urllib.error
import urllib.request

API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

# (accession, short_name, species_dir, note)
TARGETS = [
    ("P21648", "MrkD", "3_Klebsiella_pneumoniae", "K. pneumoniae type 3 fimbrial adhesin — no PDB in dataset"),
    ("G3XDA1", "ExoS", "4_Pseudomonas_aeruginosa", "P. aeruginosa T3SS effector — no PDB in dataset"),
    ("Q9I788", "ExoT", "4_Pseudomonas_aeruginosa", "ExoT catalytic domains — 6JNP/4JMF cover residues 23-79 only"),
    ("O34208", "ExoU", "4_Pseudomonas_aeruginosa", "P. aeruginosa phospholipase A2 — no PDB in dataset"),
    ("P18640", "BoNT_C1", "2_Clostridium_botulinum", "BoNT serotype C1 — sequence only"),
    ("P19321", "BoNT_D", "2_Clostridium_botulinum", "BoNT serotype D — sequence only"),
    ("Q00496", "BoNT_E", "2_Clostridium_botulinum", "BoNT serotype E — sequence only"),
    ("P30996", "BoNT_F", "2_Clostridium_botulinum", "BoNT serotype F — sequence only"),
    ("Q60393", "BoNT_G", "2_Clostridium_botulinum", "BoNT serotype G — sequence only"),
    ("O05442", "CpnT_full", "1_Mycobacterium_tuberculosis", "full-length CpnT — 4QLP covers the TNT domain only"),
    ("A0A0H3H2I8", "FimH_full", "3_Klebsiella_pneumoniae", "full-length FimH — 9AT9 covers the lectin domain only"),
]


def query(accession):
    """Return the AlphaFold DB record for an accession, or None if there is none."""
    try:
        with urllib.request.urlopen(API.format(acc=accession), timeout=30) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None
        raise
    if not isinstance(payload, list) or not payload:
        return None
    return payload[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accessions", nargs="*", help="limit to these accessions")
    parser.add_argument("--check", action="store_true", help="report availability only")
    args = parser.parse_args()

    wanted = [t for t in TARGETS if not args.accessions or t[0] in args.accessions]
    if not wanted:
        sys.exit(f"No known target matches {args.accessions}. "
                 f"Known: {', '.join(t[0] for t in TARGETS)}")

    missing = []
    for accession, name, species_dir, note in wanted:
        record = query(accession)
        if record is None:
            print(f"[no model] {accession:12s} {name:12s} {note}")
            missing.append((accession, name))
            continue
        url = record["pdbUrl"]
        version = record.get("latestVersion", "?")
        residues = record.get("uniprotEnd", "?")
        if args.check:
            print(f"[available] {accession:12s} {name:12s} v{version}, {residues} residues")
            continue
        out_dir = species_dir
        os.makedirs(out_dir, exist_ok=True)
        destination = os.path.join(out_dir, f"MODEL_{accession}_{name}_AF.pdb")
        urllib.request.urlretrieve(url, destination)
        size = os.path.getsize(destination)
        print(f"[saved] {destination}  (v{version}, {residues} residues, {size // 1024} kB)")

    if missing:
        print("\nAlphaFold DB has no model for:")
        for accession, name in missing:
            print(f"  {name} ({accession})")
        print("Build these by template-based modelling instead — the dataset already "
              "contains\nclose homologues to use as templates (BoNT/B 1EPW or BoNT/A "
              "3BTA for any BoNT\nserotype). See Step 8d of the PyMOL workshop manual.")


if __name__ == "__main__":
    main()
