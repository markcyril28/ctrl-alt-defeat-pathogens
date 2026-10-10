#!/usr/bin/env python3
"""Recover a crystal ligand's chemistry, and keep its pose as the answer key.

    database_ligand.py <crystal_ligand.pdb> <code> <out.crystal.sdf> <ccd_dir>

      <crystal_ligand.pdb>  the ligand as it sat in the structure, written by
                            database_receptor.py
      <code>                the three-character residue name, e.g. MAN, TAD
      <out.crystal.sdf>     where to write the ligand's crystal pose, with bond
                            orders, for redock_rmsd.py to compare against
      <ccd_dir>             folder that keeps the downloaded chemistry, so a
                            rerun does not download it again

Prints the ligand's SMILES, alone on one line, for the calling script to hand
to smiles_to_sdf.py. Everything else goes to stderr.

A PDB file records where each atom is and almost nothing about the bonds:
single, double and aromatic are not in it, and a ligand read from one has its
bonds guessed from distances. Docked as it stands, a guessed ligand is not the
molecule that crystallised. The RCSB keeps the true chemistry of every ligand
in its Chemical Component Dictionary, one file per three-character code. This
downloads that file, and fits the crystal atoms onto it, which gives the
crystal pose with correct bond orders.

The SMILES that comes out of this is the CCD's, not the crystal's. That is
deliberate. Program I docks a fresh conformer built from the SMILES, not the
crystal geometry itself, so the redocking test starts from a shape that has
never seen the answer. Handing Vina the crystal conformer would let it keep
the ring pucker and the torsions it was supposed to find.

What the CCD does not say is the protonation state the ligand had in the
crystallisation buffer. Phosphates and amines come back as the dictionary
draws them. See the top of smiles_to_sdf.py.
"""
import os
import sys
import time
import urllib.request

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")      # RDKit's warnings go to stderr unasked

CCD_URL = "https://files.rcsb.org/ligands/download/{code}_ideal.sdf"
ATTEMPTS = 5                        # a flaky connection is the usual failure


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    crystal_pdb, code, destination, ccd_dir = argv
    code = code.upper()

    os.makedirs(ccd_dir, exist_ok=True)
    template_sdf = os.path.join(ccd_dir, f"{code}_ideal.sdf")
    if not os.path.exists(template_sdf):
        url = CCD_URL.format(code=code)
        problem = None
        for attempt in range(ATTEMPTS):
            try:
                with urllib.request.urlopen(url, timeout=30) as reply:
                    body = reply.read()
                # Written under another name and renamed, so a download that
                # dies halfway never leaves a truncated file that the next
                # run would read as the real chemistry.
                with open(template_sdf + ".part", "wb") as handle:
                    handle.write(body)
                os.replace(template_sdf + ".part", template_sdf)
                break
            except OSError as failure:
                problem = failure
                time.sleep(5)
        else:
            print(f"FAIL {code} — could not download {url}: {problem}\n"
                  f"      Offline? Save that file as {template_sdf} and run again.",
                  file=sys.stderr)
            return 1

    template = Chem.MolFromMolFile(template_sdf)
    if template is None:
        print(f"FAIL {code} — RDKit cannot read {template_sdf}", file=sys.stderr)
        return 1
    template = Chem.RemoveHs(template)

    crystal = Chem.MolFromPDBFile(crystal_pdb, removeHs=True, sanitize=False)
    if crystal is None:
        print(f"FAIL {code} — RDKit cannot read {crystal_pdb}", file=sys.stderr)
        return 1

    try:
        posed = AllChem.AssignBondOrdersFromTemplate(template, crystal)
    except ValueError as problem:
        print(f"FAIL {code} — the crystal atoms do not fit the dictionary entry "
              f"({crystal.GetNumAtoms()} atoms in the structure, "
              f"{template.GetNumAtoms()} in the CCD): {problem}\n"
              f"      A ligand with atoms missing from the density cannot be "
              f"used as an answer key.", file=sys.stderr)
        return 1

    Chem.MolToMolFile(posed, destination)

    smiles = Chem.MolToSmiles(template)
    print(f"    OK {code} crystal pose with bond orders -> {destination}",
          file=sys.stderr)
    print(smiles)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
