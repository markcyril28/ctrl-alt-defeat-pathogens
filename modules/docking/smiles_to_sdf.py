#!/usr/bin/env python3
"""Turn a SMILES string into a 3-D SDF a docking program can use.

    smiles_to_sdf.py <name> <SMILES> <out.sdf> <seed> <for_organism>

Prints one tab-separated row, with no header, for the calling script's
ligand table — so the identity can be checked against the database the
SMILES came from:

    <name>  <for_organism>  <SMILES>  <formula>  <mw>  <inchikey>
    <rotatable bonds>

<for_organism> is which organism's receptors this ligand was chosen for, and
is passed straight through to the row. A ligand has no species of its own; it
is a compound, and this column records the pairing that was intended, which
is the thing a reader of the table needs and the thing the docking program
acts on. "shared" means every receptor, which is where a control belongs.

The same values are repeated in readable form on stderr.

A SMILES string is a 2-D sketch. Docking needs coordinates, so one conformer
is generated with ETKDG and relaxed with MMFF. Three things this does not do,
each of which changes a docking result:

    It does not pick the protonation state. SMILES written on a database page
    is usually the neutral form; at pH 7 a carboxylic acid is an anion and an
    amine is a cation, and the charge is most of the electrostatics the
    scoring function sees.

    It does not search conformers. Vina rotates the torsions you give it, but
    it does not re-pucker a ring or invert a centre. A sugar handed over in
    the wrong ring conformation stays wrong.

    It does not check that you typed the right molecule. Stereochemistry is
    carried by the @ and @@ marks, and a SMILES missing them is a different
    compound that will still dock. Compare the formula and InChIKey printed
    here against PubChem or DrugBank before you trust a score.
"""
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")      # RDKit's warnings go to stderr unasked


def main(argv):
    if len(argv) != 5:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    name, smiles, destination = argv[0], argv[1], argv[2]
    seed = int(argv[3])
    for_organism = argv[4]

    flat = Chem.MolFromSmiles(smiles)
    if flat is None:
        print(f"FAIL {name} — RDKit cannot parse this SMILES: {smiles}", file=sys.stderr)
        return 1

    # The identity is taken from the SMILES, before any 3-D work. Once a
    # conformer exists, InChI reads stereochemistry off the coordinates, so a
    # SMILES with an undefined centre comes back with a stereo layer ETKDG
    # invented — a key that can never match the database page you are being
    # told to compare it against.
    try:
        key = Chem.MolToInchiKey(flat) or "-"
    except Exception:                       # RDKit built without InChI
        key = "-"

    # An undefined stereocentre is not an error, but it does mean one
    # arbitrary enantiomer gets docked. Say so rather than letting the
    # coordinates decide quietly.
    undefined = [
        centre for centre, label in
        Chem.FindMolChiralCenters(flat, includeUnassigned=True, useLegacyImplementation=False)
        if label == "?"
    ]

    # Hydrogens must be explicit before embedding, or the geometry is wrong
    # in a way nothing downstream reports.
    molecule = Chem.AddHs(flat)

    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed            # same seed, same conformer
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        print(f"FAIL {name} — no conformer could be embedded", file=sys.stderr)
        return 1

    # MMFF covers most drug-like molecules; UFF is the fallback for the rest.
    if AllChem.MMFFHasAllMoleculeParams(molecule):
        AllChem.MMFFOptimizeMolecule(molecule)
        force_field = "MMFF"
    else:
        AllChem.UFFOptimizeMolecule(molecule)
        force_field = "UFF"

    molecule.SetProp("_Name", name)

    with Chem.SDWriter(destination) as writer:
        writer.write(molecule)

    formula = rdMolDescriptors.CalcMolFormula(molecule)
    weight = Descriptors.MolWt(molecule)

    rotatable = rdMolDescriptors.CalcNumRotatableBonds(molecule)

    print(f"{name}\t{for_organism}\t{smiles}\t{formula}\t{weight:.2f}\t{key}"
          f"\t{rotatable}")
    print(f"    OK {name} formula={formula} mw={weight:.2f} inchikey={key} "
          f"rotatable={rotatable} forcefield={force_field}", file=sys.stderr)

    if undefined:
        centres = ", ".join(str(atom) for atom in undefined)
        print(f"NOTE {name} has undefined stereochemistry at atom(s) {centres} — "
              f"one enantiomer was embedded arbitrarily. Write the stereochemistry "
              f"into the SMILES if the molecule has it.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
