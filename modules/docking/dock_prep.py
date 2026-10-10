#!/usr/bin/env python
"""dock_prep.py — docking preparation helper for the Ologist Workshop toxin dataset.

Lives in modules/docking/. Run it from WORKING_FOLDER/INPUT_DATASETS/from_Database/ (paths below are relative to it)
inside the `protein_modeling` conda environment
(made by setup_protein_modeling_conda_envs.sh: pymol-open-source, vina, meeko, pdbfixer, rdkit).

Subcommands
-----------
  split     Split a PDB into a receptor and a reference ligand; report the
            ligand centre and any residues with missing side-chain atoms.
  receptor  Rebuild missing side chains (pdbfixer) and write receptor PDBQT
            plus a Vina box file (meeko).
  ligand    Build a ligand PDBQT with correct bond orders, taken from the RCSB
            chemical component dictionary rather than perceived from geometry.
  rmsd      Symmetry-aware, in-place RMSD of docked poses against a reference.

Typical run (FimH / mannose redocking control):

  python ../../../modules/docking/dock_prep.py split    3_Klebsiella_pneumoniae/PDB_9AT9_FimH_lectin_mannose.pdb MAN --out fimh
  python ../../../modules/docking/dock_prep.py receptor fimh/receptor.pdb --ref-ligand fimh/ligand_ref.pdb --out fimh/fimh
  python ../../../modules/docking/dock_prep.py ligand   fimh/ligand_ref.pdb MAN --out fimh/ligand.pdbqt
  vina --receptor fimh/fimh.pdbqt --ligand fimh/ligand.pdbqt \\
       --config fimh/fimh.box.txt --exhaustiveness 16 --seed 42 --out fimh/poses.pdbqt
  mk_export.py fimh/poses.pdbqt -s fimh/poses.sdf
  python ../../../modules/docking/dock_prep.py rmsd fimh/poses.sdf fimh/ligand_correct.sdf
"""

import argparse
import math
import os
import sys
import urllib.request

# Heavy-atom count of each complete amino acid residue. A residue with fewer
# atoms than this has a disordered side chain that was not modelled in the
# crystal, which meeko's residue templates reject.
HEAVY_ATOMS = {
    "ALA": 5, "ARG": 11, "ASN": 8, "ASP": 8, "CYS": 6, "GLN": 9, "GLU": 9,
    "GLY": 4, "HIS": 10, "ILE": 8, "LEU": 8, "LYS": 9, "MET": 8, "PHE": 11,
    "PRO": 7, "SER": 6, "THR": 7, "TRP": 14, "TYR": 12, "VAL": 7,
}

CCD_URL = "https://files.rcsb.org/ligands/download/{code}_ideal.sdf"


# ---------------------------------------------------------------- utilities

def read_atoms(path):
    """Yield (record, name, resn, chain, resi, x, y, z, element) per atom line."""
    atoms = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                atoms.append((
                    line[:6].strip(), line[12:16].strip(), line[17:20].strip(),
                    line[21], int(line[22:26]),
                    float(line[30:38]), float(line[38:46]), float(line[46:54]),
                    (line[76:78].strip() or line[12:16].strip()[:1]),
                ))
    return atoms


def centroid(atoms):
    return [round(sum(a[5 + i] for a in atoms) / len(atoms), 3) for i in range(3)]


def distance(a, b):
    return math.dist(a[5:8], b[5:8])


def start_pymol():
    """Start PyMOL headless and return its cmd module."""
    import pymol
    pymol.finish_launching(["pymol", "-qc"])
    return pymol.cmd


# -------------------------------------------------------------------- split

def cmd_split(args):
    """Extract receptor and reference ligand, and audit the receptor."""
    os.makedirs(args.out, exist_ok=True)
    receptor_pdb = os.path.join(args.out, "receptor.pdb")
    ligand_pdb = os.path.join(args.out, "ligand_ref.pdb")

    keep = "polymer"
    if args.keep_metals:
        keep += " or resn ZN+MG+MN+CA+NA+FE+CU+NI+CO"

    cmd = start_pymol()
    cmd.load(args.pdb, "src")
    cmd.create("rec", f"src and ({keep}) and not resn {args.ligand}")
    cmd.create("lig", f"src and resn {args.ligand}")
    if cmd.count_atoms("lig") == 0:
        sys.exit(f"ERROR: no residue named {args.ligand} in {args.pdb}")
    cmd.remove("rec and hydro")
    cmd.save(receptor_pdb, "rec")
    cmd.save(ligand_pdb, "lig")
    print(f"receptor : {receptor_pdb}  ({cmd.count_atoms('rec')} heavy atoms)")
    print(f"ligand   : {ligand_pdb}  ({cmd.count_atoms('lig')} atoms)")

    ligand = [a for a in read_atoms(ligand_pdb) if a[8] != "H"]
    print(f"ligand centre (box centre): {centroid(ligand)}")

    # Which residues are incomplete, and do any of them line the pocket?
    receptor = read_atoms(receptor_pdb)
    by_residue = {}
    for atom in receptor:
        by_residue.setdefault((atom[3], atom[4]), []).append(atom)

    pocket = {key for key, ats in by_residue.items()
              if any(distance(a, l) < args.pocket_radius for a in ats for l in ligand)}
    print(f"\npocket residues within {args.pocket_radius} A of the ligand: "
          + ", ".join(f"{by_residue[k][0][2]}{k[1]}" for k in sorted(pocket, key=lambda k: k[1])))

    incomplete = []
    for key, ats in sorted(by_residue.items(), key=lambda kv: kv[0][1]):
        resn = ats[0][2]
        expected = HEAVY_ATOMS.get(resn)
        if expected and len(ats) < expected:
            nearest = min(distance(a, l) for a in ats for l in ligand)
            incomplete.append((key, resn, len(ats), expected, nearest, key in pocket))

    if not incomplete:
        print("\nAll residues are complete — no side-chain rebuilding needed.")
        return
    print(f"\n{len(incomplete)} residue(s) with missing side-chain atoms:")
    print(f"  {'residue':>12} {'atoms':>6} {'expected':>9} {'dist to ligand':>15}")
    for key, resn, have, expected, nearest, in_pocket in incomplete:
        flag = "  <-- IN POCKET: rebuild, do not delete" if in_pocket else ""
        print(f"  {resn + str(key[1]) + ' ' + key[0]:>12} {have:6d} {expected:9d} {nearest:12.2f} A{flag}")
    if any(row[5] for row in incomplete):
        print("\nAt least one incomplete residue lines the binding site. Rebuild it "
              "(the `receptor` subcommand does this) — deleting it would remove part\n"
              "of the pocket and invalidate the docking result.")
    else:
        print("\nNone of the incomplete residues line the binding site.")


# ----------------------------------------------------------------- receptor

def cmd_receptor(args):
    """pdbfixer rebuild -> strip OXT -> meeko receptor PDBQT + Vina box."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    stem = args.out
    os.makedirs(os.path.dirname(stem) or ".", exist_ok=True)
    fixed_pdb = stem + "_fixed.pdb"
    ready_pdb = stem + "_ready.pdb"

    print(f"rebuilding missing heavy atoms in {args.receptor} ...")
    fixer = PDBFixer(filename=args.receptor)
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # do not build unresolved loops de novo
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    with open(fixed_pdb, "w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    # pdbfixer adds a chemically correct C-terminal OXT, but meeko's residue
    # templates reject it. Dropping it does not affect the binding site.
    with open(fixed_pdb) as src, open(ready_pdb, "w") as dst:
        for line in src:
            if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "OXT":
                continue
            if line.startswith(("ATOM", "HETATM")) and line[76:78].strip() == "H":
                continue
            dst.write(line)
    before = len([a for a in read_atoms(args.receptor)])
    after = len([a for a in read_atoms(ready_pdb)])
    print(f"heavy atoms: {before} -> {after} (rebuilt side chains, OXT removed)")

    if args.ref_ligand:
        ligand = [a for a in read_atoms(args.ref_ligand) if a[8] != "H"]
        center = centroid(ligand)
    else:
        center = [args.center_x, args.center_y, args.center_z]
    print(f"box centre {center}, box size {args.box_size} A")

    import subprocess
    call = [
        "mk_prepare_receptor.py", "--read_pdb", ready_pdb, "-o", stem,
        "-p", "-v",
        "--box_center", *[str(c) for c in center],
        "--box_size", *[str(args.box_size)] * 3,
    ]
    result = subprocess.run(call, capture_output=True, text=True)
    sys.stdout.write(result.stdout[-1500:])
    if result.returncode != 0:
        sys.stderr.write(result.stderr[-2500:])
        sys.exit(
            "\nmeeko rejected the receptor. Inspect the residues it names. If they are "
            "surface\nresidues far from the pocket (check with `split`), re-run adding "
            "--allow-bad-res."
        )
    print(f"\nwrote {stem}.pdbqt, {stem}.box.txt, {stem}.box.pdb")
    print(f"visualise the box in PyMOL with:  load {stem}.box.pdb")


# ------------------------------------------------------------------- ligand

def cmd_ligand(args):
    """Build a ligand PDBQT whose bond orders come from the RCSB CCD."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import subprocess

    code = args.resname.upper()
    template_sdf = os.path.join(os.path.dirname(args.out) or ".", f"{code}_ideal.sdf")
    if not os.path.exists(template_sdf):
        url = CCD_URL.format(code=code)
        print(f"downloading reference chemistry: {url}")
        urllib.request.urlretrieve(url, template_sdf)

    template = Chem.MolFromMolFile(template_sdf)
    if template is None:
        sys.exit(f"ERROR: could not read {template_sdf}")
    template = Chem.RemoveHs(template)

    crystal = Chem.MolFromPDBFile(args.ligand_pdb, removeHs=True, sanitize=False)
    if crystal is None:
        sys.exit(f"ERROR: could not read {args.ligand_pdb}")

    molecule = AllChem.AssignBondOrdersFromTemplate(template, crystal)
    molecule = Chem.AddHs(molecule, addCoords=True)
    reference_sdf = os.path.join(os.path.dirname(args.out) or ".", "ligand_correct.sdf")
    Chem.MolToMolFile(molecule, reference_sdf)
    print(f"reference SMILES : {Chem.MolToSmiles(Chem.RemoveHs(molecule))}")
    print(f"wrote {reference_sdf}  (use this as the RMSD reference)")

    result = subprocess.run(
        ["mk_prepare_ligand.py", "-i", reference_sdf, "-o", args.out],
        capture_output=True, text=True,
    )
    sys.stdout.write(result.stdout[-800:])
    if result.returncode != 0:
        sys.stderr.write(result.stderr[-1500:])
        sys.exit("ERROR: mk_prepare_ligand.py failed")
    print(f"wrote {args.out}")


# --------------------------------------------------------------------- rmsd

def cmd_rmsd(args):
    """Symmetry-aware RMSD without superposition — the redocking test."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    reference = Chem.MolFromMolFile(args.reference, removeHs=True)
    if reference is None:
        sys.exit(f"ERROR: could not read {args.reference}")

    print(f"{'pose':>5} {'RMSD (A)':>10}   verdict")
    verdicts = []
    for index, pose in enumerate(Chem.SDMolSupplier(args.poses, removeHs=True), 1):
        if pose is None:
            print(f"{index:5d} {'unreadable':>10}")
            continue
        value = rdMolAlign.CalcRMS(pose, reference)
        verdict = "reproduces the crystal pose" if value <= 2.0 else "different pose"
        verdicts.append(value)
        print(f"{index:5d} {value:10.2f}   {verdict}")
        if index >= args.top:
            break
    if verdicts:
        best = verdicts[0]
        print()
        if best <= 2.0:
            print(f"PASS: top-ranked pose is {best:.2f} A from the crystal ligand "
                  "(<= 2.0 A).\nThe protocol reproduces a known answer, so scores from "
                  "it are worth interpreting.")
        else:
            print(f"FAIL: top-ranked pose is {best:.2f} A from the crystal ligand "
                  "(> 2.0 A).\nDo not trust scores from this setup for this site — see "
                  "Step 11 of the manual.")


# --------------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest="command", required=True)

    split = subparsers.add_parser("split", help="split into receptor + reference ligand")
    split.add_argument("pdb")
    split.add_argument("ligand", help="residue name of the ligand, e.g. MAN")
    split.add_argument("--out", default="dock", help="output folder")
    split.add_argument("--keep-metals", action="store_true",
                       help="keep catalytic metals in the receptor")
    split.add_argument("--pocket-radius", type=float, default=5.0)
    split.set_defaults(func=cmd_split)

    receptor = subparsers.add_parser("receptor", help="receptor PDBQT + Vina box")
    receptor.add_argument("receptor")
    receptor.add_argument("--out", required=True, help="output basename")
    receptor.add_argument("--ref-ligand", help="PDB of the reference ligand (box centre)")
    receptor.add_argument("--center-x", type=float, dest="center_x", default=0.0)
    receptor.add_argument("--center-y", type=float, dest="center_y", default=0.0)
    receptor.add_argument("--center-z", type=float, dest="center_z", default=0.0)
    receptor.add_argument("--box-size", type=float, default=20.0)
    receptor.set_defaults(func=cmd_receptor)

    ligand = subparsers.add_parser("ligand", help="ligand PDBQT with correct chemistry")
    ligand.add_argument("ligand_pdb")
    ligand.add_argument("resname", help="PDB chemical component code, e.g. MAN, DM2")
    ligand.add_argument("--out", required=True, help="output .pdbqt")
    ligand.set_defaults(func=cmd_ligand)

    rmsd = subparsers.add_parser("rmsd", help="RMSD of poses vs a reference")
    rmsd.add_argument("poses", help="SDF written by mk_export.py")
    rmsd.add_argument("reference", help="ligand_correct.sdf from the `ligand` step")
    rmsd.add_argument("--top", type=int, default=9)
    rmsd.set_defaults(func=cmd_rmsd)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
