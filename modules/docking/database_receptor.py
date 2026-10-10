#!/usr/bin/env python3
"""Cut one chain out of an experimental structure and make it a receptor.

    database_receptor.py <in.pdb> <chain> <crystal_ligand> <metal> \
                         <box_size> <out_stem>

      <chain>           chain to dock into, e.g. B. "" keeps every chain.
      <crystal_ligand>  residue name of the ligand that crystallised in the
                        site, e.g. MAN. "" for none.
      <metal>           metal to keep in the receptor, e.g. ZN. "" for none.
      <box_size>        box edge in angstroms, only used to check the fit.
      <out_stem>        path and name stem for the three files below.

Writes, for the calling script:

    <out_stem>.clean.pdb            the receptor
    <out_stem>.centre.tsv           x  y  z  what the centre was taken from
    <out_stem>.crystal_ligand.pdb   the crystal ligand, if one was named

The centre file has the same four columns prep_receptor.py writes for the
predicted models, so receptor_report.py reads it unchanged.

Why this is a different program from prep_receptor.py. A predicted model has
no ligand and usually no metal, so its box has to come from a guess. An
experimental structure can hold the answer: the ligand that crystallised in
the site is where the site is. The box is centred on, in order:

    the crystal ligand   the strongest basis there is, because it is also the
                         answer key. Program I redocks that ligand and asks
                         whether the pose comes back.
    the metal            if the site has one and no ligand was named.
    the whole protein    a last resort, labelled whole_protein_blind so it
                         cannot be mistaken for a chosen site later.

What is removed from the receptor, and why it matters here:

    other chains   The file holds the whole asymmetric unit. Chain B of a
                   dimer is a second copy of the site, or a peptide, or a
                   chaperone, and a box centred on one chain has to be docked
                   into that chain alone.
    water          Water in the pocket is a wall the ligand cannot pass.
    other ligands  The sulfate sitting in the 1EPW catalytic site would be
                   docked around rather than displaced. A bound ligand left
                   in the site is the site, filled.
    alternate      Two conformations of one side chain at once is two side
    locations      chains in the same place. Altloc A is kept.

Missing side chains. Crystal structures leave out atoms it could not see: a
surface lysine with no electron density is modelled as a stub. Meeko has no
template for a residue with atoms missing and, told to carry on, deletes it
without saying so. If the missing residue lines the pocket, that deletion
changes the answer. So the missing atoms are rebuilt with PDBFixer when it is
installed (conda install -c conda-forge pdbfixer). When it is not, the
receptor is left as it is and receptor_report.py lists, after Meeko, every
residue that was deleted and marks those inside the box.

Notes go to stderr; stdout is one summary line.
"""
import math
import os
import sys

import numpy
from pymol import cmd

# Room, in angstroms, a ligand needs beyond its crystal length. A flexible
# ligand can be longer than the shape it crystallised in. A box that fits is
# necessary, not sufficient: a ligand with many rotatable bonds can fail
# inside a box that is large enough.
MARGIN = 4.0


def centroid(selection):
    """The geometric centre of a selection, hydrogens already removed."""
    return cmd.get_coords(selection).mean(axis=0)


def longest_span(selection):
    """The largest distance between any two atoms of the selection."""
    coords = cmd.get_coords(selection)
    gaps = coords[:, None, :] - coords[None, :, :]
    return float(numpy.sqrt((gaps ** 2).sum(axis=2)).max())


def rebuild(path):
    """Add the heavy atoms the crystal left out. Returns a note for stderr."""
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError:
        return ("NOTE   PDBFixer is not installed, so missing side-chain atoms "
                "were not rebuilt. Meeko will delete any incomplete residue; "
                "read the list receptor_report.py prints.")

    fixer = PDBFixer(filename=path)
    fixer.findMissingResidues()
    fixer.missingResidues = {}      # a gap in the chain is not filled de novo
    fixer.findMissingAtoms()
    added = sum(len(atoms) for atoms in fixer.missingAtoms.values())
    fixer.addMissingAtoms()
    with open(path, "w", encoding="utf-8") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    # PDBFixer adds a chemically correct C-terminal OXT, which Meeko's residue
    # templates reject. It is at the chain end, nowhere near a pocket.
    with open(path, encoding="utf-8") as handle:
        lines = [line for line in handle if line[12:16].strip() != "OXT"]
    with open(path, "w", encoding="utf-8") as handle:
        handle.writelines(lines)
    return f"rebuilt {added} missing heavy atom(s) with PDBFixer"


def main(argv):
    if len(argv) != 6:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    source, chain, ligand, metal, box_size, stem = argv
    ligand, metal, box_size = ligand.upper(), metal.upper(), float(box_size)

    cmd.delete("all")
    cmd.load(source, "src")

    # Altloc A where a choice was recorded, then clear the label itself, as
    # prep_receptor.py does: the survivors would still announce an alternative
    # that no longer exists.
    cmd.remove("src and hydro")
    cmd.remove("src and not alt ''+A")
    cmd.alter("src", "alt=''")
    cmd.sort("src")

    here = f"src and chain {chain}" if chain else "src"
    if cmd.count_atoms(f"{here} and polymer") == 0:
        print(f"no polymer in chain '{chain}' of {source}", file=sys.stderr)
        return 1

    keep = "polymer"
    kept_metal = False
    if metal:
        kept_metal = cmd.count_atoms(f"{here} and resn {metal}") > 0
        if kept_metal:
            keep += f" or resn {metal}"
        else:
            print(f"NOTE   no {metal} in chain '{chain}' of {source}",
                  file=sys.stderr)
    cmd.create("rec", f"{here} and ({keep})")

    removed = sorted({atom.resn for atom in
                      cmd.get_model(f"{here} and not solvent and not ({keep})").atom})
    if removed:
        print(f"    left out of the receptor: {', '.join(removed)}", file=sys.stderr)

    out_dir = os.path.dirname(os.path.abspath(stem))
    os.makedirs(out_dir, exist_ok=True)

    if ligand:
        cmd.create("lig", f"{here} and resn {ligand}")
        atoms = cmd.get_model("lig").atom
        if not atoms:
            print(f"no residue {ligand} in chain '{chain}' of {source}",
                  file=sys.stderr)
            return 1
        # Several copies in one chain: take the first, so the answer key is
        # one molecule and not two on top of each other.
        first = atoms[0].resi
        cmd.remove(f"lig and not resi {first}")
        cmd.save(f"{stem}.crystal_ligand.pdb", "lig")

        centre, basis = centroid("lig"), f"ligand:{ligand}"
        span = longest_span("lig")
        needed = math.ceil(span + MARGIN)
        if box_size < needed:
            print(f"WARNING {ligand} spans {span:.1f} A and the box is "
                  f"{box_size:.0f} A — too tight for it to move. Raise "
                  f"BOX_SIZE to at least {needed}.", file=sys.stderr)
    elif kept_metal:
        centre, basis = centroid(f"rec and resn {metal}"), f"metal:{metal}"
    else:
        centre, basis = centroid("rec"), "whole_protein_blind"

    cmd.save(f"{stem}.clean.pdb", "rec")
    print(f"    {rebuild(stem + '.clean.pdb')}", file=sys.stderr)

    with open(f"{stem}.centre.tsv", "w", encoding="utf-8") as handle:
        handle.write(f"{centre[0]:.3f}\t{centre[1]:.3f}\t{centre[2]:.3f}"
                     f"\t{basis}\n")

    name = os.path.basename(stem)
    print(f"--- {name}: centre {centre[0]:.3f} {centre[1]:.3f} "
          f"{centre[2]:.3f}   basis {basis}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
