#!/usr/bin/env python3
"""Clean a receptor for docking and print where the box should go.

    prep_receptor.py <in.pdb> <out.pdb> <centre.tsv> <metal> <site> \
                     <keep_ligands>

      <metal>         keep this metal and centre the box on it, e.g. ZN.
                      "" for none.
      <site>          centre the box on a PyMOL selection instead, e.g.
                      "resi 231+235+292". "" for none.
      <keep_ligands>  "true" retains non-metal heteroatoms, "false" removes them

Writes one tab-separated line to <centre.tsv>, with no header, so the calling
script can read a field out of it with cut:

    <x>  <y>  <z>  <what the centre was taken from>

What is removed, and why each one matters for docking:

    solvent       Water in the pocket is a wall the ligand cannot pass.
    hydrogens     Meeko adds its own; the ones in the file are in the way.
    alternate     Two conformations of one side chain docked at once is two
    locations     side chains in the same place.
    other ligands A bound ligand left in the site is the site, filled.

The metal is the exception, and the reason this script exists. A zinc
metalloprotease without its zinc is not the enzyme you meant to dock into:
the charge that organises the whole active site is simply gone, and the run
will still produce scores. Predicted models usually arrive with no metal at
all, so if the enzyme needs one you must put it back and be able to say where
it came from — the template, the homologue, or an assumption.

Choosing the centre, in the order this script tries:

    the metal       the strongest basis available, if the site has one
    --site          residues you named, from the literature or a homologue
    the whole       a last resort. A box around an entire protein is blind
    protein         docking: the search space is mostly solvent, the
                    exhaustiveness no longer means what it meant, and the top
                    pose is as likely to be a surface dent as the real site.
"""
import os
import sys

from pymol import cmd

SOLVENT = "solvent or resn HOH+WAT+DOD"


def main(argv):
    if len(argv) != 6:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    source, destination, centre_file = argv[0], argv[1], argv[2]
    metal = argv[3].upper()
    site = argv[4]
    keep_ligands = argv[5] == "true"

    cmd.delete("all")
    cmd.load(source, "rec")

    cmd.remove(f"rec and ({SOLVENT})")
    cmd.remove("rec and hydro")
    # Keep altloc A where a choice was recorded, drop the rest — then clear
    # the label itself. Removing the other conformers is not enough: the
    # survivors still carry "A" in the altLoc column, and a receptor that
    # still announces alternate locations sends Meeko down its altloc path
    # for atoms that no longer have an alternative.
    cmd.remove("rec and not alt ''+A")
    cmd.alter("rec", "alt=''")
    cmd.sort("rec")

    kept_metal = metal and cmd.count_atoms(f"rec and resn {metal}") > 0
    if not keep_ligands:
        exclude = "polymer"
        if kept_metal:
            exclude += f" or resn {metal}"
        cmd.remove(f"rec and not ({exclude})")

    if cmd.count_atoms("rec and polymer") == 0:
        print(f"no polymer left in {source}", file=sys.stderr)
        return 1

    # The centre, from the best basis available.
    if kept_metal:
        selection, basis = f"rec and resn {metal}", f"metal:{metal}"
    elif site:
        selection, basis = f"rec and ({site})", f"site:{site}"
        if cmd.count_atoms(selection) == 0:
            print(f"--site matched no atoms: {site}", file=sys.stderr)
            return 1
    else:
        selection, basis = "rec and polymer", "whole_protein_blind"

    centre = cmd.centerofmass(selection)

    os.makedirs(os.path.dirname(os.path.abspath(destination)), exist_ok=True)
    cmd.save(destination, "rec")

    with open(centre_file, "w", encoding="utf-8") as handle:
        handle.write(f"{centre[0]:.3f}\t{centre[1]:.3f}\t{centre[2]:.3f}"
                     f"\t{basis}\n")

    name = os.path.basename(destination)
    print(f"--- {name[:-10] if name.endswith('.clean.pdb') else name}: "
          f"centre {centre[0]:.3f} {centre[1]:.3f} {centre[2]:.3f}"
          f"   basis {basis}")
    if metal and not kept_metal:
        print(f"NOTE   no {metal} in this model — the box is not centred on a metal",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
