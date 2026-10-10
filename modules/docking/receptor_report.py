#!/usr/bin/env python3
"""Record a prepared receptor: its box, and the residues Meeko deleted.

    receptor_report.py <out_dir> <name> <box_size> <organism>

Reads <out_dir>/<name>.centre.tsv, .clean.pdb and .pdbqt, then:

    writes <name>.dropped.txt   the deleted residues, if there were any
    prints one TSV row          for the calling script's receptors.tsv

<out_dir> is the organism's folder, and <organism> is its name repeated into
the table: the receptors are prepared one organism at a time, so the row says
which one without anyone having to read the path in the last column. The
calling script passes the name — this helper prepares one receptor and has no
view of the pipeline that grouped them.

The box itself is not written here: Meeko already writes <name>.box.txt in
the config format Vina reads, and program I points at that file rather than
being told the centre again. One box, written once, by the tool that built
it.

Meeko types a receptor by matching each residue against a template, and it
has no template for a residue with atoms missing. Told to carry on anyway, it
deletes those residues — and does not report which. So the two files are
compared: a residue present in the cleaned structure and absent from the
PDBQT was deleted.

Read that list when it appears. A residue deleted from the surface changes
the file. A residue deleted from the pocket changes the docking result, and
the score comes back looking exactly as plausible as it would have otherwise.
"""
import os
import sys

HEADER = ["receptor", "species", "centre_x", "centre_y", "centre_z",
          "box_size", "basis", "dropped_residues", "pdbqt"]


def residues_of(path, centre=None, half=0.0):
    """Every residue in a PDB-format file, and which ones lie in the box.

    A residue is identified by the columns PDB reserves for it: residue name,
    chain and sequence number, which is what makes "HIS A 235" unique.
    """
    present, in_box = set(), set()

    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            # The name is compared without its padding: PyMOL right-aligns a
            # metal's name (" ZN") and Meeko left-aligns it ("ZN "), and a key
            # that includes the padding reports a kept zinc as deleted.
            residue = line[17:20].strip().rjust(3) + line[20:26]
            present.add(residue)

            if centre:
                try:
                    position = (float(line[30:38]), float(line[38:46]),
                                float(line[46:54]))
                except ValueError:
                    continue
                if all(abs(position[axis] - centre[axis]) <= half
                       for axis in (0, 1, 2)):
                    in_box.add(residue)

    return present, in_box


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    out_dir, name, box_size, species = argv
    stem = os.path.join(out_dir, name)

    with open(f"{stem}.centre.tsv", encoding="utf-8") as handle:
        field = handle.readline().rstrip("\n").split("\t")
    centre = (float(field[0]), float(field[1]), float(field[2]))
    basis = field[3]

    before, in_box = residues_of(f"{stem}.clean.pdb", centre,
                                 float(box_size) / 2.0)
    after, _ = residues_of(f"{stem}.pdbqt")
    dropped = sorted(before - after)

    if dropped:
        with open(f"{stem}.dropped.txt", "w", encoding="utf-8") as record:
            for residue in dropped:
                marker = "\tIN THE BOX" if residue in in_box else ""
                record.write(f"{residue}{marker}\n")

        print(f"    Meeko could not type {len(dropped)} residue(s) and "
              f"removed them:", file=sys.stderr)
        for residue in dropped:
            print(f"      {residue}"
                  f"{'   <- IN THE BOX' if residue in in_box else ''}",
                  file=sys.stderr)
        if any(residue in in_box for residue in dropped):
            print("    ^ at least one was inside the docking box — this "
                  "changes the result", file=sys.stderr)
    elif os.path.exists(f"{stem}.dropped.txt"):
        os.remove(f"{stem}.dropped.txt")

    print("\t".join([name, species, f"{centre[0]:.3f}", f"{centre[1]:.3f}",
                     f"{centre[2]:.3f}", box_size, basis, str(len(dropped)),
                     f"{stem}.pdbqt"]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
