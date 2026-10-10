#!/usr/bin/env python3
"""Read the best score out of a Vina pose file, with its box.

    pose_score.py <receptors.tsv> <receptor> <ligand> <pose.pdbqt>

Prints one TSV row for the calling script's score table:

    receptor  ligand  best_kcal_mol  rmsd_lb  rmsd_ub  basis  box_size
    pose_file

Vina writes each pose's score into the pose file itself, on a line reading
"REMARK VINA RESULT:". The first one is the best, and the two RMSDs are that
pose's distance from it, so they are zero for mode 1 by definition.

The basis and the box edge are copied from program G's receptors.tsv rather
than passed in again, so the row says which box produced the number. That
matters more than it looks: a score is comparable to another score only
within one receptor and one box, and a table of scores with no box recorded
invites exactly the comparison that cannot be made.
"""
import os
import sys

NA = "-"


def box_of(table, receptor):
    """The basis and box edge program G recorded for this receptor.

    Columns of receptors.tsv, which receptor_report.py writes:
    receptor, species, centre x, y, z, box_size, basis, dropped, pdbqt.
    """
    with open(table, encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            field = line.rstrip("\n").split("\t")
            if index and len(field) >= 7 and field[0] == receptor:
                return field[6], field[5]
    return NA, NA


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    table, receptor, ligand, pose = argv
    basis, box_size = box_of(table, receptor)

    score, rmsd_lb, rmsd_ub = NA, NA, NA
    with open(pose, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("REMARK VINA RESULT"):
                # "REMARK VINA RESULT:   -5.2   0.000   0.000"
                field = line.split()
                score, rmsd_lb, rmsd_ub = field[3], field[4], field[5]
                break

    if score == NA:
        print(f"WARNING: no score in {pose}", file=sys.stderr)
    else:
        print(f"    {receptor} x {ligand}:  best {score} kcal/mol",
              file=sys.stderr)

    print("\t".join([receptor, ligand, score, rmsd_lb, rmsd_ub, basis,
                     box_size, os.path.join("poses", os.path.basename(pose))]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
