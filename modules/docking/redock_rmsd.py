#!/usr/bin/env python3
"""Ask whether a docked pose is where the crystal ligand was.

    redock_rmsd.py <poses.pdbqt> <crystal.sdf>

Prints one TSV row, with no header, to be appended to the score row:

    rmsd_top  rmsd_best  redock

      rmsd_top   RMSD in angstroms of the top-ranked pose from the crystal pose
      rmsd_best  the smallest RMSD of any pose Vina kept
      redock     pass, ranked-wrong or fail, read as below

The two columns separate two different failures:

    pass          the top pose is within 2.0 A of the crystal pose. The
                  protocol found a known answer, so its other scores are worth
                  reading.
    ranked-wrong  the top pose is not, but a lower-ranked one is. Vina searched
                  well enough to find the pose and the scoring function ranked
                  it below something wrong. The search is fine; the scores of
                  this pair should not be trusted for ranking.
    fail          no pose Vina kept is within 2.0 A. The pose was never found:
                  the box, the receptor preparation or the ligand's chemistry
                  is wrong, or the search was too short. Do not interpret any
                  score from this receptor until that is fixed.

If <crystal.sdf> does not exist, this ligand was not the one that crystallised
in this receptor, there is no answer to compare against, and the row is three
dashes.

The comparison is symmetry-aware and does no superposition: the pose is
already in the crystal's coordinate frame, so a ligand that is the right shape
in the wrong place scores badly, as it should. Hydrogens are ignored.

2.0 A is the conventional line, not a law. A small, rigid ligand such as
mannose that misses by 2.5 A has still missed; a 40-atom dinucleotide that
lands at 2.2 A has found the site and mis-placed a flexible tail.
"""
import os
import subprocess
import sys
import tempfile

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdMolAlign

RDLogger.DisableLog("rdApp.*")      # RDKit's warnings go to stderr unasked

NONE = "-\t-\t-"
CUTOFF = 2.0


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    poses, crystal = argv
    if not os.path.exists(crystal):
        print(NONE)
        return 0

    reference = Chem.MolFromMolFile(crystal, removeHs=True)
    if reference is None:
        print(f"WARNING: cannot read {crystal}", file=sys.stderr)
        print(NONE)
        return 0

    # Vina writes PDBQT, which carries no bond orders. Meeko wrote the
    # molecule's real chemistry into the ligand's PDBQT as a remark, and
    # mk_export.py reads it back and writes the poses out as an SDF.
    with tempfile.TemporaryDirectory() as folder:
        exported = os.path.join(folder, "poses.sdf")
        result = subprocess.run(["mk_export.py", poses, "-s", exported],
                                capture_output=True, text=True)
        if result.returncode != 0 or not os.path.exists(exported):
            print(f"WARNING: mk_export.py failed on {poses}: "
                  f"{result.stderr.strip()[-300:]}", file=sys.stderr)
            print(NONE)
            return 0

        values = []
        for pose in Chem.SDMolSupplier(exported, removeHs=True):
            try:
                values.append(rdMolAlign.CalcRMS(pose, reference))
            except (RuntimeError, ValueError, TypeError) as problem:
                print(f"WARNING: no RMSD for a pose of {poses}: {problem}",
                      file=sys.stderr)
                print(NONE)
                return 0

    if not values:
        print(NONE)
        return 0

    top, best = values[0], min(values)
    if top <= CUTOFF:
        verdict = "pass"
    elif best <= CUTOFF:
        verdict = "ranked-wrong"
    else:
        verdict = "fail"

    print(f"    redock {os.path.basename(poses)}: top {top:.2f} A, "
          f"best {best:.2f} A -> {verdict}", file=sys.stderr)
    print(f"{top:.2f}\t{best:.2f}\t{verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
