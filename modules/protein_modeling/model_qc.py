#!/usr/bin/env python3
"""Measure every model in a folder, and say what its B-factor column means.

    model_qc.py <models_dir> [--conf-min N] [--convert]

      --conf-min N   confidence below this counts as low (default 70, on the
                     0-100 scale)
      --convert      also write a .pdb beside every .cif, for the tools later
                     in the pipeline that do not read mmCIF

Prints a TSV to stdout, one row per model.

The B-factor column is the trap this script exists for. In a crystal
structure it is a crystallographic B-factor. In an AlphaFold model it is
pLDDT, per-residue confidence from 0 to 100, where high is good — the
opposite direction to a B-factor. In a SWISS-MODEL file it is usually
QMEANDisCo local, from 0 to 1. Three different quantities in the same column,
and nothing in the file announces which one you have.

So the scale is detected from the values and reported next to the mean, and
the low-confidence fraction is computed after putting everything on the 0-100
scale. A single mean over a whole chain still hides the thing that matters
most for docking: confidence over the site you intend to dock into. Open the
model and colour by B-factor before you trust any number here.
"""
import glob
import os
import sys

from pymol import cmd

METALS = {"ZN", "MG", "MN", "CA", "FE", "CU", "NI", "CO", "NA", "K", "FE2", "FES"}
SOLVENT = {"HOH", "WAT", "DOD"}

HEADER = [
    "model", "source", "format", "chains", "residues", "atoms", "metals",
    "ligands", "mean_conf", "conf_scale", "pct_low_conf", "verdict",
]


def _source(stem):
    """Where the model came from, by the naming convention programs E and F set."""
    upper = stem.upper()
    if upper.startswith("AF3"):
        return "alphafold3"
    if upper.startswith("SWISSMODEL"):
        return "swiss-model"
    return "other"


def _confidence(values, source):
    """(mean, scale label, scale factor to 0-100) for a list of B-factors.

    The values alone cannot settle this. A crystallographic B-factor usually
    falls between 10 and 80, and so does pLDDT — identical ranges, opposite
    meanings, same column. So the source vouches for the scale, and where it
    cannot, the label says "unverified" rather than picking the flattering
    reading. A real structure read as pLDDT comes out looking like a terrible
    model; a model read as a B-factor comes out looking like a good crystal.
    """
    if not values:
        return None, "none", 1.0

    mean = sum(values) / len(values)
    high = max(values)

    if high > 100.0:
        # Out of range for either confidence score, so it is a B-factor.
        return mean, "b_factor", 0.0

    if high <= 1.5:
        # 0 to 1: QMEANDisCo local, which is what SWISS-MODEL writes.
        return mean, "0-1_qmean" if source != "alphafold3" else "0-1_unexpected", 100.0

    # 0 to 100: pLDDT if an AlphaFold job produced it. From anywhere else —
    # a crystal structure dropped into a models folder, a file renamed by
    # hand — it is not safe to call this confidence at all.
    if source == "alphafold3":
        return mean, "0-100_plddt", 1.0
    return mean, "0-100_unverified", 1.0


def main(argv):
    if not argv:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    models_dir = argv[0]
    conf_min, convert = 70.0, False

    rest = argv[1:]
    while rest:
        flag = rest.pop(0)
        if flag == "--conf-min":
            conf_min = float(rest.pop(0))
        elif flag == "--convert":
            convert = True
        else:
            print(f"unknown option: {flag}", file=sys.stderr)
            return 2

    paths = sorted(
        glob.glob(os.path.join(models_dir, "*.pdb"))
        + glob.glob(os.path.join(models_dir, "*.cif"))
    )

    print("\t".join(HEADER))

    # A .cif converted on an earlier run leaves a .pdb of the same name, and
    # measuring both would double every row.
    stems_seen = set()

    # Two things worth saying out loud after the table, rather than leaving
    # them to be spotted in a column.
    low_confidence, no_metal = [], []

    for path in paths:
        stem, extension = os.path.splitext(os.path.basename(path))
        if stem in stems_seen:
            continue
        stems_seen.add(stem)

        cmd.delete("all")
        try:
            cmd.load(path, "m")
        except Exception as error:                      # a truncated download
            print(f"{stem}\t{_source(stem)}\t{extension.lstrip('.')}\t"
                  f"-\t-\t-\t-\t-\t-\t-\t-\tunreadable", file=sys.stdout)
            print(f"could not load {path}: {error}", file=sys.stderr)
            continue

        residues = cmd.count_atoms("m and polymer and name CA")
        atoms = cmd.count_atoms("m")
        chains = len(cmd.get_chains("m"))

        hetero = set()
        cmd.iterate("m and not polymer", "hetero.add(resn)", space={"hetero": hetero})
        hetero -= SOLVENT
        metals = sorted(hetero & METALS)
        ligands = sorted(hetero - METALS)

        values = []
        cmd.iterate("m and polymer and name CA", "values.append(b)",
                    space={"values": values})
        mean, scale, factor = _confidence(values, _source(stem))

        if mean is None:
            mean_text, pct_low, verdict = "-", "-", "no_confidence_column"
        elif scale == "b_factor":
            mean_text, pct_low = f"{mean:.1f}", "-"
            verdict = "experimental_b_factor"
        else:
            scaled = [v * factor for v in values]
            low = sum(1 for v in scaled if v < conf_min)
            mean_text = f"{mean * factor:.1f}"
            pct_low = f"{100.0 * low / len(scaled):.1f}"
            if scale.endswith("_unverified") or scale.endswith("_unexpected"):
                # The numbers are reported, but nothing vouches for what they
                # mean, so they must not be read as a verdict on the model.
                verdict = "check_confidence_scale"
            else:
                verdict = "ok" if low / len(scaled) < 0.25 else "low_confidence"

        if convert and extension == ".cif":
            # Overwritten, not skipped. The .pdb is derived from the .cif, and
            # a re-downloaded model keeps its name: skipping an existing file
            # would leave the QC table describing the new model while programs
            # H and J went on docking the old one.
            cmd.save(os.path.join(models_dir, f"{stem}.pdb"), "m")

        print("\t".join([
            stem, _source(stem), extension.lstrip("."), str(chains),
            str(residues), str(atoms),
            ",".join(metals) or "-", ",".join(ligands) or "-",
            mean_text, scale, pct_low, verdict,
        ]))

        if verdict == "low_confidence":
            low_confidence.append(f"{stem} ({pct_low}% below {conf_min:g})")
        if not metals and _source(stem) != "other":
            no_metal.append(stem)

    # To stderr, so stdout stays a clean TSV the calling script can redirect.
    if low_confidence:
        print("\nlow confidence over a quarter of the chain:", file=sys.stderr)
        for note in low_confidence:
            print(f"  {note}", file=sys.stderr)
    if no_metal:
        print("\nno metal in the model — if the enzyme needs one, program G "
              "must put it back:", file=sys.stderr)
        for note in no_metal:
            print(f"  {note}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
