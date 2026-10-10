#!/usr/bin/env python3
"""Gather the downloaded models into one folder, under one naming scheme.

    collect_models.py <models_dir> <targets_dir> <swissmodel_drop> <af3_drop>

Copies every .pdb and .cif from the two drop folders into <models_dir> as
SWISSMODEL_<target>.<ext> or AF3_<target>.<ext>, and prints what it did.

The drop folder is what identifies the service, which is why you do not have
to rename anything: drop SWISS-MODEL's model_01.pdb and AlphaFold's
fold_exot_hjbccn_007355_model_0.cif in as they download. That identification
matters more than it looks — the two services write different quantities into
the same B-factor column, QMEANDisCo from 0 to 1 against pLDDT from 0 to 100,
and nothing inside either file says which one you have.

Inside a drop folder there is one folder per reference organism, the same
folders and the same names as program C's, because that is where programs D
and E said to put each download. <models_dir> is flat on purpose: it is what
programs G and I dock, the organism is already in every report, and a docking
run cares which model it has, not which organism's folder it arrived in.

Which target a download belongs to is recovered from the file name, because
it is the only link back to the gene. The services mangle the job name their
own way — AlphaFold lowercases it and wraps it in fold_..._model_0 — so the
test is whether a target name appears anywhere inside the file name. A file
no target matches is still copied, under its own name, and said out loud:
usually it is a download from a job this pipeline did not create.

<models_dir> is derived and yours to delete. A list of what this script put
there last time is kept in .collected, and only those files are cleared
before copying — so a download you have since deleted stops being docked,
while a model you placed in the folder by hand is left alone.
"""
import glob
import os
import shutil
import sys

SERVICES = (("SWISSMODEL", 2), ("AF3", 3))      # prefix, position in argv


def target_names(targets_dir):
    """Every prepared target, longest first so the most specific name wins."""
    names = [os.path.basename(path)[:-4]
             for path in sorted(glob.glob(os.path.join(targets_dir, "*", "*.faa")))]
    return sorted(names, key=len, reverse=True)


def species_names(targets_dir):
    """The organism folders program C grouped the targets into."""
    return sorted(os.path.basename(path)
                  for path in glob.glob(os.path.join(targets_dir, "*"))
                  if os.path.isdir(path))


def downloads(drop, species):
    """Every model file in a drop folder: one organism deep, and loose.

    Two levels and no further, which is a choice worth stating. An unzipped
    AlphaFold job folder holds five models of the same target and its template
    hits besides; SWISS-MODEL's mmcif/ holds the model that was already
    collected as .pdb. Searching all the way down would pool every one of them
    and hand the lot to the docking programs. A file left loose in the drop
    folder is still taken, because that is where downloads landed before there
    were folders per organism.
    """
    paths = []
    for folder in [drop] + [os.path.join(drop, name) for name in species]:
        paths += glob.glob(os.path.join(folder, "*.pdb"))
        paths += glob.glob(os.path.join(folder, "*.cif"))
    return sorted(paths)


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    models_dir, targets_dir = argv[0], argv[1]
    manifest = os.path.join(models_dir, ".collected")
    names = target_names(targets_dir)
    species = species_names(targets_dir)

    if os.path.exists(manifest):
        with open(manifest, encoding="utf-8") as handle:
            for line in handle:
                stale = line.strip()
                if stale:
                    # Only inside models_dir, and only a plain file name: this
                    # list is read back from disk and used to delete.
                    path = os.path.join(models_dir, os.path.basename(stale))
                    if os.path.isfile(path):
                        os.remove(path)

    collected = []
    with open(manifest, "w", encoding="utf-8") as record:
        for prefix, position in SERVICES:
            for path in downloads(argv[position], species):
                base = os.path.basename(path)
                stem, extension = os.path.splitext(base)

                target = ""
                for name in names:
                    if name.lower() in stem.lower():
                        target = name
                        break

                if target:
                    destination = f"{prefix}_{target}{extension}"
                else:
                    destination = f"{prefix}_{stem}{extension}"
                    print(f"  no target matched: {base} — kept as {destination}")

                shutil.copyfile(path, os.path.join(models_dir, destination))
                record.write(destination + "\n")
                collected.append(destination)
                print(f"  {base}  ->  {destination}")

    print(f"collected {len(collected)} file(s) into {models_dir}/")

    # Targets that went out and have not come back. Worth naming, because an
    # empty return from SWISS-MODEL is a finding about the PDB and is easy to
    # mistake for a step you forgot to run.
    present = [os.path.basename(path)
               for path in glob.glob(os.path.join(models_dir, "*"))]
    missing = [name for name in sorted(names)
               if not any(name in file_name for file_name in present)]
    if missing:
        print("\nsubmitted but no model collected yet:")
        for name in missing:
            print(f"  {name}")
        print("  (a SWISS-MODEL target with no template returns nothing — "
              "record that)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
