#!/usr/bin/env python3
"""Generate WORKING_FOLDER/INPUT_DATASETS/MANIFEST.tsv from the files on disk.

The manifest is an inventory of the curated dataset: one row per data file,
with its length and a checksum. Deriving it from the filesystem is the point —
a hand-maintained inventory drifts from the directory it describes, and this
one did: it was missing a file that was on disk, and its paths were written
relative to a directory the species folders had since moved out of.

PATHS ARE RELATIVE TO THE MANIFEST
----------------------------------
Every ``path`` resolves from the directory MANIFEST.tsv sits in, so
``from_Database/1_Mycobacterium_tuberculosis/FASTA_O05442_CpnT.fasta`` is
openable as-is from WORKING_FOLDER/INPUT_DATASETS/. That is what makes a
manifest usable by a script that is handed only the manifest.

WHAT IS COUNTED
---------------
The curated files — FASTA_ (protein sequences), GENE_ (nucleotide CDS) and
PDB_ (experimental structures). ``MODEL_*.pdb`` files are predicted structures,
not curated downloads, and are deliberately excluded; so are the ``.orig``
backups that fetch_gene_cds.py leaves behind.

COLUMN CONVENTIONS — these match the manifest as it was first written
---------------------------------------------------------------------
length/units   FASTA -> residues, "aa"; GENE -> bases, "nt";
               PDB   -> lines beginning "ATOM", "ATOM_records" (HETATM
               records, which are ligands and waters, are not counted).
sha256_file    the first 16 hex characters of the SHA-256 of the file's bytes.

Usage (from the repository root):

    python modules/datasets/write_manifest.py            # show what would change
    python modules/datasets/write_manifest.py --write    # rewrite MANIFEST.tsv

Standard library only - no extra conda packages required.
"""

import argparse
import hashlib
import os
import re
import sys

COLUMNS = ("path", "category", "molecule", "length", "units", "sha256_file")

# Prefix -> (category, molecule, units). Order is the order rows are grouped
# within a species folder; MODEL_ is absent on purpose, see the module docstring.
KINDS = {
    "FASTA_": ("sequence", "protein", "aa"),
    "GENE_": ("sequence", "DNA", "nt"),
    "PDB_": ("structure", "PDB", "ATOM_records"),
}

SPECIES = re.compile(r"^\d+_")


def digest(path):
    """First 16 hex characters of the file's SHA-256."""
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()[:16]


def sequence_length(path):
    """Residues or bases in a one-record FASTA, ignoring the header."""
    with open(path, encoding="utf-8") as handle:
        return sum(len(line.strip()) for line in handle
                   if not line.startswith(">"))


def atom_records(path):
    """Lines beginning ATOM. HETATM is not an ATOM record."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.startswith("ATOM"))


def kind_of(name):
    """The KINDS entry for a filename, or None if it is not a curated file."""
    for prefix, spec in KINDS.items():
        if name.startswith(prefix):
            return spec
    return None


def rows(datasets):
    """One row per curated file, sorted by path, as a list of tuples."""
    found = []
    root = os.path.join(datasets, "from_Database")
    for species in sorted(os.listdir(root)):
        folder = os.path.join(root, species)
        if not SPECIES.match(species) or not os.path.isdir(folder):
            continue
        for name in sorted(os.listdir(folder)):
            if name.endswith(".orig"):
                continue
            spec = kind_of(name)
            if spec is None:
                continue
            category, molecule, units = spec
            full = os.path.join(folder, name)
            length = (atom_records(full) if units == "ATOM_records"
                      else sequence_length(full))
            # Relative to the manifest's own directory, see the docstring.
            path = "/".join(("from_Database", species, name))
            found.append((path, category, molecule, str(length), units,
                          digest(full)))
    return sorted(found)


def existing(manifest):
    """The manifest's current data rows, or [] if it is not there yet."""
    if not os.path.isfile(manifest):
        return []
    with open(manifest, encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    return [tuple(line.split("\t")) for line in lines[1:] if line.strip()]


def describe(before, after):
    """Print the difference between two row lists; return True if identical."""
    was = {row[0]: row for row in before}
    now = {row[0]: row for row in after}
    added = sorted(set(now) - set(was))
    removed = sorted(set(was) - set(now))
    changed = sorted(p for p in set(was) & set(now) if was[p] != now[p])

    for path in added:
        print(f"  + {path}")
    for path in removed:
        print(f"  - {path}")
    for path in changed:
        old_fields = dict(zip(COLUMNS, was[path]))
        new_fields = dict(zip(COLUMNS, now[path]))
        diffs = ", ".join(f"{c}: {old_fields[c]} -> {new_fields[c]}"
                          for c in COLUMNS[1:]
                          if old_fields.get(c) != new_fields.get(c))
        print(f"  ~ {path}  ({diffs})")

    print(f"\n{len(before)} rows before, {len(after)} after — "
          f"{len(added)} added, {len(removed)} removed, {len(changed)} changed")
    return not (added or removed or changed)


def main():
    parser = argparse.ArgumentParser(
        description="Generate MANIFEST.tsv from the dataset files on disk.")
    parser.add_argument("--write", action="store_true",
                        help="rewrite MANIFEST.tsv (default: show the diff)")
    parser.add_argument("--datasets", default="WORKING_FOLDER/INPUT_DATASETS",
                        help="directory holding MANIFEST.tsv and from_Database/")
    args = parser.parse_args()

    manifest = os.path.join(args.datasets, "MANIFEST.tsv")
    after = rows(args.datasets)
    if not after:
        print(f"No curated files found under {args.datasets}/from_Database/.")
        return 1

    counts = {}
    for row in after:
        counts[row[2]] = counts.get(row[2], 0) + 1
    summary = ", ".join(f"{n} {molecule}" for molecule, n in sorted(counts.items()))
    print(f"{len(after)} curated files on disk ({summary})\n")

    identical = describe(existing(manifest), after)

    if not args.write:
        if identical:
            print("MANIFEST.tsv is up to date.")
            return 0
        print("Re-run with --write to apply.")
        return 0

    with open(manifest, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in after:
            handle.write("\t".join(row) + "\n")
    print(f"\nwrote {manifest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
