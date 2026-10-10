#!/usr/bin/env python3
"""Write the SWISS-MODEL uploads, one per target, and the manifest.

    swissmodel_upload.py <prep_dir> <out_dir> <drop_dir> <warn_length>

Writes <out_dir>/<organism>/<target>.fasta for every prepared target, plus
submission_manifest.tsv with one row each:

    target  species  length  reference  locus_tag  pct_identity
    upload_file  drop_model_into

The folders mirror program C's: one per reference organism, named the same
way. SWISS-MODEL is a web form, so this part of the work is done by hand —
fifteen uploads, then fifteen downloads, then the numbers copied off each
results page. Doing that an organism at a time is the only way it stays
straight, and a folder per organism is what makes that possible. Each row's
drop_model_into is the matching folder under models/, so a download goes back
where its upload came from.

The upload is a plain FASTA with the target name as its only header. Program
C's header is self-documenting and useful on disk, but SWISS-MODEL shows the
header back to you as the project title, and a line of key=value pairs makes
an unreadable one.

The evidence that made this a target is carried into the manifest instead, so
the table you fill in later already says where the sequence came from.

Two sequence identities will be in play and they are not the same number. The
pct_identity here is your metagenome sequence against the reference gene,
from program A's BLAST — that is what makes this a target. The identity
SWISS-MODEL reports is your sequence against its template, and that is what
decides whether the model is any good. Keep them apart.
"""
import glob
import os
import sys

import fasta

HEADER = ["target", "species", "length", "reference", "locus_tag",
          "pct_identity", "upload_file", "drop_model_into"]


def evidence(report):
    """Reference gene, locus tag and BLAST identity per target, from program C."""
    found = {}
    with open(report, encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            field = line.rstrip("\n").split("\t")
            if index and len(field) >= 6:
                found[field[0]] = (field[4], field[3], field[5])
    return found


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    prep_dir, out_dir, drop_dir, warn_length = argv
    warn_length = int(warn_length)

    known = evidence(os.path.join(prep_dir, "prep_report.tsv"))

    # One organism folder deep, which is how program C writes them. The
    # organism comes from the folder rather than the report, so the uploads
    # are grouped exactly the way the sequences on disk are.
    paths = sorted(glob.glob(os.path.join(prep_dir, "targets", "*", "*.faa")))
    if not paths:
        print(f"ERROR: no prepared sequences at {prep_dir}/targets/*/",
              file=sys.stderr)
        print("Run:  bash C_protein_prep_python.sh", file=sys.stderr)
        return 1

    manifest = open(os.path.join(out_dir, "submission_manifest.tsv"), "w",
                    encoding="utf-8")
    manifest.write("\t".join(HEADER) + "\n")
    notes = []

    for path in paths:
        target = os.path.basename(path)[:-4]
        species = os.path.basename(os.path.dirname(path))

        with open(path, encoding="utf-8", errors="replace") as handle:
            sequence = "".join(text for _, text in fasta.records(handle))

        # The folder the model comes back to is made now, empty, rather than
        # when the download arrives: an empty folder per organism is the
        # instruction, and program F reads the same folders.
        os.makedirs(os.path.join(out_dir, species), exist_ok=True)
        os.makedirs(os.path.join(drop_dir, species), exist_ok=True)

        upload = os.path.join(species, f"{target}.fasta")
        with open(os.path.join(out_dir, upload), "w",
                  encoding="utf-8") as handle:
            handle.write(fasta.record(target, sequence))

        length = len(sequence)
        reference, locus, identity = known.get(target, ("-", "-", "-"))

        manifest.write("\t".join([target, species, str(length), reference,
                                  locus, identity, upload,
                                  f"{drop_dir}/{species}/"]) + "\n")

        if length > warn_length:
            notes.append(f"NOTE: {target} is {length} residues — expect "
                         f"partial coverage; consider submitting one domain")

    manifest.close()
    print(f"{len(paths)} uploads written to {out_dir}/, one folder per organism")
    for note in notes:
        print(note)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
