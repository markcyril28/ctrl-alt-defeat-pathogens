#!/usr/bin/env python3
"""Split program B's grouped proteins into one clean sequence per target.

    prep_proteins.py <extraction_dir> <out_dir> <min_aa> <warn_aa> \\
                     <drop_duplicates>

Writes into <out_dir>, grouped by the organism each reference gene came from:

    targets/<organism>/<target>.faa  one sequence each, which is what the
                                     services want, in the organism's folder
    by_species/<organism>.faa        that organism's targets in one file, which
                                     is what program E batches per job upload
    all_targets.faa                  every organism in one file, for grepping a
                                     locus tag by hand or for a tool that takes
                                     one input
    prep_report.tsv                  what was changed, per target, and the
                                     verdict

Why the folders. These targets come from four reference organisms, and a flat
folder of fifteen files mixes them: the next thing anyone does by hand is
submit, download and compare one organism at a time. The organism is already
in the reference gene's name, "3_Klebsiella_pneumoniae__FimH", so the grouping
needs no lookup and no guessing. The rule itself lives in species_of.py,
because programs G and J have to arrive at the same answer for a model whose
name no longer says which gene it came from.

What it does not mean is that the sequence came from that organism. The
organism folder here says which reference gene the sequence matched. A locus
tag belongs to the community, not to a species — this pipeline does no
binning, so a hit to a Klebsiella gene says the sequence is in the sample, not
which organism carries it.

What it cleans, and why each one is worth doing:

    Stop characters. A "*" inside a sequence is not a formatting quirk, it is
    a stop codon in the middle of the open reading frame — the frame is wrong,
    or the gene really is truncated. Both services reject the character, so it
    is removed, but the count is reported: a sequence that needed several
    removed is a sequence to look at before you model it.

    Non-standard residues. Anything outside the twenty standard amino acids
    becomes X. Bakta writes X where the DNA had an ambiguous base, and an X is
    honest — it says the residue is unknown. What it is not is modellable:
    both services will place something there anyway, and that something is an
    invention.

    Duplicates. Identical sequences are collapsed to one, because paying a
    service twice for the same answer tells you nothing twice. This is not
    hypothetical: a metagenome assembly often carries the same region on two
    contigs, and program B reports both — correctly, since they are two
    separate pieces of evidence. One of them is enough to model.

    Length. Shorter than min_aa is not a protein, it is a fragment of one, and
    a fold predicted from it is a guess with nothing to constrain it. Longer
    than warn_aa is flagged rather than dropped: it will model, slowly, and
    usually better one domain at a time.

The verdict column says whether a sequence can be submitted, and nothing at
all about whether it is the gene you want.
"""
import glob
import os
import sys

import fasta
import species_of

STANDARD = "ACDEFGHIKLMNPQRSTVWY"

UNASSIGNED = species_of.UNASSIGNED

HEADER = [
    "target", "species", "assembly", "locus_tag", "reference", "pct_identity",
    "aa_in", "aa_out", "stops_removed", "nonstandard", "duplicate_of",
    "verdict", "product",
]


def read_targets(path):
    """Program B's table, in order, as one dictionary per target."""
    targets = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            field = line.rstrip("\n").split("\t")
            if index == 0 or len(field) < 15:
                continue
            targets.append({
                "target": field[0], "assembly": field[1], "locus": field[2],
                "reference": field[3], "pident": field[4],
                "species": species_of.of_reference(field[3]),
                "ref_range": f"{field[6]}-{field[7]}",
                "location": f"{field[9]}:{field[10]}-{field[11]}({field[12]})",
                "product": field[14],
            })
    return targets


def read_sequences(folder):
    """Every sequence in every .faa of a folder, keyed by its first header word."""
    sequences = {}
    for path in sorted(glob.glob(os.path.join(folder, "*.faa"))):
        name = None
        with open(path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.rstrip("\r\n")
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    sequences.setdefault(name, "")
                elif name:
                    sequences[name] += line.replace(" ", "").replace("\t", "")
    return sequences


def clean(raw):
    """Upper case, stops removed, unknown residues as X. Returns the counts."""
    raw = raw.upper()
    stops = raw.count("*")
    raw = raw.replace("*", "")

    residues, nonstandard = [], 0
    for residue in raw:
        if residue not in STANDARD:
            nonstandard += 1
            residue = "X"
        residues.append(residue)
    return "".join(residues), stops, nonstandard


def write_grouped(out_dir, prepared):
    """The three views of the same sequences: per target, per species, pooled.

    Three files for one sequence looks like duplication, and it is — each one
    is there because something downstream reads that shape. Program D uploads
    one target at a time, program E sends one organism per batch, and a single
    pooled file is what you grep when you want to find a locus tag by hand.
    """
    for species, target, record in prepared:
        folder = os.path.join(out_dir, "targets", species)
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, f"{target}.faa"), "w",
                  encoding="utf-8") as one:
            one.write(record)

    by_species = {}
    for species, _, record in prepared:
        by_species.setdefault(species, []).append(record)

    os.makedirs(os.path.join(out_dir, "by_species"), exist_ok=True)
    pooled = []
    for species in sorted(by_species):
        records = by_species[species]
        # A blank line between records and none at the end, which is how every
        # multi-record file in this pipeline is written.
        text = "\n".join(records)
        with open(os.path.join(out_dir, "by_species", f"{species}.faa"), "w",
                  encoding="utf-8") as group:
            group.write(text)
        pooled.append(text)

    with open(os.path.join(out_dir, "all_targets.faa"), "w",
              encoding="utf-8") as everything:
        everything.write("\n".join(pooled))

    return by_species


def main(argv):
    if len(argv) != 5:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    extraction, out_dir, min_aa, warn_aa, drop_duplicates = argv
    min_aa, warn_aa = int(min_aa), int(warn_aa)
    drop_duplicates = drop_duplicates == "true"

    targets = read_targets(os.path.join(extraction, "targets.tsv"))
    sequences = read_sequences(os.path.join(extraction, "proteins"))

    if not sequences:
        print(f"ERROR: {extraction}/proteins/ holds no sequences", file=sys.stderr)
        print("Run:  bash B_gene_protein_extraction_bakta.sh", file=sys.stderr)
        return 1

    report = open(os.path.join(out_dir, "prep_report.tsv"), "w",
                  encoding="utf-8")
    report.write("\t".join(HEADER) + "\n")

    first_seen, notes, kept = {}, [], 0
    prepared = []       # (species, target, record) for everything submittable
    missing, duplicates, too_short, unknown, too_long = [], [], [], [], []

    for target in targets:
        name = target["target"]
        raw = sequences.get(name, "")

        if not raw:
            missing.append(name)
            report.write("\t".join([
                name, target["species"], target["assembly"], target["locus"],
                target["reference"], target["pident"], "0", "0", "-", "-", "-",
                "no_sequence", target["product"],
            ]) + "\n")
            continue

        sequence, stops, nonstandard = clean(raw)

        # Identical sequence, already seen: the first one keeps the slot. The
        # comparison is across every organism, not within one, because two
        # reference genes can be close enough that one contig matches both.
        duplicate_of = "-"
        if sequence in first_seen and drop_duplicates:
            duplicate_of = first_seen[sequence]
        elif sequence not in first_seen:
            first_seen[sequence] = name

        if duplicate_of != "-":
            verdict = "duplicate"
            duplicates.append(f"  {name} is identical to {duplicate_of}")
        elif len(sequence) < min_aa:
            verdict = "too_short"
            too_short.append(f"  {name} ({len(sequence)} aa)")
        else:
            verdict = "ok"
            kept += 1

            if len(sequence) > warn_aa:
                too_long.append(f"  {name} ({len(sequence)} aa)")
            if nonstandard:
                unknown.append(f"  {name} ({nonstandard} of {len(sequence)} "
                               f"unknown)")

            # The header carries its own provenance, so a file that ends up
            # somewhere on its own still says where it came from — the folder
            # it sits in is not part of the file. Program D replaces the header
            # with the bare name, because SWISS-MODEL shows the header back as
            # the project title.
            header = (f"{name} locus={target['locus']}"
                      f" assembly={target['assembly']}"
                      f" species={target['species']}"
                      f" reference={target['reference']}"
                      f" pident={target['pident']}"
                      f" ref_range={target['ref_range']}"
                      f" location={target['location']}"
                      f" aa={len(sequence)}")
            prepared.append((target["species"], name,
                             fasta.record(header, sequence)))

        report.write("\t".join([
            name, target["species"], target["assembly"], target["locus"],
            target["reference"], target["pident"], str(len(raw)),
            str(len(sequence)), str(stops),
            str(nonstandard) if nonstandard else "-", duplicate_of, verdict,
            target["product"],
        ]) + "\n")

    report.close()
    by_species = write_grouped(out_dir, prepared)

    print(f"{len(targets)} targets from program B, {kept} ready to submit")
    for species in sorted(by_species):
        print(f"  {species:<32} {len(by_species[species])}")

    # The notes go to stderr so the table the calling script prints stays a
    # clean TSV.
    if missing:
        notes.append("no sequence found for: " + ", ".join(missing))
    for title, lines in (("identical sequences, one kept", duplicates),
                         (f"dropped as too short (under {min_aa} aa)", too_short),
                         ("contains unknown residues, written as X", unknown),
                         (f"longer than {warn_aa} aa — consider one domain per "
                          f"job", too_long)):
        if lines:
            notes.append(title + ":\n" + "\n".join(lines))
    if UNASSIGNED in by_species:
        notes.append(f"grouped under {UNASSIGNED}/: the reference gene name "
                     f"carries no organism, so there was nothing to group on")
    if not kept:
        notes.append("nothing passed, so there is nothing to submit.")

    for note in notes:
        print(note, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
