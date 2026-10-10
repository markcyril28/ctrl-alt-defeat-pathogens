#!/usr/bin/env python3
"""Pull the metagenome's own copy of each matched gene out of the annotation.

    extract_targets.py <catalog_dir> <annotation_dir> <out_dir> \\
                       <references> <min_aa>

      <references>  reference genes to keep, space separated, by gene name
                    ("TNT_CpnT FimH") or by the full name program A gave them
                    ("1_Mycobacterium_tuberculosis__TNT_CpnT"). "" keeps all.
      <min_aa>      drop proteins shorter than this many residues

Writes into <out_dir>:

    targets.tsv              one row per target, with all of its evidence
    genes/<reference>.fna    the DNA, from Bakta's .ffn
    proteins/<reference>.faa the protein, from Bakta's .faa

Records in the two FASTA files are separated by a blank line, which is for
reading them: several genes under one reference, each with a long header, run
together otherwise.

The three Bakta files are all needed and none is optional: .ffn has the
nucleotide sequence of each feature, .faa has the translation, and the .tsv is
the only one that says which contig the gene sits on and which way round it
runs. All three are keyed by locus tag, which is what joins them.

That locus tag is also the weak point, and the reason this script refuses to
run on a stale catalog. Bakta assigns tags per run, so re-running program 09
renumbers every gene. The tags in program A's matches would still resolve —
to different genes — and nothing further down the pipeline could notice. So
the size and modification time of each annotation table, as program A
recorded them, are checked before anything is extracted.

Within each output file the genes are in the order they land on the
reference. When consecutive locus tags on one contig tile one reference — one
covering the start, the next the middle — they are pieces of one gene split by
frameshifts in the assembly, not several genes. The contig and the
coordinates in targets.tsv are how you tell.
"""
import os
import sys

import fasta

HEADER = [
    "target", "assembly", "locus_tag", "reference", "pct_identity",
    "aln_length", "ref_start", "ref_end", "pct_ref_covered",
    "contig", "start", "stop", "strand", "aa_length", "product",
]


def rows_of(path):
    """The data rows of a TSV, without its header, as lists of strings."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            if index and line.strip():
                yield line.rstrip("\n").split("\t")


def check_not_stale(catalog_dir):
    """Stop if an assembly was re-annotated after program A searched it."""
    source = os.path.join(catalog_dir, "source.tsv")
    for assembly, table, size, mtime in rows_of(source):
        if (not os.path.isfile(table)
                or str(os.path.getsize(table)) != size
                or str(int(os.path.getmtime(table))) != mtime):
            print(f"ERROR: {assembly} was re-annotated after program A ran "
                  f"({table} changed).", file=sys.stderr)
            print("Locus tags are assigned per Bakta run, so A's matches no "
                  "longer name these genes.", file=sys.stderr)
            print("Rebuild:  bash A_gene_catalog_blastn.sh", file=sys.stderr)
            return False
    return True


def safe_name(text):
    """Keep only characters a shell and a web upload form both accept."""
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in text)


def select(catalog_dir, references):
    """Program A's matches worth carrying forward, in reference order."""
    wanted = references.split()
    chosen = []

    for row in rows_of(os.path.join(catalog_dir, "matches.tsv")):
        reference = row[2]
        gene = reference.split("__", 1)[-1]
        if wanted and reference not in wanted and gene not in wanted:
            continue

        # A minus-strand match runs backwards along the reference. Stored low
        # to high so the pieces of one gene can be sorted along it.
        low, high = sorted((int(row[8]), int(row[9])))

        chosen.append({
            "target": safe_name(f"{gene}_{row[1]}"),
            "assembly": row[0], "locus": row[1], "reference": reference,
            "pident": row[4], "aln_length": row[5],
            "low": low, "high": high, "pct_ref_covered": row[12],
        })

    chosen.sort(key=lambda t: (t["reference"], t["assembly"], t["low"]))
    return chosen


def read_fasta(path, wanted):
    """The sequence of each wanted locus tag in a Bakta .ffn or .faa."""
    sequences, tag = {}, None
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                tag = line[1:].split()[0]
                tag = tag if tag in wanted else None
            elif tag:
                sequences[tag] = sequences.get(tag, "") + line
    return sequences


def read_locations(path, wanted):
    """Contig, start, stop, strand and product of each wanted locus tag."""
    locations = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("#"):
            continue
        field = line.rstrip("\n").split("\t")
        if len(field) >= 8 and field[5] in wanted:
            locations[field[5]] = (field[0], field[2], field[3], field[4],
                                   field[7])
    return locations


def main(argv):
    if len(argv) != 5:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    catalog_dir, annotation_dir, out_dir, references, min_aa = argv
    min_aa = int(min_aa)

    if not check_not_stale(catalog_dir):
        return 1

    targets = select(catalog_dir, references)
    total = sum(1 for _ in rows_of(os.path.join(catalog_dir, "matches.tsv")))

    # One read of each assembly's three Bakta files, not one per target.
    dna, protein, where = {}, {}, {}
    for assembly in sorted({t["assembly"] for t in targets}):
        wanted = {t["locus"] for t in targets if t["assembly"] == assembly}
        stem = os.path.join(annotation_dir, assembly, assembly)

        for extension in ("tsv", "ffn", "faa"):
            if not os.path.isfile(f"{stem}.{extension}"):
                print(f"ERROR: {stem}.{extension} is missing", file=sys.stderr)
                return 1

        dna[assembly] = read_fasta(f"{stem}.ffn", wanted)
        protein[assembly] = read_fasta(f"{stem}.faa", wanted)
        where[assembly] = read_locations(f"{stem}.tsv", wanted)

    table = open(os.path.join(out_dir, "targets.tsv"), "w", encoding="utf-8")
    table.write("\t".join(HEADER) + "\n")
    kept = 0
    started = set()     # output files already holding a record

    for target in targets:
        assembly, locus = target["assembly"], target["locus"]
        sequence = dna[assembly].get(locus)
        residues = protein[assembly].get(locus)

        if not sequence or not residues:
            print(f"WARNING: {assembly} {locus} is not in the annotation - "
                  f"skipped", file=sys.stderr)
            continue
        if len(residues) < min_aa:
            continue

        contig, start, stop, strand, product = where[assembly].get(
            locus, ("-", "-", "-", "-", "-"))

        # The header carries the evidence, because a bare sequence in a folder
        # two weeks later is not a result you can defend.
        header = (f"{target['target']} locus={locus} assembly={assembly}"
                  f" reference={target['reference']} pident={target['pident']}"
                  f" ref_range={target['low']}-{target['high']}"
                  f" location={contig}:{start}-{stop}({strand})")

        for folder, extension, text in (("genes", "fna", sequence),
                                        ("proteins", "faa", residues)):
            path = os.path.join(out_dir, folder,
                                f"{target['reference']}.{extension}")
            with open(path, "a", encoding="utf-8") as output:
                # A blank line between records, so a file holding a dozen
                # pieces of one gene can be read by eye. A program reading
                # FASTA skips the gap: a blank line is not a header and has
                # no residues on it.
                if path in started:
                    output.write("\n")
                output.write(fasta.record(header, text))
            started.add(path)

        table.write("\t".join([
            target["target"], assembly, locus, target["reference"],
            target["pident"], target["aln_length"], str(target["low"]),
            str(target["high"]), target["pct_ref_covered"],
            contig, start, stop, strand, str(len(residues)), product,
        ]) + "\n")
        kept += 1

    table.close()

    print(f"{total} matches from program A, {len(targets)} kept by the "
          f"REFERENCES setting, {kept} targets of {min_aa} residues or longer")

    # How many genes landed under each reference. Several pieces of one
    # reference gene is the normal result here, not a surprise: see the note
    # about frameshifts above.
    per_reference = {}
    for target in targets:
        per_reference[target["reference"]] = \
            per_reference.get(target["reference"], 0) + 1
    for reference in sorted(per_reference):
        print(f"  {per_reference[reference]:>3}  {reference}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
