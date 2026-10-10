#!/usr/bin/env python3
"""Build program A's BLAST inputs, and turn its output into tables.

    gene_catalog.py references <reference_dir> <out_dir>
    gene_catalog.py queries    <annotation_dir> <out_dir>
    gene_catalog.py tables     <out_dir> <min_identity> <min_length>
    gene_catalog.py matched    <out_dir>

Four steps, called in that order, with makeblastdb and blastn run between
them by the calling script.

references  Pools every WORKING_FOLDER/INPUT_DATASETS/from_Database/*/GENE_*.gene.fna
            into one FASTA, renaming each sequence to its folder and gene,
            "3_Klebsiella_pneumoniae__FimH", and keeping the original header
            in reference_genes.tsv. The renaming matters: downloaded headers
            are long and full of spaces and pipes, and BLAST cuts a sequence
            name at the first space, so two genes whose headers begin alike
            become one name in the results.

queries     Pools the DNA of every hypothetical CDS from every assembly.
            Bakta lists the CDS it could not name in .hypotheticals.faa, but
            that file holds proteins and blastn compares nucleotides — so the
            .hypotheticals.faa is read for its locus tags and the .ffn for
            the sequences. Also records which annotation this was built from,
            in source.tsv, because Bakta assigns locus tags per run.

tables      hits.tsv, matches.tsv, reference_coverage.tsv and
            search_summary.tsv. hits.tsv is the evidence and matches.tsv is
            the opinion: nothing is filtered out of hits.tsv beyond the
            e-value blastn was given, so the thresholds can be changed and
            argued about without running the search again. matches.tsv is
            grouped by bacterial species and ranked inside a species by
            e-value, then percent identity; hits.tsv keeps bitscore order.

matched     blast/matched_<species>.fna, the DNA of the matched loci,
            split into one query file per bacterial species, for a second
            blastn run per species that writes the pairwise alignments.
            Pulling the subset out first is the point: asking blastn to draw
            every alignment over the whole query file would bury the handful
            that matched under a "No hits found" block for each of the
            thousands that did not. The split is on the query side, and each
            run still searches the whole reference database, so an e-value in
            an alignment is the same number the tables report. Searching a
            species against its own genes alone would shrink the search space
            and quietly restate every e-value.

Two things in the tables that are easy to get wrong by hand:

    A minus-strand hit runs backwards along the reference, so ref_start is
    larger than ref_end. The span is the distance between them either way.

    Coverage is counted base by base, as a union. One gene broken into three
    overlapping fragments by frameshifts does not cover 150% of a reference,
    and a column that says it does invites exactly the wrong conclusion.
"""
import glob
import os
import re
import sys

import fasta

# evalue and pct_identity sit fourth and fifth, beside the reference they
# judge: they are what the table is sorted on and what gets read first. The
# rows are built, sorted and read by position, so this order is the one the
# row[...] indexes in tables() and in extract_targets.py follow.
HITS_HEADER = [
    "assembly", "locus_tag", "reference", "evalue", "pct_identity",
    "aln_length", "locus_start", "locus_end", "ref_start", "ref_end",
    "locus_length", "ref_length", "pct_ref_covered", "bitscore",
]


def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(cell) for cell in row) + "\n")


def lines_of(path):
    """Every line of a text file, without its line ending.

    Stripped of "\\r" as well: the reference files were saved on Windows.
    """
    with open(path, encoding="utf-8", errors="replace") as handle:
        return [line.rstrip("\r\n") for line in handle]


# ── references ──────────────────────────────────────────────────────

def references(reference_dir, out_dir):
    pool = os.path.join(out_dir, "reference_genes.fna")
    rows, organisms, unsearched = [], set(), []
    written = 0         # records in the pool, so far

    with open(pool, "w", encoding="utf-8") as output:
        for folder in sorted(glob.glob(os.path.join(reference_dir, "*"))):
            if not os.path.isdir(folder):
                continue
            organism = os.path.basename(folder)
            files = sorted(glob.glob(os.path.join(folder, "GENE_*.gene.fna")))

            # An organism folder with no gene FASTA is not searched. Said out
            # loud, so an empty result for it is not mistaken for an absence.
            if not files:
                unsearched.append(organism)
                continue
            organisms.add(organism)

            for path in files:
                for header, sequence in fasta.records(lines_of(path)):
                    name = re.sub(r"[^A-Za-z0-9._-]", "_",
                                  f"{organism}__{header.split('|')[0]}")
                    sequence = re.sub(r"[ \t]", "", sequence)
                    if written:                 # blank line between records
                        output.write("\n")
                    output.write(fasta.record(name, sequence))
                    written += 1
                    rows.append([name, len(sequence), path, header])

    write_tsv(os.path.join(out_dir, "reference_genes.tsv"),
              ["reference", "length", "file", "header"], rows)

    if not rows:
        print(f"ERROR: no reference genes (GENE_*.gene.fna) under "
              f"{reference_dir}", file=sys.stderr)
        return 1

    print(f"{len(rows)} reference genes from {len(organisms)} organism(s), "
          f"listed in {out_dir}/reference_genes.tsv")
    for organism in unsearched:
        print(f"  {organism} — no GENE_*.gene.fna, not searched")
    return 0


# ── queries ─────────────────────────────────────────────────────────

def queries(annotation_dir, out_dir):
    source, searched = [], 0
    pool = open(os.path.join(out_dir, "hypotheticals.fna"), "w",
                encoding="utf-8")
    proteins = open(os.path.join(out_dir, "proteins.faa"), "w",
                    encoding="utf-8")
    in_proteins, in_pool = 0, 0     # records in each pool, so far

    for folder in sorted(glob.glob(os.path.join(annotation_dir, "*"))):
        if not os.path.isdir(folder):
            continue
        name = os.path.basename(folder)
        stem = os.path.join(folder, name)

        if not (os.path.isfile(f"{stem}.tsv") and os.path.isfile(f"{stem}.faa")):
            print(f"skipping {name} — no Bakta .tsv/.faa pair in {folder}")
            continue

        source.append([name, f"{stem}.tsv", os.path.getsize(f"{stem}.tsv"),
                       int(os.path.getmtime(f"{stem}.tsv"))])

        # Every protein Bakta translated, pooled into one file. No later
        # program reads this — program B goes back to the per-assembly .faa —
        # so it is here to be searched by hand: one grep to find what a locus
        # tag translates to without opening the annotation. The assembly name
        # goes into the ID because two assemblies annotated in separate runs
        # can carry the same locus tag prefix.
        for header, sequence in fasta.records(lines_of(f"{stem}.faa")):
            if in_proteins:                     # blank line between records
                proteins.write("\n")
            proteins.write(fasta.record(f"{name}__{header}", sequence))
            in_proteins += 1

        hypotheticals = f"{stem}.hypotheticals.faa"
        if (not os.path.isfile(hypotheticals)
                or os.path.getsize(hypotheticals) == 0
                or not os.path.isfile(f"{stem}.ffn")):
            print(f"skipping {name} — no {name}.hypotheticals.faa or "
                  f"{name}.ffn to search")
            continue

        wanted = {line[1:].split()[0] for line in lines_of(hypotheticals)
                  if line.startswith(">")}
        found = 0
        for header, sequence in fasta.records(lines_of(f"{stem}.ffn")):
            tag = header.split()[0]
            if tag in wanted:
                found += 1
                if in_pool:                     # blank line between records
                    pool.write("\n")
                pool.write(fasta.record(f"{name}__{tag}", sequence))
                in_pool += 1

        searched += found
        # A tag on the hypothetical list with no sequence in the .ffn is a
        # gene that will not be searched, and an unsearched gene looks exactly
        # like a gene with no match.
        if found != len(wanted):
            print(f"WARNING: {name} lists {len(wanted)} hypothetical CDS but "
                  f"its .ffn holds {found} of them", file=sys.stderr)

    pool.close()
    proteins.close()
    write_tsv(os.path.join(out_dir, "source.tsv"),
              ["assembly", "annotation_tsv", "bytes", "mtime"], source)

    if not searched:
        print(f"ERROR: no hypothetical CDS found under {annotation_dir}",
              file=sys.stderr)
        return 1

    print(f"{searched} hypothetical CDS to search, from "
          f"{len(source)} assembly/assemblies")
    return 0


# ── tables ──────────────────────────────────────────────────────────

def tables(out_dir, min_identity, min_length):
    min_identity, min_length = float(min_identity), float(min_length)

    rows = []
    for line in lines_of(os.path.join(out_dir, "blast", "blastn.out")):
        field = line.split("\t")
        if len(field) < 12:
            continue

        # The query name is "<assembly>__<locus_tag>", joined by queries().
        assembly, locus = field[0].split("__", 1)
        ref_start, ref_end = int(field[6]), int(field[7])
        ref_length = int(field[9])
        span = abs(ref_end - ref_start) + 1

        rows.append([assembly, locus, field[1], field[10], field[2],
                     field[3], field[4], field[5], ref_start, ref_end,
                     field[8], ref_length,
                     f"{100.0 * span / ref_length:.1f}", field[11]])

    # Strongest first, by bitscore, so the first row of a locus is its best.
    rows.sort(key=lambda row: float(row[13]), reverse=True)
    write_tsv(os.path.join(out_dir, "hits.tsv"), HITS_HEADER, rows)

    matches, seen = [], set()
    for row in rows:
        locus_key = (row[0], row[1])
        if locus_key in seen:
            continue
        if float(row[4]) >= min_identity and float(row[5]) >= min_length:
            seen.add(locus_key)
            matches.append(row)

    # Read as a table, matches.tsv is wanted grouped by species rather than
    # by bitscore, so the ordering is applied here, after the dedup that
    # needs hits.tsv's strongest-first order. A species is the part of a
    # reference name before "__", and the 1-4 prefix it carries keeps the
    # four datasets in their own order rather than an alphabetical one.
    matches.sort(key=lambda row: (row[0], row[2].split("__", 1)[0],
                                  float(row[3]), -float(row[4])))
    write_tsv(os.path.join(out_dir, "matches.tsv"), HITS_HEADER, matches)

    covered, loci, lengths, order = {}, {}, {}, []
    for row in matches:
        key = (row[0], row[2])
        if key not in covered:
            covered[key], loci[key] = set(), 0
            order.append(key)
        loci[key] += 1
        lengths[key] = int(row[11])
        low, high = sorted((int(row[8]), int(row[9])))
        covered[key].update(range(low, high + 1))

    write_tsv(os.path.join(out_dir, "reference_coverage.tsv"),
              ["assembly", "reference", "ref_length", "loci",
               "pct_ref_covered"],
              [[key[0], key[1], lengths[key], loci[key],
                f"{100.0 * len(covered[key]) / lengths[key]:.1f}"]
               for key in order])

    searched = {}
    for line in lines_of(os.path.join(out_dir, "hypotheticals.fna")):
        if line.startswith(">"):
            assembly = line[1:].split("__")[0]
            searched[assembly] = searched.get(assembly, 0) + 1

    with_a_hit, matched = {}, {}
    for assembly, locus in {(row[0], row[1]) for row in rows}:
        with_a_hit[assembly] = with_a_hit.get(assembly, 0) + 1
    for row in matches:
        matched[row[0]] = matched.get(row[0], 0) + 1

    write_tsv(os.path.join(out_dir, "search_summary.tsv"),
              ["assembly", "hypotheticals", "with_a_hit", "matched"],
              [[assembly, searched[assembly], with_a_hit.get(assembly, 0),
                matched.get(assembly, 0)] for assembly in searched])

    print(f"{len(rows)} hits, {len(matches)} of them matches at "
          f"identity >= {min_identity:g}% and alignment >= {min_length:g} bp")
    if not matches:
        print("No hypothetical CDS matched a reference gene.")
    return 0


# ── matched ──────────────────────────────────────────────

def matched(out_dir):
    """The DNA of the matched loci, one query file per bacterial species.

    Split into files here, before the search, because one blastn run writes
    one output file. A species is the part of a reference name before "__",
    the same grouping matches.tsv is sorted by, and within a species the loci
    keep their matches.tsv order: blastn works through a query file from the
    top, so the alignments come out ranked the way the table is.

    A locus is filed under the species of the reference it matched. Its
    alignments against any other species' genes still appear in that file,
    because the run is against the whole database — the evidence for a locus
    stays in one place rather than being split across files.
    """
    rows = [line.split("\t")
            for line in lines_of(os.path.join(out_dir, "matches.tsv"))[1:]
            if line.strip()]

    species_of = {}
    for row in rows:
        species = row[2].split("__", 1)[0]
        species_of.setdefault(species, []).append(f"{row[0]}__{row[1]}")

    sequences = {}
    for header, sequence in fasta.records(
            lines_of(os.path.join(out_dir, "hypotheticals.fna"))):
        sequences[header.split()[0]] = sequence

    # A species that matched last time and not this time would otherwise
    # leave its old query file behind for the search to pick up again.
    for stale in glob.glob(os.path.join(out_dir, "blast", "matched_*.fna")):
        os.remove(stale)

    missing, written = [], 0
    for species in sorted(species_of):
        path = os.path.join(out_dir, "blast", f"matched_{species}.fna")
        with open(path, "w", encoding="utf-8") as output:
            in_file = 0
            for name in species_of[species]:
                if name not in sequences:
                    missing.append(name)
                    continue
                if in_file:                     # blank line between records
                    output.write("\n")
                output.write(fasta.record(name, sequences[name]))
                in_file += 1
        written += in_file
        print(f"{in_file:3d} matched loci to align  {path}")

    # hypotheticals.fna is what was searched, so every matched locus has to
    # be in it. One that is not means the two files are from different runs.
    if missing:
        print(f"ERROR: {len(missing)} matched loci are not in "
              f"hypotheticals.fna, starting with {missing[0]} — rerun the "
              f"whole program, the files are from different runs",
              file=sys.stderr)
        return 1

    if not written:
        print("no matches, so no alignments to draw")
        return 0

    print(f"{written} matched loci across {len(species_of)} species")
    return 0


COMMANDS = {"references": references, "queries": queries, "tables": tables,
            "matched": matched}
ARGUMENTS = {"references": 2, "queries": 2, "tables": 3, "matched": 1}


def main(argv):
    if not argv or argv[0] not in COMMANDS:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    command, rest = argv[0], argv[1:]
    if len(rest) != ARGUMENTS[command]:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    return COMMANDS[command](*rest)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
