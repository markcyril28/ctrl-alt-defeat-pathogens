#!/usr/bin/env python3
"""Join the AMR/virulence screening hits to the annotated genes they sit in.

Programs 09 and 10 of the metagenomics pipeline looked at the same assemblies
and wrote coordinates in two different name spaces:

    Bakta (09)   renamed every contig to contig_1, contig_2, ... in input
                 order, and dropped contigs below --min-contig-length
    AMRFinderPlus / ABRicate (10)   kept the assembler's own contig names

So "contig_3" in the annotation and "contig_3" in the screening table are not
the same piece of DNA, and in this dataset they genuinely are not: Bakta's
contig_3 is MetaFlye's contig_11. Joining the two tables by contig name
produces hits attached to the wrong genes, silently.

This script joins them by the one thing that cannot be renamed — the contig
sequence itself. Each Bakta contig is matched to its assembly contig by MD5 of
the uppercased sequence, and the screening hit's coordinates are then compared
against the CDS coordinates on that contig. The CDS with the largest overlap
wins.

    map_screening_hits.py <metagenomics_results_dir>

Prints one TSV row per screening hit to stdout. A hit with no annotated CDS
is reported as a row with locus_tag "-" and a reason, not dropped: a hit on a
contig too short to annotate is a result about your assembly, not a gap to
pass over.
"""
import csv
import glob
import hashlib
import os
import sys

NA = "-"

HEADER = [
    "assembly", "source", "hit_gene", "hit_product", "pct_identity",
    "pct_coverage", "element", "contig", "hit_start", "hit_stop",
    "locus_tag", "bakta_gene", "bakta_product", "aa_length", "overlap_bp",
    "note",
]


def _read_fasta(path):
    """Return {sequence_id: sequence} for a FASTA file."""
    sequences = {}
    name, chunks = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def _fingerprint(sequence):
    return hashlib.md5(sequence.upper().encode()).hexdigest()


def _contig_map(assembly_fasta, bakta_fna):
    """Map assembler contig name -> Bakta contig name, by sequence.

    A sequence appearing twice in one file is left out of the map rather than
    guessed at: an ambiguous join is worse than a missing one.
    """
    def by_fingerprint(path):
        counts, first = {}, {}
        for name, sequence in _read_fasta(path).items():
            key = _fingerprint(sequence)
            counts[key] = counts.get(key, 0) + 1
            first.setdefault(key, name)
        return {k: first[k] for k, v in counts.items() if v == 1}

    assembly = by_fingerprint(assembly_fasta)
    bakta = by_fingerprint(bakta_fna)
    return {name: bakta[key] for key, name in assembly.items() if key in bakta}


def _bakta_cds(bakta_tsv):
    """Return {contig: [(start, stop, locus_tag, gene, product), ...]}."""
    cds = {}
    with open(bakta_tsv, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[1] != "cds":
                continue
            cds.setdefault(f[0], []).append(
                (int(f[2]), int(f[3]), f[5], f[6] or NA, f[7] or NA)
            )
    return cds


def _aa_lengths(bakta_faa):
    lengths = {}
    locus = None
    with open(bakta_faa, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                locus = line[1:].split()[0]
                lengths[locus] = 0
            elif locus:
                lengths[locus] += len(line.strip())
    return lengths


def _best_overlap(cds_list, start, stop):
    """The CDS sharing the most bases with [start, stop], or None."""
    best, best_bp = None, 0
    for cds in cds_list:
        overlap = min(cds[1], stop) - max(cds[0], start) + 1
        if overlap > best_bp:
            best, best_bp = cds, overlap
    return best, best_bp


# ── the two screening formats ───────────────────────────────────────
# Column positions are 0-based and taken from the files programs 10 writes.
# Both tools report start < stop regardless of strand, which is what the
# overlap test above assumes.
#
# A data row is recognised by having integers where the coordinates belong,
# not by matching the header text. AMRFinderPlus renames its columns between
# versions — "Gene symbol" became "Element symbol", "Protein identifier"
# became "Protein id" — while keeping the column order, so a parser that
# matched the header stopped recognising the header and tried to read it as a
# hit. The coordinates are the stable part.

def _is_data_row(row, start_index, stop_index, minimum_columns):
    if len(row) < minimum_columns:
        return False
    try:
        int(row[start_index])
        int(row[stop_index])
    except (ValueError, IndexError):
        return False
    return True


def _amrfinder_hits(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not _is_data_row(row, 2, 3, 17):
                continue
            element = row[8] if len(row) > 8 else NA
            subtype = row[9] if len(row) > 9 else ""
            yield {
                "contig": row[1], "start": row[2], "stop": row[3],
                "gene": row[5], "product": row[6],
                "coverage": row[15], "identity": row[16],
                "element": f"{element}/{subtype}" if subtype else element,
            }


def _abricate_hits(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not _is_data_row(row, 2, 3, 14):
                continue
            yield {
                "contig": row[1], "start": row[2], "stop": row[3],
                "gene": row[5], "product": row[13],
                "coverage": row[9], "identity": row[10],
                "element": row[11],
            }


def _screens(screening_dir, assembly):
    """(source label, path, parser) for every screening table of an assembly."""
    found = []
    amrfinder = os.path.join(screening_dir, f"{assembly}.amrfinder.tsv")
    if os.path.exists(amrfinder):
        found.append(("amrfinder", amrfinder, _amrfinder_hits))
    for path in sorted(glob.glob(os.path.join(screening_dir, f"{assembly}.abricate.*.tsv"))):
        database = os.path.basename(path).split(".")[-2]
        found.append((f"abricate_{database}", path, _abricate_hits))
    return found


def main(argv):
    if len(argv) != 1:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    results = argv[0]
    annotation_dir = os.path.join(results, "09_Annotation")
    screening_dir = os.path.join(results, "10_AMR_Virulence_Screening")
    assemblies_dir = os.path.join(results, "Assemblies")

    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(HEADER)

    for path in sorted(glob.glob(os.path.join(annotation_dir, "*", "*.tsv"))):
        assembly = os.path.basename(os.path.dirname(path))
        if os.path.basename(path) != f"{assembly}.tsv":
            continue            # skip metaspades.hypotheticals.tsv and friends

        assembly_fasta = os.path.join(assemblies_dir, f"{assembly}.fasta")
        bakta_fna = os.path.join(annotation_dir, assembly, f"{assembly}.fna")
        bakta_faa = os.path.join(annotation_dir, assembly, f"{assembly}.faa")
        if not (os.path.exists(assembly_fasta) and os.path.exists(bakta_fna)):
            print(f"skipping {assembly}: no assembly/annotation FASTA pair",
                  file=sys.stderr)
            continue

        name_map = _contig_map(assembly_fasta, bakta_fna)
        cds_by_contig = _bakta_cds(path)
        aa_length = _aa_lengths(bakta_faa) if os.path.exists(bakta_faa) else {}

        screens = _screens(screening_dir, assembly)
        if not screens:
            print(f"skipping {assembly}: no screening tables in {screening_dir}",
                  file=sys.stderr)
            continue

        for source, screen_path, parse in screens:
            for hit in parse(screen_path):
                start, stop = sorted((int(hit["start"]), int(hit["stop"])))
                bakta_contig = name_map.get(hit["contig"])
                row = [
                    assembly, source, hit["gene"] or NA, hit["product"] or NA,
                    hit["identity"], hit["coverage"], hit["element"] or NA,
                    hit["contig"], start, stop,
                ]

                if bakta_contig is None:
                    # Below --min-contig-length in program 09, or a contig
                    # whose sequence appears twice. Either way: no gene call.
                    writer.writerow(row + [NA, NA, NA, NA, 0, "contig_not_annotated"])
                    continue

                cds, overlap = _best_overlap(cds_by_contig.get(bakta_contig, []), start, stop)
                if cds is None:
                    writer.writerow(row + [NA, NA, NA, NA, 0, "no_cds_overlap"])
                    continue

                writer.writerow(row + [
                    cds[2], cds[3], cds[4], aa_length.get(cds[2], 0), overlap, "ok",
                ])

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
