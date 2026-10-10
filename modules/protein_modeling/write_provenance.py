#!/usr/bin/env python3
"""Append the provenance of a modelling run to the version record.

    write_provenance.py <results_dir> <upstream_versions.txt>

Prints three sections to stdout, for the calling script to append to
versions.txt:

    the reference genes program A searched
    blanks for the two web services, to fill in by hand
    the metagenomics pipeline's own version record

Which reference gene files were in the folder on the day decides which genes
could be found at all, so the list is recorded rather than described.

The services are not installed here, so their versions cannot be read. They
are left as blanks to fill in rather than omitted, because a structure
prediction depends on the service's model version and the day it ran —
neither of which is in the file it returns, and the AlphaFold Server is not
the AlphaFold3 you install. A model without them cannot be reproduced or
compared to anyone else's, including your own from last month.

The upstream record is copied in rather than referenced because it is the
reason these targets exist at all: a particular Bakta database called these
open reading frames, and a particular BLAST run matched them to a reference
gene.
"""
import os
import sys

BLANKS = """===== web services — fill these in by hand
SWISS-MODEL      date submitted: ____________  template PDB/chain: ____________
                 sequence identity: ______  coverage: ______  QMEANDisCo: ______
AlphaFold Server date submitted: ____________  model version reported: ____________
AlphaFold3 local version: ____________  weights obtained: ____________
                 database snapshot date: ____________"""


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    results, upstream = argv
    genes = os.path.join(results, "A_Gene_Catalog", "reference_genes.tsv")

    print("===== reference genes searched (program A)")
    if os.path.isfile(genes):
        with open(genes, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                field = line.rstrip("\n").split("\t")
                print("\t".join(field[:2]))
    else:
        print(f"({genes} not found — run: bash A_gene_catalog_blastn.sh)")

    print(BLANKS)

    print("===== upstream (metagenomics pipeline, program 11)")
    if os.path.isfile(upstream):
        with open(upstream, encoding="utf-8", errors="replace") as handle:
            print(handle.read().rstrip("\n"))
    else:
        print(f"({upstream} not found — run: bash 11_report_python.sh)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
