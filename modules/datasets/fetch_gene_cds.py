#!/usr/bin/env python3
"""Fetch the reference gene CDS files — the GENE_*.gene.fna files under
WORKING_FOLDER/INPUT_DATASETS/from_Database/.

This script is the only thing that should write those files. The coordinates
below are the record of where each gene is; the sequences are downloaded from
NCBI so no sequence is ever edited by hand.

COORDINATE CONVENTION — 1-based, both ends included
---------------------------------------------------
``start`` and ``stop`` are 1-based positions in the accession, and both are
part of the gene: a 531 nt gene has ``abs(stop - start) + 1 == 531``. Position
1 is the first base of the record. ``start > stop`` means the gene is on the
minus strand and is read from ``start`` down to ``stop``.

This is the convention NCBI uses in a GenBank feature table, in the efetch
``seq_start``/``seq_stop`` parameters this script sends, and in the FASTA
headers it writes. It is NOT the convention of a BED file, of a GFF ``start``
read as an array index, or of a Python slice — all of which count the first
base as 0. Taking a 0-based start and sending it to efetch unchanged fetches a
window that begins one base too early: the first codon is then the last base
before the gene plus the first two bases of ATG, every codon after it is read
one base out of frame, and the translation is garbage even though the file has
exactly the right number of bases.

That is a real mistake, not a hypothetical: six of the seven files in the
dataset were first extracted that way. A shifted window is easy to miss
because the length still looks right, so verify() below checks the things a
shift actually breaks — the start codon, the stop codon, and the absence of
internal stops — and nothing is written unless every gene passes.

Usage (from the repository root):

    conda activate protein_modeling
    python modules/datasets/fetch_gene_cds.py            # download and verify
    python modules/datasets/fetch_gene_cds.py --offline   # audit the files on disk
    python modules/datasets/fetch_gene_cds.py --write      # rewrite the files

--write keeps the previous file as <name>.orig if no .orig is there yet,
because WORKING_FOLDER/INPUT_DATASETS/ is in .gitignore and an overwrite
cannot be undone with git.

Standard library only - no extra conda packages required.
"""

import argparse
import hashlib
import os
import shutil
import sys
import time
import urllib.parse
import urllib.request

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL = "ologist-workshop-docs"
USER_AGENT = TOOL + "/1.0"

# One entry per gene file.
#   path        - relative to WORKING_FOLDER/INPUT_DATASETS/from_Database/
#   header      - FASTA header, without the trailing coordinates
#   accession   - RefSeq nucleotide record the coordinates belong to
#   start, stop - 1-based and inclusive, see the convention above
#   nt          - expected length, checked against the coordinates and the download
#   aa          - expected protein length, so nt is not just self-consistent
#   protein     - sibling FASTA_*.fasta to compare the translation against, or
#                 None where the dataset's protein is from a different strain
GENES = [
    {
        "path": "1_Mycobacterium_tuberculosis/GENE_IFT_immunity.gene.fna",
        "header": ">IFT_immunity|Rv3902c|Mycobacterium tuberculosis H37Rv",
        "accession": "NC_000962.3", "start": 4387895, "stop": 4387365,
        "nt": 531, "aa": 176,
        "protein": "1_Mycobacterium_tuberculosis/FASTA_O05443_IFT_immunity.fasta",
    },
    {
        "path": "1_Mycobacterium_tuberculosis/GENE_TNT_CpnT.gene.fna",
        "header": ">TNT_CpnT|cpnT/Rv3903c|Mycobacterium tuberculosis H37Rv",
        "accession": "NC_000962.3", "start": 4390432, "stop": 4387892,
        "nt": 2541, "aa": 846,
        "protein": "1_Mycobacterium_tuberculosis/FASTA_O05442_CpnT.fasta",
    },
    {
        "path": "3_Klebsiella_pneumoniae/GENE_FimH.gene.fna",
        "header": ">FimH|fimH|Klebsiella pneumoniae MGH 78578",
        "accession": "NC_012731.1", "start": 4386276, "stop": 4387184,
        "nt": 909, "aa": 302,
        # A0A0H3H2I8 is ST11, this gene is MGH 78578: same length, not the same strain.
        "protein": None,
    },
    {
        "path": "3_Klebsiella_pneumoniae/GENE_MrkD.gene.fna",
        "header": ">MrkD|mrkD|Klebsiella pneumoniae MGH 78578",
        "accession": "NC_012731.1", "start": 4369882, "stop": 4368887,
        "nt": 996, "aa": 331,
        # P21648 is 321 aa, from another strain — see
        # WORKING_FOLDER/INPUT_DATASETS/README.md.
        "protein": None,
    },
    {
        "path": "4_Pseudomonas_aeruginosa/GENE_ExoA_toxA.gene.fna",
        "header": ">ExoA_toxA|toxA|Pseudomonas aeruginosa PAO1",
        "accession": "NC_002516.2", "start": 1242500, "stop": 1240584,
        "nt": 1917, "aa": 638,
        "protein": "4_Pseudomonas_aeruginosa/FASTA_P11439_ExoA.fasta",
    },
    {
        "path": "4_Pseudomonas_aeruginosa/GENE_ExoS.gene.fna",
        "header": ">ExoS|exoS|Pseudomonas aeruginosa PAO1",
        "accession": "NC_002516.2", "start": 4304502, "stop": 4303141,
        "nt": 1362, "aa": 453,
        "protein": "4_Pseudomonas_aeruginosa/FASTA_G3XDA1_ExoS_PAE.fasta",
    },
    {
        "path": "4_Pseudomonas_aeruginosa/GENE_ExoT.gene.fna",
        "header": ">ExoT|exoT|Pseudomonas aeruginosa PAO1",
        "accession": "NC_002516.2", "start": 58786, "stop": 60159,
        "nt": 1374, "aa": 457,
        "protein": "4_Pseudomonas_aeruginosa/FASTA_Q9I788_ExoT_PAE.fasta",
    },
]

# Bacterial genetic code (NCBI translation table 11). Table 11 differs from the
# standard code only in which codons may start a gene, which START_CODONS holds,
# so the codon-to-amino-acid mapping below is the standard one.
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    first + second + third: AMINO[i * 16 + j * 4 + k]
    for i, first in enumerate(BASES)
    for j, second in enumerate(BASES)
    for k, third in enumerate(BASES)
}
START_CODONS = ("ATG", "GTG", "TTG")
STOP_CODONS = ("TAA", "TAG", "TGA")


def coordinates(gene):
    """The gene's coordinates as they are written in the FASTA header."""
    return f"{gene['accession']}:{gene['start']}-{gene['stop']}"


def span(start, stop):
    """Bases in a 1-based window with both ends included, either strand."""
    return abs(stop - start) + 1


def fetch(gene):
    """Download one window from NCBI. start > stop means the minus strand."""
    minus = gene["start"] > gene["stop"]
    low, high = sorted((gene["start"], gene["stop"]))
    query = urllib.parse.urlencode({
        "db": "nuccore", "id": gene["accession"],
        "rettype": "fasta", "retmode": "text",
        # efetch is 1-based and inclusive, the same as the table above, so the
        # coordinates go out as they are. Lower bound first either way; the
        # strand, not the order, is what reverse-complements the result.
        "seq_start": str(low), "seq_stop": str(high),
        "strand": "2" if minus else "1", "tool": TOOL,
    })
    request = urllib.request.Request(EFETCH + "?" + query,
                                     headers={"User-Agent": USER_AGENT})
    text = urllib.request.urlopen(request, timeout=90).read().decode()
    return sequence_of(text.splitlines())


def sequence_of(lines):
    """The sequence of a one-record FASTA, upper-cased, as a single string."""
    return "".join(line.strip() for line in lines
                   if not line.startswith(">")).upper()


def translate(dna):
    """Translate a CDS, reading its first codon as methionine."""
    protein = [CODON_TABLE.get(dna[i:i + 3], "X")
               for i in range(0, len(dna) - len(dna) % 3, 3)]
    if protein and dna[:3] in START_CODONS:
        protein[0] = "M"
    return "".join(protein)


def verify(gene, dna):
    """Problems with a sequence for this gene; an empty list means it is good.

    The first four checks are what a one-base shift breaks.
    """
    problems = []
    if dna[:3] not in START_CODONS:
        problems.append(f"begins {dna[:3]}, not a start codon "
                        f"({'/'.join(START_CODONS)})")
    if dna[-3:] not in STOP_CODONS:
        problems.append(f"ends {dna[-3:]}, not a stop codon "
                        f"({'/'.join(STOP_CODONS)})")
    if len(dna) % 3:
        problems.append(f"{len(dna)} nt is not a whole number of codons")
    protein = translate(dna)[:-1]            # drop the stop
    if "*" in protein:
        problems.append(f"{protein.count('*')} internal stop codon(s), "
                        f"first at codon {protein.index('*') + 1}")
    if len(dna) != gene["nt"]:
        problems.append(f"{len(dna)} nt, expected {gene['nt']}")
    if len(protein) != gene["aa"]:
        problems.append(f"translates to {len(protein)} aa, "
                        f"expected {gene['aa']}")
    return problems


def check_table(gene):
    """Problems with the entry itself, before anything is downloaded."""
    problems = []
    counted = span(gene["start"], gene["stop"])
    if counted != gene["nt"]:
        problems.append(f"{coordinates(gene)} spans {counted} nt, "
                        f"but nt says {gene['nt']} — an off-by-one in the table")
    if gene["nt"] % 3:
        problems.append(f"nt {gene['nt']} is not a whole number of codons")
    elif gene["nt"] // 3 - 1 != gene["aa"]:
        problems.append(f"nt {gene['nt']} is {gene['nt'] // 3 - 1} codons "
                        f"before the stop, but aa says {gene['aa']}")
    return problems


def identity(one, two):
    """Percent identity of two equal-length strings, or None if they differ."""
    if len(one) != len(two) or not one:
        return None
    same = sum(1 for a, b in zip(one, two) if a == b)
    return 100.0 * same / len(one)


def against_protein(gene, dna, root):
    """One informational line comparing the translation to the dataset protein."""
    if not gene["protein"]:
        return None
    path = os.path.join(root, gene["protein"])
    if not os.path.isfile(path):
        return f"  (no {gene['protein']} to compare against)"
    with open(path, encoding="utf-8") as handle:
        expected = sequence_of(handle.read().splitlines())
    percent = identity(translate(dna)[:-1], expected)
    name = os.path.basename(gene["protein"])
    if percent is None:
        return f"  translation is a different length from {name}"
    flag = "" if percent >= 95.0 else "   <-- check this"
    return f"  {percent:.1f}% identical to {name}{flag}"


def report(gene, dna, root):
    """Print one gene's verdict; return True if it verified."""
    problems = verify(gene, dna)
    label = os.path.basename(gene["path"])
    if problems:
        print(f"FAILED  {label:28s} {len(dna):>5} nt  {coordinates(gene)}")
        for problem in problems:
            print(f"  {problem}")
        return False
    print(f"ok      {label:28s} {len(dna):>5} nt  {coordinates(gene)}  "
          f"{dna[:3]}...{dna[-3:]}")
    note = against_protein(gene, dna, root)
    if note:
        print(note)
    return True


def offline(root):
    """Audit the files on disk, without touching the network.

    Checks each file against this table: same coordinates in the header, same
    sequence length, and a sequence that reads as a CDS. This is what catches a
    file that was written by some other route than this script.
    """
    good = True
    for gene in GENES:
        path = os.path.join(root, gene["path"])
        label = os.path.basename(gene["path"])
        if not os.path.isfile(path):
            print(f"MISSING {label:28s} {path}")
            good = False
            continue
        with open(path, encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        on_disk = lines[0].split("|")[-1].strip() if lines else ""
        if on_disk != coordinates(gene):
            print(f"FAILED  {label:28s} header says {on_disk}, "
                  f"table says {coordinates(gene)}")
            good = False
            continue
        if not report(gene, sequence_of(lines), root):
            good = False
    return good


def write(gene, dna, root):
    """Write one gene file, keeping any previous version as <name>.orig."""
    path = os.path.join(root, gene["path"])
    backup = path + ".orig"
    if os.path.isfile(path) and not os.path.isfile(backup):
        shutil.copy2(path, backup)
        print(f"kept    {os.path.basename(backup)}")
    header = f"{gene['header']}|{coordinates(gene)}"
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(header + "\n" + dna + "\n")
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]
    print(f"wrote   {path}  sha256 {digest}")
    return digest


def main():
    parser = argparse.ArgumentParser(
        description="Fetch and verify the GENE_*.gene.fna reference CDS files.")
    parser.add_argument("--write", action="store_true",
                        help="rewrite the gene files (default: verify only)")
    parser.add_argument("--offline", action="store_true",
                        help="audit the files on disk, no download")
    parser.add_argument("--root", default="WORKING_FOLDER/INPUT_DATASETS/from_Database",
                        help="directory holding the per-species folders")
    args = parser.parse_args()

    # The table before the network: coordinates that do not span the length
    # they claim are the bug this script exists to prevent, and no download
    # can tell you about them.
    broken = False
    for gene in GENES:
        for problem in check_table(gene):
            print(f"TABLE   {os.path.basename(gene['path'])}: {problem}")
            broken = True
    if broken:
        print("\nFix the coordinates in GENES first; nothing was fetched.")
        return 1

    if args.offline:
        if offline(args.root):
            print(f"\nAll {len(GENES)} files on disk match the coordinates "
                  f"in this script and read as a CDS.")
            return 0
        print("\nThe files on disk do not match this script. Re-run with "
              "--write to replace them from NCBI.")
        return 1

    fetched, good = {}, True
    for gene in GENES:
        dna = fetch(gene)
        fetched[gene["path"]] = dna
        good = report(gene, dna, args.root) and good
        time.sleep(0.5)             # NCBI asks for no more than 3 requests/second

    if not good:
        print("\nOne or more sequences did not verify; nothing was written.")
        return 1
    if not args.write:
        print(f"\nAll {len(GENES)} sequences verified. "
              f"Re-run with --write to replace the files on disk.")
        return 0

    print()
    digests = [(gene["path"], write(gene, fetched[gene["path"]], args.root))
               for gene in GENES]
    print("\nUpdate the sha256_file column in WORKING_FOLDER/INPUT_DATASETS/MANIFEST.tsv:")
    for path, digest in digests:
        print(f"  {path}\t{digest}")
    print("The coordinates in WORKING_FOLDER/INPUT_DATASETS/README.md name these files too.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
