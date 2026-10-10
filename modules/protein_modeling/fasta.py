"""FASTA as every program here writes it: 80 characters to a sequence line.

    import fasta
    handle.write(fasta.record("FimH_HJBCCN_009237 locus=...", sequence))

Only the sequence is wrapped. The header stays on one line however long it
gets, because a header broken across lines is read as sequence.

Bakta writes each sequence on a single line, thousands of characters long for
a big gene. Every tool reads that, but a person cannot: it runs off the side
of the terminal. 80 fits a terminal, and it is the width of the NCBI genomes
in WORKING_FOLDER/INPUT_DATASETS/. Wrapping changes nothing for a program
reading the file: the lines are joined back together.

Where a file holds more than one record, the programs here put a blank line
between them and none at the end — also for reading it, since a run of
records with long headers otherwise runs together. record() returns one
record without that blank line: the writer adds it, because only the writer
knows whether a record is the first in its file. A blank line is neither a
header nor residues, so a program reading FASTA skips it; blastn and
makeblastdb (2.17.0) give identical results either way.
"""

WIDTH = 80


def wrap(sequence):
    """The sequence cut into lines of WIDTH characters, the last one shorter."""
    return [sequence[at:at + WIDTH] for at in range(0, len(sequence), WIDTH)]


def record(header, sequence):
    """One record, ready to write: ">" + header, then the wrapped sequence."""
    return "\n".join([f">{header}"] + wrap(sequence)) + "\n"


def records(lines):
    """(header, sequence) for each record in some lines of FASTA.

    The header comes back without its ">", and the sequence joined into one
    string however it was wrapped. Anything before the first header belongs
    to no record and is skipped.
    """
    header, pieces = None, []
    for line in lines:
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(pieces)
            header, pieces = line[1:], []
        elif header is not None:
            pieces.append(line.strip())
    if header is not None:
        yield header, "".join(pieces)
