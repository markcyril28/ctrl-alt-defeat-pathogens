#!/usr/bin/env python3
"""Which reference organism a gene, a target or a model belongs to.

    species_of.py <prep_report.tsv> <target or model name>

Prints the organism folder name, so a shell loop can ask for one without
holding the rule itself:

    SPECIES=$(python3 modules/protein_modeling/species_of.py \\
        "$PREP/prep_report.tsv" "$NAME")

Four programs need this answer and they would otherwise each carry their own
copy of it, which is how two folders of the same organism end up spelled two
ways. There are two rules, used at two different points in the pipeline:

    from a reference gene name   Program A renamed every reference gene to
                                 "<organism>__<gene>" when it pooled them, so
                                 the organism is in the name. This is how
                                 program C groups the sequences in the first
                                 place.

    from a model or receptor     By the time a model comes back from a web
                                 service the name is SWISSMODEL_<target> or
                                 AF3_<target>, and the services mangle it
                                 further — AlphaFold lowercases the job name
                                 and wraps it in fold_..._model_0. So the
                                 target is recovered from the file name and
                                 looked up in program C's report, which is
                                 the table that recorded the grouping.

What the answer does not mean. The organism is the reference gene the
sequence matched. These contigs are a community and nothing in this pipeline
bins them, so a hit to a Klebsiella gene says the sequence is in the sample,
not that Klebsiella carries it.
"""
import os
import sys

UNASSIGNED = "unassigned"


def of_reference(reference):
    """The organism folder a reference gene came from, from its name.

    basename() because this becomes a folder name: the organism comes from a
    folder under WORKING_FOLDER/INPUT_DATASETS/, but a name that is about to
    be joined to a path is worth keeping to one path component.
    """
    organism = os.path.basename(reference.split("__")[0].strip())
    if "__" not in reference or organism in ("", ".", ".."):
        return UNASSIGNED
    return organism


def by_target(prep_report):
    """Program C's grouping, as {target: organism}.

    Read from the report rather than from the folder names so that a target
    C dropped — too short, a duplicate — still has an organism here.
    """
    grouped = {}
    if not os.path.exists(prep_report):
        return grouped
    with open(prep_report, encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            field = line.rstrip("\n").split("\t")
            if index and len(field) >= 2 and field[0]:
                grouped[field[0]] = field[1]
    return grouped


def of_model(name, grouped):
    """The organism of the target a model's file name carries.

    Longest target first, so a target whose name is contained in another one
    cannot claim it. A name no target matches is unassigned rather than an
    error: usually it is a download from a job this pipeline did not create,
    and it is still worth keeping and saying so.
    """
    for target in sorted(grouped, key=len, reverse=True):
        if target.lower() in name.lower():
            return grouped[target]
    return UNASSIGNED


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    report, name = argv
    grouped = by_target(report)

    # A name that is itself a target is answered directly; anything else is a
    # model or receptor name with a target inside it.
    print(grouped.get(name) or of_model(name, grouped))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
