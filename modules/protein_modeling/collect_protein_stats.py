#!/usr/bin/env python3
"""Turn the protein pipeline's scattered reports into the worksheet's tables.

Every program in this pipeline writes its own report in its own format. This
script rearranges them and nothing else: it reads only files the pipeline
already wrote, and it computes no statistics of its own.

    collect_protein_stats.py targets <results_dir>   gene -> sequence -> model
    collect_protein_stats.py models  <results_dir>   what came back, and how good
    collect_protein_stats.py docking <results_dir>   best ligand per receptor
    collect_protein_stats.py selfcheck               run the built-in checks

<results_dir> is WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset.
Each table prints as TSV to stdout.
A stage you have not run yet is reported as a row, not as a crash: the models
arrive from web services on their own schedule, and a half-finished pipeline
is the normal state of this work for days at a time.

Which file each table comes from, so a renamed program folder has one place to
be fixed:

    targets   B_Gene_Protein_Extraction/targets.tsv  + C_Protein_Prep/prep_report.tsv
    models    F_Model_QC/model_qc.tsv                + C_Protein_Prep/prep_report.tsv
    docking   I_Docking/docking_scores.tsv

The species column in all three tables is the reference organism the gene
matched, which is the folder program C grouped the sequence into. It is read
from program C's report, by species_of.py, rather than worked out again here.
It is not a statement about which organism in the sample carries the
sequence: these contigs are a community and nothing here bins them.
"""
import contextlib
import io
import os
import sys

import species_of

NA = "-"

EXTRACTION = ("B_Gene_Protein_Extraction", "targets.tsv")
PREP = ("C_Protein_Prep", "prep_report.tsv")
MODEL_QC = ("F_Model_QC", "model_qc.tsv")
DOCKING = ("I_Docking", "docking_scores.tsv")


def _rows(path, skip_header=True):
    """Yield the non-empty rows of a TSV as lists of strings."""
    if not os.path.exists(path):
        return
    with open(path, encoding="utf-8", errors="replace") as fh:
        for index, line in enumerate(fh):
            if skip_header and index == 0:
                continue
            line = line.rstrip("\n")
            if line.strip():
                yield line.split("\t")


def _table(path, column_count):
    """Rows of a TSV, padded so a short row cannot raise IndexError."""
    for row in _rows(path):
        yield row + [NA] * (column_count - len(row))


def _print(header, rows):
    print("\t".join(header))
    for row in rows:
        print("\t".join(str(cell) for cell in row))


# ── the organism each target was grouped under, from C's report ─────

def _prepared(prep):
    """C's row per target, keyed by target name."""
    return {row[0]: row for row in _table(prep, 13)}


# ── targets: B's extraction joined to C's verdict and F's models ────

def targets(results):
    extraction = os.path.join(results, *EXTRACTION)
    prep = os.path.join(results, *PREP)
    qc = os.path.join(results, *MODEL_QC)

    header = ["target", "species", "assembly", "locus_tag", "reference",
              "pct_identity", "aln_length", "pct_ref_covered", "aa_length",
              "prep_verdict", "models", "product"]

    if not os.path.exists(extraction):
        _print(header, [["(B_Gene_Protein_Extraction not run)"] + [NA] * 11])
        return

    # C's row per target: the organism it was grouped under, and whether the
    # sequence was fit to submit at all.
    prepared = _prepared(prep)

    # How many models came back for each target. The model name carries the
    # target inside it — SWISSMODEL_<target> or AF3_<target> — which is the
    # only link between a download and the gene it came from.
    model_count = {}
    for row in _table(qc, 12):
        for target in prepared:
            if target.lower() in row[0].lower():
                model_count[target] = model_count.get(target, 0) + 1

    rows = []
    for row in _table(extraction, 15):
        target = row[0]
        prep_row = prepared.get(target)
        rows.append([
            target,
            # The organism is in the reference gene's name either way, so a
            # target C has not reached yet still has one.
            prep_row[1] if prep_row else species_of.of_reference(row[3]),
            row[1], row[2], row[3], row[4], row[5], row[8], row[13],
            prep_row[11] if prep_row else "(not prepped)",
            model_count.get(target, 0),
            row[14],
        ])
    _print(header, rows)


# ── models: F's measurements, one row per model ─────────────────────

def models(results):
    qc = os.path.join(results, *MODEL_QC)
    prep = os.path.join(results, *PREP)
    header = ["model", "species", "source", "chains", "residues", "metals",
              "ligands", "mean_conf", "conf_scale", "pct_low_conf", "verdict"]

    if not os.path.exists(qc):
        _print(header, [["(no models yet — F_Model_QC not run)"] + [NA] * 10])
        return

    # WORKING_FOLDER/RESULTS/protein_modeling/Models/ is one flat folder, so
    # the organism is not in the path of a model the way it is for the uploads.
    # It is put back here, from the target the file name carries.
    grouped = species_of.by_target(prep)

    rows = [[row[0], species_of.of_model(row[0], grouped),
             row[1], row[3], row[4], row[6], row[7], row[8], row[9],
             row[10], row[11]] for row in _table(qc, 12)]
    _print(header, rows)


# ── docking: I's scores, and the ranking within each receptor ───────

def docking(results):
    scores = os.path.join(results, *DOCKING)
    prep = os.path.join(results, *PREP)
    header = ["receptor", "species", "ligand", "best_kcal_mol",
              "rank_in_receptor", "basis", "box_size"]

    if not os.path.exists(scores):
        _print(header, [["(no docking runs yet — I_Docking not run)"] + [NA] * 6])
        return

    # A receptor is named after the model it was built from, so the organism
    # comes back the same way it does for a model.
    grouped = species_of.by_target(prep)

    # Rank within a receptor, because that is the only axis along which these
    # numbers can be compared at all.
    by_receptor = {}
    for row in _table(scores, 8):
        by_receptor.setdefault(row[0], []).append(row)

    rows = []
    for receptor in sorted(by_receptor):
        runs = by_receptor[receptor]
        try:
            runs.sort(key=lambda r: float(r[2]))
        except ValueError:
            pass                        # a failed run leaves a non-numeric cell
        for rank, row in enumerate(runs, start=1):
            rows.append([row[0], species_of.of_model(row[0], grouped),
                         row[1], row[2], rank, row[5], row[6]])
    _print(header, rows)


# ── selfcheck ───────────────────────────────────────────────────────

def selfcheck(_=None):
    failures = []

    padded = list(_table(os.devnull, 5))
    if padded:
        failures.append("_table should yield nothing for an empty file")

    if list(_rows("/nonexistent/path.tsv")):
        failures.append("_rows should yield nothing for a missing file")

    # A missing stage has to print a header plus one row, not raise: the
    # report is run while most of the pipeline is still unfinished. The table
    # goes to a buffer rather than the terminal so the check reports only
    # pass or fail.
    for table in (targets, models, docking):
        buffer = io.StringIO()
        try:
            with contextlib.redirect_stdout(buffer):
                table("/nonexistent/results")
        except Exception as error:
            failures.append(f"{table.__name__} raised on a missing stage: {error!r}")
            continue
        if len(buffer.getvalue().splitlines()) != 2:
            failures.append(f"{table.__name__} should print a header and one row")

    if failures:
        for failure in failures:
            print(f"FAIL {failure}")
        return 1
    print("OK all checks passed")
    return 0


COMMANDS = {
    "targets": targets,
    "models": models,
    "docking": docking,
    "selfcheck": selfcheck,
}


def main(argv):
    if not argv or argv[0] not in COMMANDS:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    command = argv[0]
    if command == "selfcheck":
        return selfcheck()

    if len(argv) < 2:
        print(f"usage: collect_protein_stats.py {command} <results_dir>", file=sys.stderr)
        return 2

    COMMANDS[command](argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
