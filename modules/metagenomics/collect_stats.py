#!/usr/bin/env python3
"""Turn the pipeline's scattered tool reports into the worksheet's tables.

Every tool in the metagenomics pipeline writes its own report in its own
format. The trainee worksheet asks for the same numbers arranged by stage.
This script does that rearranging and nothing else: it reads only files the
pipeline already wrote, and it computes no statistics of its own.

    collect_stats.py reads      <results_dir>   worksheet 2.4 and 5  — read stats
    collect_stats.py annotation <results_dir>   worksheet 13         — feature counts
    collect_stats.py screening  <results_dir>   worksheet 14         — AMR/virulence hits
    collect_stats.py selfcheck                  run the built-in checks

Each prints a TSV to stdout. A missing input is reported as a row, not as a
crash: a stage you have not run yet is a normal state, not an error.
"""
import glob
import json
import os
import sys

NA = "-"


def _rows(path, sep="\t"):
    """Yield the non-empty rows of a delimited file as lists of strings."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.strip():
                yield line.split(sep)


def _print(header, rows):
    print("\t".join(header))
    for row in rows:
        print("\t".join(str(c) for c in row))


# ── reads: fastp/fastplong JSON + NanoStat ──────────────────────────
# fastp reports before_filtering and after_filtering in one file, so a single
# JSON gives both the raw and the cleaned row of the worksheet's comparison.

def _fastp_rows(path, label):
    if not os.path.exists(path):
        return [[label, "(not run)"] + [NA] * 5]
    with open(path, encoding="utf-8") as fh:
        summary = json.load(fh)["summary"]
    out = []
    for when, key in (("raw", "before_filtering"), ("cleaned", "after_filtering")):
        s = summary[key]
        out.append([
            label, when,
            s.get("total_reads", NA),
            s.get("total_bases", NA),
            s.get("read1_mean_length", s.get("mean_length", NA)),
            s.get("q20_rate", NA),
            s.get("q30_rate", NA),
        ])
    return out


def _nanostat(path):
    """NanoStat prints 'Metric:   1,234.5' lines. Return {metric: value}."""
    stats = {}
    if not os.path.exists(path):
        return stats
    for line in open(path, encoding="utf-8", errors="replace"):
        if ":" in line:
            key, _, value = line.partition(":")
            value = value.strip().replace(",", "")
            if value:
                stats[key.strip().lower()] = value
    return stats


def cmd_reads(results):
    clean = os.path.join(results, "03_Read_Cleaning")
    rows = _fastp_rows(os.path.join(clean, "fastp.json"), "short_reads")
    rows += _fastp_rows(os.path.join(clean, "fastplong.json"), "long_reads")
    _print(["dataset", "state", "reads", "bases", "mean_length", "q20_rate", "q30_rate"], rows)

    # NanoStat adds what a JSON summary cannot: median length, N50, read quality.
    print()
    nano = [
        ("raw", os.path.join(results, "02_Read_QC", "raw", "long_reads_nanostat.txt")),
        ("cleaned", os.path.join(results, "02_Read_QC", "cleaned", "long_reads_nanostat.txt")),
    ]
    wanted = ["number of reads", "total bases", "mean read length",
              "median read length", "read length n50", "mean read quality",
              "median read quality"]
    rows = []
    for state, path in nano:
        stats = _nanostat(path)
        rows.append([state] + [stats.get(k, NA) for k in wanted])
    _print(["long_reads"] + [k.replace(" ", "_") for k in wanted], rows)


# ── annotation: Bakta's own summary file ────────────────────────────

def cmd_annotation(results):
    wanted = ["Length", "Count", "CDSs", "tRNAs", "tmRNAs", "rRNAs",
              "ncRNAs", "ncRNA regions", "CRISPR arrays", "pseudogenes",
              "hypotheticals", "signal peptides"]
    rows = []
    for summary in sorted(glob.glob(os.path.join(results, "09_Annotation", "*", "*.txt"))):
        assembly = os.path.basename(os.path.dirname(summary))
        stats = {}
        for line in open(summary, encoding="utf-8", errors="replace"):
            key, _, value = line.partition(":")
            if value.strip():
                stats[key.strip()] = value.strip()
        rows.append([assembly] + [stats.get(k, NA) for k in wanted])
    if not rows:
        rows = [["(program 09 has not run)"] + [NA] * len(wanted)]
    _print(["assembly"] + [k.replace(" ", "_").lower() for k in wanted], rows)


# ── screening: AMRFinderPlus and ABRicate, side by side ─────────────
# Both tools report a gene, an identity and a coverage; they disagree on
# column names, units (fraction vs percent) and order. Normalise, do not
# reinterpret: a hit here means sequence similarity, nothing more.

def _amrfinder_hits(path, assembly):
    rows = list(_rows(path))
    if len(rows) < 2:
        return []
    index = {name: i for i, name in enumerate(rows[0])}

    def col(row, *names):
        for name in names:
            i = index.get(name)
            if i is not None and i < len(row):
                return row[i]
        return NA

    out = []
    for row in rows[1:]:
        out.append([
            assembly, "amrfinderplus",
            col(row, "Scope", "scope"),
            col(row, "Gene symbol", "Element symbol"),
            col(row, "% Identity to reference sequence", "% Identity to reference"),
            col(row, "% Coverage of reference sequence", "% Coverage of reference"),
            col(row, "Element type"),
            col(row, "Sequence name", "Element name"),
        ])
    return out


def _abricate_hits(path, assembly):
    rows = list(_rows(path))
    if len(rows) < 2:
        return []
    index = {name.lstrip("#"): i for i, name in enumerate(rows[0])}

    def col(row, name):
        i = index.get(name)
        return row[i] if i is not None and i < len(row) else NA

    out = []
    for row in rows[1:]:
        out.append([
            assembly, "abricate",
            col(row, "DATABASE"),
            col(row, "GENE"),
            col(row, "%IDENTITY"),
            col(row, "%COVERAGE"),
            col(row, "RESISTANCE"),
            col(row, "PRODUCT"),
        ])
    return out


def cmd_screening(results):
    screen = os.path.join(results, "10_AMR_Virulence_Screening")
    rows = []
    for path in sorted(glob.glob(os.path.join(screen, "*.amrfinder.tsv"))):
        rows += _amrfinder_hits(path, os.path.basename(path).split(".")[0])
    for path in sorted(glob.glob(os.path.join(screen, "*.abricate.*.tsv"))):
        rows += _abricate_hits(path, os.path.basename(path).split(".")[0])
    if not rows:
        rows = [["(program 10 found no hits, or has not run)"] + [NA] * 7]
    _print(["assembly", "tool", "database", "gene", "identity", "coverage",
            "type", "product"], rows)


# ── self-check ──────────────────────────────────────────────────────

def cmd_selfcheck(_args=None):
    """Build one of every input file, then check the three readers."""
    import io
    import tempfile
    import contextlib

    with tempfile.TemporaryDirectory() as tmp:
        def write(rel, text):
            path = os.path.join(tmp, rel)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(text)

        write("03_Read_Cleaning/fastp.json", json.dumps({"summary": {
            "before_filtering": {"total_reads": 50000, "total_bases": 15050000,
                                 "read1_mean_length": 301, "q20_rate": 0.97, "q30_rate": 0.93},
            "after_filtering": {"total_reads": 48000, "total_bases": 14000000,
                                "read1_mean_length": 292, "q20_rate": 0.99, "q30_rate": 0.96}}}))
        write("02_Read_QC/raw/long_reads_nanostat.txt",
              "General summary:\nNumber of reads:\t99,597.0\nTotal bases:\t758,531,352.0\n"
              "Mean read length:\t7,616.0\nMedian read length:\t5,120.0\n"
              "Read length N50:\t11,284.0\nMean read quality:\t12.3\nMedian read quality:\t13.1\n")
        write("09_Annotation/metaspades/metaspades.txt",
              "Length: 20123456\nCount: 812\nCDSs: 19877\ntRNAs: 142\nrRNAs: 9\n"
              "pseudogenes: 31\nhypotheticals: 4201\n")
        write("10_AMR_Virulence_Screening/metaspades.amrfinder.tsv",
              "Name\tGene symbol\t% Identity to reference sequence\t"
              "% Coverage of reference sequence\tElement type\tSequence name\tScope\n"
              "ctg1\tblaSHV-11\t99.88\t100.00\tAMR\tclass A beta-lactamase\tcore\n")
        write("10_AMR_Virulence_Screening/metaspades.abricate.vfdb.tsv",
              "#FILE\tSEQUENCE\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tPRODUCT\tRESISTANCE\n"
              "ctg.fa\tctg1\tfimA\t98.5\t99.1\tvfdb\ttype 1 fimbriae\t.\n")

        def capture(fn):
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                fn(tmp)
            return buf.getvalue()

        out = capture(cmd_reads)
        assert "short_reads\traw\t50000" in out, out
        assert "short_reads\tcleaned\t48000" in out, out
        assert "long_reads\t(not run)" in out, out           # fastplong.json absent
        assert "raw\t99597.0\t758531352.0" in out, out        # commas stripped
        assert "cleaned\t-\t-" in out, out                    # stage 03 absent

        out = capture(cmd_annotation)
        assert "metaspades\t20123456\t812\t19877\t142" in out, out
        assert "\t-\t" in out, out                            # tmRNAs missing -> "-"

        out = capture(cmd_screening)
        assert "metaspades\tamrfinderplus\tcore\tblaSHV-11\t99.88\t100.00\tAMR" in out, out
        assert "metaspades\tabricate\tvfdb\tfimA\t99.1\t98.5" in out, out  # columns reordered
        assert out.count("\n") == 3, out                      # header + 2 hits

    # An empty results directory must produce tables, not tracebacks.
    with tempfile.TemporaryDirectory() as empty:
        for fn in (cmd_reads, cmd_annotation, cmd_screening):
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                fn(empty)
            assert buf.getvalue().strip(), fn.__name__

    print("collect_stats.py: all checks passed")


if __name__ == "__main__":
    COMMANDS = {
        "reads": cmd_reads,
        "annotation": cmd_annotation,
        "screening": cmd_screening,
        "selfcheck": cmd_selfcheck,
    }
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit("usage: collect_stats.py {reads|annotation|screening} <results_dir>\n"
                 "       collect_stats.py selfcheck")
    if sys.argv[1] == "selfcheck":
        cmd_selfcheck()
    elif len(sys.argv) < 3:
        sys.exit("collect_stats.py %s needs the results directory" % sys.argv[1])
    else:
        COMMANDS[sys.argv[1]](sys.argv[2])
