#!/usr/bin/env python3
"""Write AlphaFold3 job files from a protein FASTA.

AlphaFold3 is driven by JSON, and there are two dialects that are not
interchangeable:

    alphafoldserver   what the AlphaFold Server web form accepts. A file is a
                      JSON *array* of jobs, chains have no IDs, and an empty
                      modelSeeds list lets the server pick the seeds.
    alphafold3        what a local AlphaFold3 install accepts. A file is a
                      single job object, every chain carries an explicit id,
                      and the seeds are yours to state.

    make_af3_json.py <fasta-or-folder> <outdir> <dialect> [options]

      --seeds N         model seeds to request (default 1). 0 writes an empty
                        list, which only the server dialect accepts.
      --ions "LIST"     add an ion chain for each code in a space-separated
                        list, e.g. --ions "ZN MG". "" adds none.
      --ligands "LIST"  add a ligand chain for each CCD code in a
                        space-separated list, e.g. --ligands "NAD". "" adds none.
      --batch-size N    also write batch files of N jobs each (server dialect
                        only, where one upload can carry several jobs).
      --batch-name NAME stem for those batch files (default AF3_server_batch).

One FASTA, or a folder of them. Given a folder, every .faa in it is taken as
one organism — the file's name is the organism's — and its jobs are written to
<outdir>/<organism>/, with any batch file named after the organism too. That
is how program E calls it, from program C's by_species/ folder: a server batch
then carries one organism per upload, so what comes back is already grouped,
and a job that fails belongs to an organism rather than to a pile of fifteen.

Given a single FASTA the jobs go straight into <outdir> and the species column
reads "-", which is what you want when you are asking about one sequence.

A TSV manifest goes to stdout: nothing is printed that the caller has to parse
out of a log.

Note on what this does and does not do: it writes inputs. It does not predict
anything. The ions and ligands named here are a hypothesis you are asking
AlphaFold3 to place — asking for a zinc does not mean the protein binds one,
and leaving it out does not mean it does not.
"""
import glob
import json
import os
import sys

# Chain IDs for the local dialect: the protein is A, extra chains follow.
CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

HEADER = ["name", "species", "dialect", "length", "ions", "ligands", "file"]

NO_SPECIES = "-"


def _read_fasta(path):
    """Yield (id, sequence) pairs. The ID is the header up to the first space."""
    name, chunks = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def _groups(source):
    """(species, fasta path, outdir suffix) for each input to write jobs from.

    A folder is one group per .faa inside it; a file is one group with no
    species and nothing added to the output path.
    """
    if os.path.isdir(source):
        return [(os.path.basename(path)[:-4], path, os.path.basename(path)[:-4])
                for path in sorted(glob.glob(os.path.join(source, "*.faa")))]
    return [(NO_SPECIES, source, "")]


def _server_job(name, sequence, seeds, ions, ligands):
    sequences = [{"proteinChain": {"sequence": sequence, "count": 1}}]
    for code in ions:
        # An ion is named by its bare code in this dialect: {"ion": "ZN"}.
        sequences.append({"ion": {"ion": code, "count": 1}})
    for code in ligands:
        # A ligand is not. The server wants the CCD code prefixed — NAD is
        # written "CCD_NAD" — and rejects the upload without it. The local
        # dialect's ccdCodes list takes the bare code instead, which is why
        # the same ligand is spelled two ways in this file.
        prefixed = code if code.startswith("CCD_") else f"CCD_{code}"
        sequences.append({"ligand": {"ligand": prefixed, "count": 1}})
    return {
        "name": name,
        "modelSeeds": seeds,
        "sequences": sequences,
        "dialect": "alphafoldserver",
        "version": 1,
    }


def _local_job(name, sequence, seeds, ions, ligands):
    sequences = [{"protein": {"id": "A", "sequence": sequence}}]
    # Each ion and ligand is its own chain, so each needs its own ID. Both are
    # "ligand" entries here, and both take the bare CCD code.
    for offset, code in enumerate(list(ions) + list(ligands)):
        chain = CHAIN_IDS[(offset + 1) % len(CHAIN_IDS)]
        bare = code[4:] if code.startswith("CCD_") else code
        sequences.append({"ligand": {"id": chain, "ccdCodes": [bare]}})
    return {
        "name": name,
        "modelSeeds": seeds or [1],
        "sequences": sequences,
        "dialect": "alphafold3",
        "version": 1,
    }


def _write(path, payload):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
        fh.write("\n")


def main(argv):
    if len(argv) < 3:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    source, outdir, dialect = argv[0], argv[1], argv[2]
    if dialect not in ("alphafoldserver", "alphafold3"):
        print(f"unknown dialect: {dialect}", file=sys.stderr)
        return 2

    seed_count, ions, ligands = 1, [], []
    batch_size, batch_name = 0, "AF3_server_batch"

    rest = argv[3:]
    while rest:
        flag = rest.pop(0)
        if flag == "--seeds":
            seed_count = int(rest.pop(0))
        elif flag == "--ions":
            ions.extend(rest.pop(0).split())
        elif flag == "--ligands":
            ligands.extend(rest.pop(0).split())
        elif flag == "--batch-size":
            batch_size = int(rest.pop(0))
        elif flag == "--batch-name":
            batch_name = rest.pop(0)
        else:
            print(f"unknown option: {flag}", file=sys.stderr)
            return 2

    if seed_count == 0 and dialect == "alphafold3":
        print("a local AlphaFold3 job needs at least one seed", file=sys.stderr)
        return 2

    groups = _groups(source)
    if not groups:
        print(f"no .faa files in {source}", file=sys.stderr)
        return 1

    seeds = list(range(1, seed_count + 1))
    print("\t".join(HEADER))

    for species, fasta_path, suffix in groups:
        folder = os.path.join(outdir, suffix) if suffix else outdir
        os.makedirs(folder, exist_ok=True)

        jobs = []
        for name, sequence in _read_fasta(fasta_path):
            if not sequence:
                print(f"skipping {name}: no sequence", file=sys.stderr)
                continue

            if dialect == "alphafoldserver":
                job = _server_job(name, sequence, seeds, ions, ligands)
                payload = [job]      # even one job is an array in this dialect
            else:
                job = _local_job(name, sequence, seeds, ions, ligands)
                payload = job

            path = os.path.join(folder, f"{name}.json")
            _write(path, payload)
            jobs.append(job)
            print("\t".join([
                name, species, dialect, str(len(sequence)),
                ",".join(ions) or "-", ",".join(ligands) or "-", path,
            ]))

        # One upload can carry several jobs, which is the difference between
        # pasting twelve forms and pasting one. The organism is in the file
        # name as well as the folder, because this is a file you pick out of a
        # file dialog and then see again in the server's job list.
        if batch_size > 0 and dialect == "alphafoldserver" and jobs:
            stem = f"{batch_name}_{species}" if suffix else batch_name
            for start in range(0, len(jobs), batch_size):
                chunk = jobs[start:start + batch_size]
                number = start // batch_size + 1
                path = os.path.join(folder, f"{stem}_{number:02d}.json")
                _write(path, chunk)
                print("\t".join([
                    f"{stem}_{number:02d}", species, dialect, "-",
                    str(len(chunk)) + " jobs", "-", path,
                ]))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
