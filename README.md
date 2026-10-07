# Ologist Bioinformatics Workshop

A hands-on, beginner-first bioinformatics workshop: you start with a laptop that has
nothing installed and finish able to work in a Linux terminal, manage software with
conda, inspect sequencing data, model and visualise proteins, and assemble and annotate
a metagenome.

Everything here assumes **no prior command-line experience**. Each module in the setup
guide ends with a checkpoint you can show an instructor before moving on.

---

## Start here

Two documents carry the whole workshop. Read them in this order.

| Document | What it is | When to use it |
| --- | --- | --- |
| [`Bioinformatics_Workshop_Setup_Instructions.pdf`](Bioinformatics_Workshop_Setup_Instructions.pdf) | 7 modules that install everything: Linux terminal, conda, the workshop environments, Git and VS Code | **Before** the workshop, at home, in order |
| [`Bioinformatics_Command_Cheatsheet.pdf`](Bioinformatics_Command_Cheatsheet.pdf) | Pocket reference: how to read a command, everyday shell commands, the workshop toolbox, PyMOL, the metagenomics pipeline, error messages | **During** the workshop, open beside the terminal |

The setup guide is organised as:

- **Part 1 — Get a Linux terminal.** Module 1 (Windows: Ubuntu via WSL) or Module 2 (Mac).
- **Part 2 — conda.** Module 3: Miniforge and the `protein_modeling` environment. Everyone.
- **Part 3 — Git, GitHub, VS Code.** Modules 4–6. Publishing to GitHub is optional.
- **Part 4 — Metagenomics.** Module 7: the `metagenomics_env` environment. Only if your session covers it.
- **Part 5 — Checkpoints and a troubleshooting table.**

Budget about **90 minutes** and ~10 GB of free disk, plus one Windows restart.
Module 7 adds roughly 40 minutes and ~20 GB more, mostly its annotation database.

> **Windows users:** the workshop runs inside **Ubuntu (WSL)**, not PowerShell and not the
> native Windows Miniforge Prompt. The tools come from bioconda, which only builds for
> Linux and macOS. Module 1 sets Ubuntu up for you; PowerShell is used for that one step
> and then left alone.

---

## Setup quickstart

For the full explanation of every line, follow the setup PDF. This is the short version
for someone who already has conda working in a Linux or macOS terminal.

```bash
git clone https://github.com/markcyril28/ctrl-alt-defeat-pathogens.git
cd ctrl-alt-defeat-pathogens

# The main environment — every session except metagenomics
bash setup_protein_modeling_conda_envs.sh

conda activate protein_modeling
python -c "import Bio, pandas, pymol; print('biopython', Bio.__version__, 'pandas', pandas.__version__)"
```

If that prints two version numbers and your prompt begins with `(protein_modeling)`,
you are ready.

Only if your session includes metagenomics:

```bash
bash setup_metagenomics_conda_envs.sh --with-db     # tools plus their databases
conda activate metagenomics_env
fastqc --version && flye --version && bakta --version
```

Both scripts are **idempotent** — they skip an environment that already exists, so
rerunning them is safe. Neither needs `sudo`, and `conda` should never be run with it.

---

## The two environments

One environment is active at a time. Switch with `conda deactivate`, then
`conda activate <name>`. A `command not found` is almost always the other environment
being active.

| | `protein_modeling` | `metagenomics_env` |
| --- | --- | --- |
| Built by | `setup_protein_modeling_conda_envs.sh` | `setup_metagenomics_conda_envs.sh` |
| Python | 3.11 | 3.10 |
| Tools | Biopython, pandas, PyMOL, FastQC, seqkit, bwa, samtools, bcftools, BLAST, AutoDock Vina, Meeko, RDKit | FastQC, NanoStat, fastp, fastplong, metaSPAdes, MetaFlye, metaQUAST, Bandage, Bakta, AMRFinderPlus, ABRicate |
| Used for | Sequence handling, alignment, variants, structure visualisation, docking | QC → cleaning → assembly → evaluation → graph → annotation → AMR screening |
| Needed by | Every session | The metagenomics session only |

They are deliberately **separate**: the assemblers pull in dependency sets large enough
to make conda downgrade packages the protein sessions need, and a toolbox that is not
also holding up everything else is far easier to delete and rebuild.

Both use the `conda-forge` and `bioconda` channels, in that priority order, and both
require Linux or macOS.

### Script options

```
setup_protein_modeling_conda_envs.sh
  --dry-run              print the conda command instead of running it
  -h, --help             usage

setup_metagenomics_conda_envs.sh
  --dry-run              print what would be created
  --with-db              also download Bakta (light), AMRFinderPlus and ABRicate/VFDB databases
  --db-dir=<path>        where to put them (default: ~/workshop/databases)
  -h, --help             usage
```

Each script verifies what it built by asking every tool for its version, and prints the
`conda activate` line to use. The metagenomics script falls back to installing the core
pipeline (QC through assembly graph) if the annotation and screening tools cannot be
solved, and exits non-zero if any check failed. OPERA-MS is **not** installed — it has no
conda package and is built from source only if a trainer asks.

A tool and its database are separate downloads, and a tool with no database fails only
when you run it. `--with-db` handles that up front; the summary prints the path to pass
to `bakta --db`.

---

## What is in this repository

```
Bioinformatics_Workshop_Setup_Instructions.pdf   install guide — read first
Bioinformatics_Command_Cheatsheet.pdf            command reference — keep open
setup_protein_modeling_conda_envs.sh             builds protein_modeling
setup_metagenomics_conda_envs.sh                 builds metagenomics_env
```

A full workshop checkout also carries the teaching material, which is distributed by the
instructors rather than tracked here:

```
Datasets/            sequences and structures for four bacterial species, plus a mock
                     metagenome community (short, long and hybrid reads) and reference
                     genomes; MANIFEST.tsv records provenance for every file
modules/             small Python helpers used in the protein-modelling and docking steps
RESULTS/             where session output lands, one folder per stage
docs/                manual sources and glossary
envs/                exported environment recipe
run_alphafold3_inputs.sh, run_alphafold3_*.toml
                     per-species AlphaFold 3 job definitions
```

Datasets carry one caveat worth repeating: `PDB_` files hold **experimental**
coordinates and `MODEL_` files hold **predicted** ones. In a predicted model the
B-factor column is per-residue pLDDT confidence, not a crystallographic B-factor. Never
mix the two in an analysis.

---

## If something goes wrong

The setup PDF (Part 5) and the cheatsheet (Part D) both end in a symptom → cause → fix
table. The four that come up most often:

| The screen says | Do this |
| --- | --- |
| `conda: command not found` | Close the terminal and open a new one — conda only exists in terminals started after it was installed |
| `command not found` for a workshop tool | Check which environment your prompt names, then `conda activate` the right one |
| `Killed` during assembly | The machine ran out of memory. Your command was not wrong — ask about a smaller input or the workshop server |
| A tool reports a missing database | The program installed, its database did not: `bakta_db list`, `amrfinder -u`, `abricate --list` |

When you ask for help, bring three things: **the exact command you typed, the full error
message, and the output of `conda info --envs`.** Copy the error before you start trying
fixes — the first error is the informative one and it is easy to lose.

---

## Scope and good habits

- These are **teaching commands**, not a validated production pipeline. They are written
  to make one tool at a time understandable.
- **Never write cleaned reads over raw ones.** Raw FASTQ is the one thing in a project you
  cannot regenerate. Keep inputs and outputs in separate folders.
- **A database hit is sequence similarity, nothing more.** It is not proof that a gene is
  expressed or that an organism is resistant, and a missing hit is not proof of absence.
  Report findings as "detected in this assembly, by this tool and database version".
- **Record versions as you go** — the tool version (`<tool> --version`) and the database
  version where one applies. The same command against a different database can give a
  different answer, so the version is part of your result.
- **Never publish** patient data, unpublished research data, credentials, `.env` files or
  large raw sequencing files to GitHub. Check your institution's data-sharing rules first.

---

## Credits

Workshop material prepared for the Ologist workshop series. The guides cite their
sources: Microsoft (WSL), conda-forge (Miniforge), Bioconda, Git, GitHub, VS Code, and
each bioinformatics tool's own documentation — see the Sources list at the end of the
setup PDF.
