# Ologist Bioinformatics Workshop

Bioinformatics from scratch — no command-line experience needed. By the end you can work
in a Linux terminal, install software with conda, assemble and annotate a metagenome, and
model and dock proteins.

## Two guides

| Guide | When |
| --- | --- |
| [Setup instructions](Bioinformatics_Workshop_Setup_Instructions.pdf) | **Before** the workshop — installs everything, step by step |
| [Command cheatsheet](Bioinformatics_Command_Cheatsheet.pdf) | **During** the workshop — keep it open beside the terminal |

Set aside about 2 hours and 30 GB of free disk space.

**Windows:** the workshop runs inside Ubuntu (WSL), not PowerShell. The setup guide
installs it in Module 1.

## Install — at home, not on the day

Follow the setup PDF; it explains every line. Short version, once conda works:

```bash
git clone https://github.com/markcyril28/ctrl-alt-defeat-pathogens.git
cd ctrl-alt-defeat-pathogens

# 1. Metagenomics tools — start here: biggest download, and the first session
bash setup_metagenomics_conda_envs.sh --with-db
conda activate meta_env
fastqc --version && flye --version && bakta --version

# 2. Protein tools
conda deactivate
bash setup_protein_modeling_conda_envs.sh
conda activate protein_modeling
python -c "import Bio, pandas, pymol; print('ok')"
```

If every command prints a version instead of an error, you are ready.

- Running a script twice is safe — it skips what already exists.
- Never use `sudo` with conda.
- Write down the Bakta database path the first script prints. Every `bakta` run needs it
  as `--db`.
- Useful flags: `--dry-run` (show, don't install), `--with-db` (download the databases),
  `--db-dir=<path>` (put them somewhere else).

## The two toolboxes

| | `meta_env` | `protein_modeling` |
| --- | --- | --- |
| Tools | FastQC, NanoStat, fastp, fastplong, metaSPAdes, MetaFlye, metaQUAST, Bandage, Bakta, AMRFinderPlus, ABRicate | Biopython, pandas, PyMOL, seqkit, bwa, samtools, bcftools, BLAST, AutoDock Vina |
| For | QC → cleaning → assembly → evaluation → graph → annotation → AMR screening | Sequences, alignment, variants, structures, docking |
| Cheatsheet | Section C6 | Sections C2–C5 |

They are kept separate on purpose: the assemblers would otherwise force conda to
downgrade the protein tools.

Only one is on at a time. Switch with `conda deactivate`, then `conda activate <name>`.
If a command is "not found", check which name your prompt shows first.

## If something breaks

| The screen says | Do this |
| --- | --- |
| `conda: command not found` | Close the terminal, open a new one |
| `command not found` for a tool | Activate the right environment |
| `Killed` during assembly | Out of memory — ask for a smaller input or the workshop server. Your command was fine |
| A missing database | `bakta_db list`, `amrfinder -u`, `abricate --list` |

More in Part 5 of the setup PDF and Part D of the cheatsheet.

Asking for help? Bring three things: the command you typed, the full error message, and
the output of `conda info --envs`.

## Good habits

- **Never overwrite raw reads.** They are the one thing you cannot regenerate.
- **A database hit is similarity, nothing more** — not proof of resistance, and a missing
  hit is not proof of absence.
- **Record versions** of both the tool and its database. They are part of your result.
- **Never put** patient data, passwords, or large raw sequencing files on GitHub.

## What's here

```
Bioinformatics_Workshop_Setup_Instructions.pdf   read first
Bioinformatics_Command_Cheatsheet.pdf            keep open
setup_metagenomics_conda_envs.sh                 builds meta_env
setup_protein_modeling_conda_envs.sh             builds protein_modeling
```

Datasets and finished outputs are not in this repo — they live in a shared folder [1].

[1]: https://drive.google.com/drive/folders/123zpKrLFz95MhzVqwrq0x5czZ7iv9Jd2?usp=drive_link
