#!/usr/bin/env bash
# run_protein_modeling.sh — scripted Part I/II of the Protein Modeling Workshop
#
# Runs the reproducible (non-GUI) half of docs/Drafts/Protein_Modeling_Workshop_Manual.md:
# load the structure dataset, export clean receptors, obtain models for the targets
# that have only a sequence, and judge how much each model can be trusted.
#
# Driven entirely by the CONTROL PANEL variables below — this script takes no
# command-line arguments. Set STEP, then run:  bash run_protein_modeling.sh
#
# Steps (set STEP to one of these):
#   session      Load every PDB, colour by chain, save a combined PyMOL session   (Step 4)
#   clean        Export docking-ready receptors, one per target recipe            (Step 6)
#   check        Report which targets have an AlphaFold DB model                  (Step 8a)
#   fetch        Download the available AlphaFold models into Datasets/MODELS     (Step 8b)
#   plddt        Mean pLDDT per model and per region of interest                  (Step 8c)
#   validate     Superpose each model on its crystal structure (super + align)    (Step 8d)
#   identity     Sequence identity of BoNT/C1 to its BoNT/A and /B templates      (Step 8e)
#   trim         Cut low-confidence tails and linkers off the models              (Step 8f)
#   zn-transfer  Give the BoNT/E and /F models the catalytic Zn they lack         (Step 10)
#   all          session, clean, check, fetch, plddt, validate, identity, trim, zn-transfer
#
# Models — fetched, trimmed and Zn-loaded — always use the MODEL_ prefix, never PDB_:
# predicted and experimental coordinates carry different kinds of error.
#
# Needs: pymol (open-source), python 3 with biopython. Docking tools are not used here —
# see run_docking.sh for those.

set -euo pipefail

# ────────────────────────────────────────────────────────────────────
# CONTROL PANEL — set the run here; the script takes no arguments.
# ────────────────────────────────────────────────────────────────────

# ── What to run ─────────────────────────────────────────────────────
# session | clean | check | fetch | plddt | validate | identity | trim | zn-transfer | all
STEP=""

# ── Environment and paths ───────────────────────────────────────────
ENV_WANTED="${OLOGIST_ENV:-}"   # conda env name; OLOGIST_ENV wins, empty = autodetect
DATASETS=""             # dataset folder; empty = <script dir>/Datasets
OUT_DIR=""              # receptors/sessions output folder; empty = <script dir>/work

# ── Behavior toggles ────────────────────────────────────────────────
DRY_RUN=false           # true = print what would run, touch nothing
# ────────────────────────────────────────────────────────────────────

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -z "$DATASETS" ]] && DATASETS="$ROOT/Datasets"
[[ -d "$DATASETS" ]] && DATASETS="$(cd "$DATASETS" && pwd)"
[[ -z "$OUT_DIR"  ]] && OUT_DIR="$ROOT/work"

# ── Output helpers (same vocabulary as setup_conda_envs.sh) ─────────
GREEN='\033[0;32m'  YELLOW='\033[0;33m'  RED='\033[0;31m'  BLUE='\033[0;34m'  NC='\033[0m'
info()  { printf "${GREEN}[OK]${NC}   %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC} %s\n" "$*"; }
fail()  { printf "${RED}[FAIL]${NC} %s\n" "$*" >&2; }
step()  { printf "\n${BLUE}━━━ %s ━━━${NC}\n" "$*"; }

# ── Validate the control-panel settings ─────────────────────────────
case "$STEP" in
    session|clean|check|fetch|plddt|validate|identity|trim|zn-transfer|all) ;;
    "") fail "STEP is not set. Edit the CONTROL PANEL at the top of this script."
        fail "Choose one of: session clean check fetch plddt validate identity trim zn-transfer all"
        exit 1 ;;
    *)  fail "unknown STEP: '$STEP'"
        fail "Choose one of: session clean check fetch plddt validate identity trim zn-transfer all"
        exit 1 ;;
esac

RECEPTORS="$OUT_DIR/receptors"
TMP=""
cleanup() { [[ -n "$TMP" && -d "$TMP" ]] && rm -rf "$TMP"; }
trap cleanup EXIT

# ── Environment resolution ──────────────────────────────────────────
# Tools are used straight from the environment's bin directory, and that directory is
# prepended to PATH so helper scripts see it too. This avoids the `conda run`
# activation bugs noted in setup_conda_envs.sh.
tools_present() {
    local dir="$1" tool
    for tool in pymol python; do
        if [[ -n "$dir" ]]; then [[ -x "$dir/$tool" ]] || return 1
        else command -v "$tool" &>/dev/null || return 1; fi
    done
}

resolve_env() {
    if [[ -z "$ENV_WANTED" ]] && tools_present ""; then
        info "using pymol/python already on PATH ($(command -v pymol))"
        return
    fi
    local candidates=() name prefix
    [[ -n "$ENV_WANTED" ]] && candidates+=("$ENV_WANTED")
    candidates+=(protein_model protein_modeling dock-workshop)
    for name in "${candidates[@]}"; do
        prefix="$(conda info --envs 2>/dev/null | awk -v n="$name" '$1==n {print $NF}')" || true
        [[ -z "$prefix" ]] && continue
        if tools_present "$prefix/bin"; then
            export PATH="$prefix/bin:$PATH"
            info "using conda environment '$name' ($prefix)"
            return
        fi
        warn "environment '$name' exists but lacks pymol or python"
    done
    fail "no environment with pymol + python found. Run: bash setup_conda_envs.sh"
    fail "then set ENV_WANTED in the CONTROL PANEL if your environment is named differently."
    exit 1
}

# Run a PyMOL or plain-Python script supplied as a quoted heredoc on stdin. Every such
# script runs with the Datasets folder as its working directory, because the dataset
# paths and the helper scripts are all relative to it.
run_pymol() {
    local script="$TMP/${1}.py"
    cat > "$script"
    if $DRY_RUN; then printf '  [dry-run] pymol -cq %s\n' "$script"; return 0; fi
    ( cd "$DATASETS" && pymol -cq "$script" )
}

run_python() {
    local script="$TMP/${1}.py"
    cat > "$script"
    if $DRY_RUN; then printf '  [dry-run] python %s\n' "$script"; return 0; fi
    ( cd "$DATASETS" && python "$script" )
}

preflight() {
    local species_dirs=("$DATASETS"/[1-9]_*/)
    [[ ${#species_dirs[@]} -gt 0 ]] || { fail "no species directories in $DATASETS"; exit 1; }
    resolve_env
    TMP="$(mktemp -d)"
    $DRY_RUN || mkdir -p "$RECEPTORS"
    export OW_DATASETS="$DATASETS" OW_RECEPTORS="$RECEPTORS"
}

# ── Step 4 — combined session ───────────────────────────────────────
do_session() {
    step "Step 4 — load the whole dataset into one session"
    run_pymol session <<'PY'
import os
from pymol import cmd
out = os.path.join(os.environ["OW_RECEPTORS"], "toxin_dataset.pse")
cmd.do("run toxin_load.py")
cmd.do("load_all_pdbs()")
cmd.do("color_by_chain_all()")
loaded = cmd.get_object_list()
cmd.save(out)
print("loaded {} structures: {}".format(len(loaded), ", ".join(sorted(loaded))))
print("session saved to", out)
PY
    info "session written to $RECEPTORS/toxin_dataset.pse"
}

# ── Step 6 — clean receptors ────────────────────────────────────────
# Each recipe is (name, pdb, keep selection, HET residues to drop, expected Zn, outfile).
# The keep selection is what survives; everything else, plus solvent and hydrogens, goes.
do_clean() {
    step "Step 6 — export docking-ready receptors"
    run_pymol clean <<'PY'
import os
from pymol import cmd

OUT = os.environ["OW_RECEPTORS"]

RECIPES = [
    # name,                pdb path (relative to Datasets/),                                    keep selection,                              drop resn,        Zn, outfile
    ("TNT",                "1_Mycobacterium_tuberculosis/PDB_4QLP_TNT_immunity.pdb",             "chain B and polymer",                       "",                0, "TNT_receptor_clean.pdb"),
    ("BoNT/A holotoxin",   "2_Clostridium_botulinum/PDB_3BTA_BoNT_A_holotoxin.pdb",              "(chain A and polymer) or resn ZN",          "",                1, "BONT_A_holo_clean.pdb"),
    ("BoNT/A LC",          "2_Clostridium_botulinum/PDB_1XTG_BoNT_A_LC_SNAP25.pdb",              "(chain A and polymer) or resn ZN",          "CL",              1, "BONT_A_LC_clean.pdb"),
    ("BoNT/B catalytic",   "2_Clostridium_botulinum/PDB_1EPW_BoNT_B.pdb",                        "polymer or resn ZN",                        "SO4",             1, "BONT_B_catalytic_clean.pdb"),
    ("BoNT/B receptor-bd", "2_Clostridium_botulinum/PDB_1I1E_BoNT_B_doxorubicin.pdb",            "polymer",                                   "DM2+SO4",         0, "BONT_B_HC_clean.pdb"),
    ("FimH lectin",        "3_Klebsiella_pneumoniae/PDB_9AT9_FimH_lectin_mannose.pdb",           "polymer",                                   "MAN",             0, "FimH_receptor_clean.pdb"),
    ("ExoA",               "4_Pseudomonas_aeruginosa/PDB_1AER_ExoA.pdb",                        "chain A and polymer",                       "TAD+TIA+AMP",     0, "ExoA_receptor_clean.pdb"),
    ("ExoT–SpcS interface","4_Pseudomonas_aeruginosa/PDB_6JNP_ExoT_SpcS_complex.pdb",           "(chain A or chain B) and polymer",          "GOL",             0, "ExoT_SpcS_clean.pdb"),
]

for name, pdb, keep, drop, want_zn, outfile in RECIPES:
    path = pdb
    if not os.path.exists(path):
        print("SKIP {}: {} not found".format(name, path))
        continue
    cmd.delete("all")
    cmd.load(path, "rec")
    before = cmd.count_atoms("rec")
    cmd.remove("rec and not ({})".format(keep))
    if drop:
        cmd.remove("rec and resn {}".format(drop))
    cmd.remove("rec and solvent")
    cmd.remove("rec and hydro")
    zinc = cmd.count_atoms("rec and resn ZN")
    target = os.path.join(OUT, outfile)
    cmd.save(target, "rec")
    print("{:22s} {:5d} -> {:5d} heavy atoms, Zn {}  ->  {}".format(
        name, before, cmd.count_atoms("rec"), zinc, outfile))
    # The manual's warning, enforced: a metalloprotease without its metal is not the
    # enzyme you meant to study.
    if zinc != want_zn:
        print("   WARNING: expected {} Zn atom(s), kept {}".format(want_zn, zinc))

print("\nNote: BoNT/B receptor-binding keeps protein only — its target site is the")
print("doxorubicin pocket, 81 A from the catalytic Zn, so the metal is not needed there.")
print("For any catalytic-site receptor the Zn must survive; the counts above are the check.")
PY
    info "receptors written to $RECEPTORS/"
}

# ── Step 8a/8b — AlphaFold models ───────────────────────────────────
do_check() {
    step "Step 8a — which targets need a model"
    if $DRY_RUN; then printf '  [dry-run] python fetch_models.py --check\n'; return; fi
    ( cd "$DATASETS" && python fetch_models.py --check )
}

do_fetch() {
    step "Step 8b — download the available AlphaFold models"
    if $DRY_RUN; then printf '  [dry-run] python fetch_models.py\n'; return; fi
    ( cd "$DATASETS" && python fetch_models.py )
    local model_count; model_count="$(find "$DATASETS" -name 'MODEL_*.pdb' | wc -l)"
    info "models across species dirs ($model_count files)"
    warn "BoNT/C1, /D and /G have no AlphaFold entry — build them by homology (Step 8e)."
}

# ── Step 8c — confidence ────────────────────────────────────────────
# pLDDT lives in the B-factor column of an AlphaFold file. Judge it over the site you
# intend to dock into, not over the whole chain.
do_plddt() {
    step "Step 8c — mean pLDDT, whole chain and per region"
    run_python plddt <<'PY'
import os, statistics

DATASETS = os.environ["OW_DATASETS"]

REGIONS = [
    # path (relative to Datasets/),                                          label,                     first, last
    ("1_Mycobacterium_tuberculosis/MODEL_O05442_CpnT_full_AF.pdb",         "whole chain",              None, None),
    ("1_Mycobacterium_tuberculosis/MODEL_O05442_CpnT_full_AF.pdb",         "TNT domain 651-846",        651,  846),
    ("1_Mycobacterium_tuberculosis/MODEL_O05442_CpnT_full_AF.pdb",         "linker 397-634",            397,  634),
    ("4_Pseudomonas_aeruginosa/MODEL_Q9I788_ExoT_AF.pdb",                  "whole chain",              None, None),
    ("4_Pseudomonas_aeruginosa/MODEL_Q9I788_ExoT_AF.pdb",                  "ADPRT domain 236-457",      236,  457),
    ("4_Pseudomonas_aeruginosa/MODEL_Q9I788_ExoT_AF.pdb",                  "GAP domain 78-235",          78,  235),
    ("4_Pseudomonas_aeruginosa/MODEL_Q9I788_ExoT_AF.pdb",                  "chaperone-binding 23-79",    23,   79),
    ("4_Pseudomonas_aeruginosa/MODEL_G3XDA1_ExoS_AF.pdb",                  "whole chain",              None, None),
    ("4_Pseudomonas_aeruginosa/MODEL_G3XDA1_ExoS_AF.pdb",                  "ADPRT domain 232-453",      232,  453),
    ("3_Klebsiella_pneumoniae/MODEL_P21648_MrkD_AF.pdb",                    "whole chain",              None, None),
    ("3_Klebsiella_pneumoniae/MODEL_P21648_MrkD_AF.pdb",                    "minus signal region 8-321",   8,  321),
    ("3_Klebsiella_pneumoniae/MODEL_A0A0H3H2I8_FimH_full_AF.pdb",          "whole chain",              None, None),
    ("3_Klebsiella_pneumoniae/MODEL_A0A0H3H2I8_FimH_full_AF.pdb",          "lectin domain 1-160",         1,  160),
    ("2_Clostridium_botulinum/MODEL_Q00496_BoNT_E_AF.pdb",                 "whole chain",              None, None),
    ("2_Clostridium_botulinum/MODEL_P30996_BoNT_F_AF.pdb",                 "whole chain",              None, None),
    ("4_Pseudomonas_aeruginosa/MODEL_O34208_ExoU_AF.pdb",                  "whole chain",              None, None),
    ("4_Pseudomonas_aeruginosa/MODEL_O34208_ExoU_AF.pdb",                  "PLA2 region 107-357",       107,  357),
]

def per_residue_plddt(path):
    scores = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("ATOM"):
                scores[int(line[22:26])] = float(line[60:66])
    return scores


def verdict(mean):
    if mean >= 90: return "very high — backbone and side chains reliable"
    if mean >= 70: return "confident backbone"
    if mean >= 50: return "low — treat with care"
    return "no structure predicted here"

print("{:55s} {:26s} {:>6s} {:>6s}  reading".format("model", "region", "pLDDT", "res"))
print("-" * 116)
missing = []
for relpath, label, first, last in REGIONS:
    path = relpath
    if not os.path.exists(path):
        if relpath not in missing:
            missing.append(relpath)
        continue
    scores = per_residue_plddt(path)
    if first is not None:
        scores = {r: b for r, b in scores.items() if first <= r <= last}
    if not scores:
        print("{:30s} {:26s} {:>6s}".format(filename, label, "n/a"))
        continue
    mean = statistics.mean(scores.values())
    print("{:55s} {:26s} {:6.1f} {:6d}  {}".format(
        os.path.basename(relpath), label, mean, len(scores), verdict(mean)))

for relpath in missing:
    print("SKIP (not downloaded yet):", relpath)
print("\nA low whole-chain mean does not condemn a model: CpnT averages ~65 because of a")
print("disordered linker you will delete anyway, while its TNT domain scores ~90.")
PY
}

# ── Step 8d — model vs crystal ──────────────────────────────────────
do_validate() {
    step "Step 8d — superpose each model on its crystal structure"
    run_pymol validate <<'PY'
import os
from pymol import cmd

PAIRS = [
    # label,        model path,                                                      model selection,  crystal path,                                                crystal selection
    ("TNT domain",  "1_Mycobacterium_tuberculosis/MODEL_O05442_CpnT_full_AF.pdb",   "resi 651-846",   "1_Mycobacterium_tuberculosis/PDB_4QLP_TNT_immunity.pdb",     "chain B and polymer"),
    ("FimH lectin", "3_Klebsiella_pneumoniae/MODEL_A0A0H3H2I8_FimH_full_AF.pdb",    "polymer",        "3_Klebsiella_pneumoniae/PDB_9AT9_FimH_lectin_mannose.pdb",  "polymer"),
]

print("{:14s} {:>12s} {:>7s} {:>12s} {:>7s}".format(
    "comparison", "super RMSD", "atoms", "align RMSD", "atoms"))
print("-" * 58)
for label, model_path, model_sel, xtal_file, xtal_sel in PAIRS:
    if not (os.path.exists(model_path) and os.path.exists(xtal_file)):
        print("SKIP {}: run the fetch step first".format(label))
        continue
    cmd.delete("all")
    cmd.load(model_path, "pred_src")
    cmd.load(xtal_file, "xtal_src")
    cmd.create("pred1", "pred_src and ({})".format(model_sel))
    cmd.create("xtal1", "xtal_src and ({})".format(xtal_sel))
    sup_rms, sup_atoms = cmd.super("pred1", "xtal1")[:2]
    # `super` is sequence-independent and rejects outliers; `align` rejects nothing.
    # Running both is the habit that catches a good number bought by discarding
    # the parts that disagreed.
    cmd.create("pred2", "pred_src and ({})".format(model_sel))
    cmd.create("xtal2", "xtal_src and ({})".format(xtal_sel))
    ali_rms, ali_atoms = cmd.align("pred2", "xtal2")[:2]
    print("{:14s} {:12.2f} {:7d} {:12.2f} {:7d}".format(
        label, sup_rms, sup_atoms, ali_rms, ali_atoms))

print("\nWhere a crystal structure of your site exists, dock into the crystal.")
print("Use models only for the sites where nothing else exists.")
PY
}

# ── Step 8e — template identity ─────────────────────────────────────
do_identity() {
    step "Step 8e — identity of BoNT/C1 to its candidate templates"
    run_python identity <<'PY'
from pathlib import Path

try:
    from Bio import Align
except ImportError:
    raise SystemExit("biopython is missing. conda install -c conda-forge biopython")


def sequence(path):
    return "".join(line.strip() for line in Path(path).read_text().splitlines()
                   if not line.startswith(">"))

aligner = Align.PairwiseAligner(scoring="blastp", mode="global")
queries = [
    ("BoNT/C1", "2_Clostridium_botulinum/FASTA_P18640_BoNT_C1.fasta"),
    ("BoNT/D",  "2_Clostridium_botulinum/FASTA_P19321_BoNT_D.fasta"),
    ("BoNT/G",  "2_Clostridium_botulinum/FASTA_Q60393_BoNT_G.fasta"),
]
templates = [
    ("BoNT/B (1EPW)", "2_Clostridium_botulinum/FASTA_P10844_BoNT_B.fasta"),
    ("BoNT/A1 (3BTA)", "2_Clostridium_botulinum/FASTA_P0DPI1_BoNT_A1.fasta"),
]

for qname, qpath in queries:
    if not Path(qpath).exists():
        print("SKIP {}: {} not found".format(qname, qpath))
        continue
    target = sequence(qpath)
    for tname, tpath in templates:
        alignment = aligner.align(target, sequence(tpath))[0]
        top, bottom = alignment[0], alignment[1]
        pairs = [(x, y) for x, y in zip(top, bottom) if x != "-" and y != "-"]
        identity = 100 * sum(x == y for x, y in pairs) / len(pairs)
        band = ("good enough to locate a binding site" if identity >= 50 else
                "locates domains; check every conclusion" if identity >= 30 else
                "fold hypothesis only")
        print("{:8s} vs {:16s} {:3.0f}% over {} aligned positions — {}".format(
            qname, tname, identity, len(pairs), band))

print("\nThese three serotypes have no AlphaFold DB entry. Build them at")
print("https://swissmodel.expasy.org/interactive using the templates above, and save")
print("the result as <species>/<ACCESSION>_<name>_SWISSMODEL.pdb with its QMEANDisCo score.")
print("Always align before estimating identity: comparing these sequences position by")
print("position without alignment reports ~9%, which is chance level and simply wrong.")
PY
}

# ── Step 8f — trim models ───────────────────────────────────────────
do_trim() {
    step "Step 8f — trim low-confidence regions off the models"
    run_pymol trim <<'PY'
import os
from pymol import cmd

# (source path, keep selection, label, output dir, outfile)
# Residues below pLDDT 50 are dropped everywhere: they are not structure,
# and they add false surface that can swallow a ligand.
TRIMS = [
    ("1_Mycobacterium_tuberculosis/MODEL_O05442_CpnT_full_AF.pdb",   "resi 651-846", "CpnT TNT domain",  "1_Mycobacterium_tuberculosis", "MODEL_CpnT_TNT_domain_trimmed.pdb"),
    ("3_Klebsiella_pneumoniae/MODEL_P21648_MrkD_AF.pdb",             "not resi 1-7", "MrkD",             "3_Klebsiella_pneumoniae",      "MODEL_MrkD_trimmed.pdb"),
    ("4_Pseudomonas_aeruginosa/MODEL_G3XDA1_ExoS_AF.pdb",           "not resi 1-96","ExoS GAP+ADPRT",   "4_Pseudomonas_aeruginosa",     "MODEL_ExoS_trimmed.pdb"),
    ("4_Pseudomonas_aeruginosa/MODEL_Q9I788_ExoT_AF.pdb",           "resi 236-457", "ExoT ADPRT",       "4_Pseudomonas_aeruginosa",     "MODEL_ExoT_ADPRT_trimmed.pdb"),
    ("4_Pseudomonas_aeruginosa/MODEL_O34208_ExoU_AF.pdb",           "resi 107-357", "ExoU PLA2 region", "4_Pseudomonas_aeruginosa",     "MODEL_ExoU_PLA2_trimmed.pdb"),
]

for srcpath, keep, label, outdir, outfile in TRIMS:
    if not os.path.exists(srcpath):
        print("SKIP {}: {} not downloaded".format(label, srcpath))
        continue
    cmd.delete("all")
    cmd.load(srcpath, "pred")
    before = cmd.count_atoms("pred and name CA")
    cmd.create("trimmed", "pred and ({})".format(keep))
    cmd.remove("trimmed and b < 50")
    kept = cmd.count_atoms("trimmed and name CA")
    first = cmd.get_model("trimmed and name CA").atom[0].resi if kept else "-"
    last = cmd.get_model("trimmed and name CA").atom[-1].resi if kept else "-"
    cmd.save(os.path.join(outdir, outfile), "trimmed")
    print("{:18s} {:4d} -> {:4d} residues (now {}-{})  ->  {}/{}".format(
        label, before, kept, first, last, outdir, outfile))

print("\nWrite down what you removed. A pose reported against residue numbering you")
print("silently altered is not reproducible.")
print("ExoU is the weakest model in the set (PLA2 region ~72) — interpret it cautiously.")
PY
}

# ── Step 10 note — Zn for the BoNT/E and /F models ──────────────────
do_zn_transfer() {
    step "Step 10 — transfer the catalytic Zn into the BoNT/E and /F models"
    run_pymol zn_transfer <<'PY'
import os
from pymol import cmd

TEMPLATE = "2_Clostridium_botulinum/PDB_1XTG_BoNT_A_LC_SNAP25.pdb"
BONT_MODELS = "2_Clostridium_botulinum"

# AlphaFold predicts the protein, not its cofactor: these models arrive with no metal,
# and a zinc metalloprotease without its zinc is the wrong receptor. Borrow the Zn from
# the BoNT/A light-chain crystal by superposing the template onto the model.
PAIRS = [
    ("BoNT/E", "MODEL_Q00496_BoNT_E_AF.pdb", "MODEL_BoNT_E_with_Zn.pdb"),
    ("BoNT/F", "MODEL_P30996_BoNT_F_AF.pdb", "MODEL_BoNT_F_with_Zn.pdb"),
]

if not os.path.exists(TEMPLATE):
    raise SystemExit("template {} not found".format(TEMPLATE))

for label, filename, outfile in PAIRS:
    path = os.path.join(BONT_MODELS, filename)
    if not os.path.exists(path):
        print("SKIP {}: {} not downloaded".format(label, filename))
        continue
    cmd.delete("all")
    cmd.load(path, "pred")
    cmd.load(TEMPLATE, "template")
    rms, atoms = cmd.super("template and polymer", "pred")[:2]
    cmd.create("with_zn", "pred or (template and resn ZN)")
    out = os.path.join(BONT_MODELS, outfile)
    cmd.save(out, "with_zn")
    zinc = cmd.count_atoms("with_zn and resn ZN")
    print("{}: template superposed at {:.2f} A over {} atoms, Zn atoms now {}  ->  {}/{}".format(
        label, rms, atoms, zinc, BONT_MODELS, outfile))
    # Confirm the borrowed metal landed in a real site. Note the two-step form: PyMOL
    # parses `byres (X) and name CA` as `byres (X and name CA)`, which selects nothing,
    # because no contact atom is itself a CA. Select first, then filter the selection.
    cmd.select("zn_shell", "byres (with_zn and polymer and not hydro within 2.6 of "
                           "(with_zn and resn ZN))")
    shell = []
    cmd.iterate("zn_shell and name CA", "shell.append((resn, resi))", space={"shell": shell})
    print("   Zn coordination shell: " + (", ".join(r + i for r, i in shell) or "none within 2.6 A"))
    if len(shell) < 3:
        print("   WARNING: fewer than three coordinating residues — check the placement")
        print("   in PyMOL before docking into this site.")
PY
    warn "A borrowed metal is a modelling assumption. Say so in any result that uses it."
}

# ── Dispatch ────────────────────────────────────────────────────────
preflight

case "$STEP" in
    session)     do_session ;;
    clean)       do_clean ;;
    check)       do_check ;;
    fetch)       do_fetch ;;
    plddt)       do_plddt ;;
    validate)    do_validate ;;
    identity)    do_identity ;;
    trim)        do_trim ;;
    zn-transfer) do_zn_transfer ;;
    all)
        do_session; do_clean; do_check; do_fetch
        do_plddt; do_validate; do_identity; do_trim; do_zn_transfer ;;
esac

step "Done"
echo "receptors + sessions : $RECEPTORS"
echo "models               : $DATASETS/*/MODEL_*.pdb"
echo "next                 : set COMMAND=control in run_docking.sh and run it"
