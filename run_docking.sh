#!/usr/bin/env bash
# run_docking.sh — scripted Part III of the Protein Modeling Workshop
#
# Prepares receptors and search boxes, then runs AutoDock Vina, following
# docs/Drafts/Protein_Modeling_Workshop_Manual.md Steps 9–11.
#
# Driven entirely by the CONTROL PANEL variables below — this script takes no
# command-line arguments. Edit the variables, then run:  bash run_docking.sh
#
# The one habit that matters: before docking anything whose answer you do not know,
# redock a ligand whose answer you do know. COMMAND=control does that and records the
# result; COMMAND=dock refuses to run until it has, unless SKIP_CONTROL_CHECK=true.
#
# Needs the docking toolkit: vina, meeko (mk_prepare_receptor.py, mk_prepare_ligand.py,
# mk_export.py), pdbfixer, rdkit, pymol. See setup_conda_envs.sh / Step 7 of the manual.

set -euo pipefail

# ────────────────────────────────────────────────────────────────────
# CONTROL PANEL — set the run here; the script takes no arguments.
# ────────────────────────────────────────────────────────────────────

# ── What to run ─────────────────────────────────────────────────────
COMMAND=""              # list | control | cautionary | prep | dock | prep-all
TARGET=""               # target key for prep/dock (COMMAND=list shows them), e.g. bont_e

# ── Vina search settings — the knobs you are most likely to change ──
EXHAUSTIVENESS=16       # search effort: higher digs harder but runs slower
SEED=42                 # random seed: fixed so a run reproduces exactly
NUM_MODES=9             # how many poses Vina reports

# ── Ligand and box overrides ────────────────────────────────────────
USER_LIGAND=""          # ligand to dock (.sdf/.mol); empty = redock the crystal ligand
CENTER_OVERRIDE=""      # box centre "X Y Z"; required for the model-based (user) targets
BOX_OVERRIDE=""         # cube edge in angstrom; empty = the recipe's value

# ── Environment and paths ───────────────────────────────────────────
ENV_WANTED="${OLOGIST_ENV:-}"   # conda env name; OLOGIST_ENV wins, empty = autodetect
DATASETS=""             # dataset folder; empty = <script dir>/Datasets
OUT_DIR=""              # run folder; empty = <script dir>/work/docking

# ── Behavior toggles ────────────────────────────────────────────────
SKIP_CONTROL_CHECK=false  # true = dock even though the FimH control has not passed here
DRY_RUN=false             # true = print the commands, run nothing
# ────────────────────────────────────────────────────────────────────

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -z "$DATASETS" ]] && DATASETS="$ROOT/Datasets"
[[ -d "$DATASETS" ]] && DATASETS="$(cd "$DATASETS" && pwd)"
[[ -z "$OUT_DIR"  ]] && OUT_DIR="$ROOT/work/docking"

GREEN='\033[0;32m'  YELLOW='\033[0;33m'  RED='\033[0;31m'  BLUE='\033[0;34m'  NC='\033[0m'
info()  { printf "${GREEN}[OK]${NC}   %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC} %s\n" "$*"; }
fail()  { printf "${RED}[FAIL]${NC} %s\n" "$*" >&2; }
step()  { printf "\n${BLUE}━━━ %s ━━━${NC}\n" "$*"; }

# ── Per-target recipes (Step 10) ────────────────────────────────────
# Fields, tab-separated:
#   key | receptor file | pre-selection | split resn | keep metals | centre | box | redock resn | note
# pre-selection  what survives before the split (chain choice); "-" keeps everything
# split resn     residue whose centroid defines the box; "-" means use the centre field
# centre         explicit box centre, or "ligand" to measure it from the split residue,
#                or "interface" to compute it from the chain A/B contact surface,
#                or "user" when only you can decide (model-based targets)
RECIPES=$(cat <<'TSV'
fimh	3_Klebsiella_pneumoniae/PDB_9AT9_FimH_lectin_mannose.pdb	-	MAN	no	ligand	20	MAN	Validated control — redock mannose first
bonta_lc	2_Clostridium_botulinum/PDB_1XTG_BoNT_A_LC_SNAP25.pdb	(chain A and polymer) or resn ZN	ZN	yes	ligand	22	-	Catalytic Zn site; SNAP-25 chain B removed
bontb_cat	2_Clostridium_botulinum/PDB_1EPW_BoNT_B.pdb	polymer or resn ZN	ZN	yes	ligand	22	-	Catalytic Zn site; the blocking sulfates are dropped
bontb_rbd	2_Clostridium_botulinum/PDB_1I1E_BoNT_B_doxorubicin.pdb	-	DM2	yes	ligand	24	DM2	Receptor-binding site, 81 A from the Zn — cautionary control
exoa	4_Pseudomonas_aeruginosa/PDB_1AER_ExoA.pdb	(chain A and polymer) or resn TAD	TAD	no	ligand	24	TAD	NAD cofactor site, chain A; TAD is a redocking control
tnt	1_Mycobacterium_tuberculosis/PDB_4QLP_TNT_immunity.pdb	chain B and polymer	-	no	-3.75 20.87 -20.40	22	-	NAD cleft at the IFT interface; chain A removed
exot_secretion	4_Pseudomonas_aeruginosa/PDB_6JNP_ExoT_SpcS_complex.pdb	(chain A or chain B) and polymer	-	no	interface	24	-	ExoT-SpcS interface; catalytic domain absent from this crystal
bont_e	2_Clostridium_botulinum/MODEL_BoNT_E_with_Zn.pdb	-	ZN	yes	ligand	22	-	Model + transferred Zn — run run_protein_modeling.sh STEP=zn-transfer first
bont_f	2_Clostridium_botulinum/MODEL_BoNT_F_with_Zn.pdb	-	ZN	yes	ligand	22	-	Model + transferred Zn — run run_protein_modeling.sh STEP=zn-transfer first
exos_cat	4_Pseudomonas_aeruginosa/MODEL_ExoS_trimmed.pdb	-	-	no	user	24	-	Model, ADPRT domain 232-453 (pLDDT 89) — locate a pocket, set CENTER_OVERRIDE
exot_cat	4_Pseudomonas_aeruginosa/MODEL_ExoT_ADPRT_trimmed.pdb	-	-	no	user	24	-	Model, ADPRT domain 236-457 (pLDDT 87) — locate a pocket, set CENTER_OVERRIDE
exou	4_Pseudomonas_aeruginosa/MODEL_ExoU_PLA2_trimmed.pdb	-	-	no	user	24	-	Model, PLA2 dyad region (pLDDT 72, weakest) — set CENTER_OVERRIDE
mrkd	3_Klebsiella_pneumoniae/MODEL_MrkD_trimmed.pdb	-	-	no	user	22	-	Model, pLDDT 89, no known ligand — find a pocket first, set CENTER_OVERRIDE
TSV
)

# ── Validate the control-panel settings ─────────────────────────────
case "$COMMAND" in
    list|control|cautionary|prep|dock|prep-all) ;;
    "") fail "COMMAND is not set. Edit the CONTROL PANEL at the top of this script."
        fail "Choose one of: list control cautionary prep dock prep-all"; exit 1 ;;
    *)  fail "unknown COMMAND: '$COMMAND'"
        fail "Choose one of: list control cautionary prep dock prep-all"; exit 1 ;;
esac

if [[ "$COMMAND" == "prep" || "$COMMAND" == "dock" ]] && [[ -z "$TARGET" ]]; then
    fail "COMMAND=$COMMAND needs TARGET set in the CONTROL PANEL (COMMAND=list shows the keys)."
    exit 1
fi

CONTROL_STAMP="$OUT_DIR/.control_passed"
TMP=""
cleanup() { [[ -n "$TMP" && -d "$TMP" ]] && rm -rf "$TMP"; }
trap cleanup EXIT

# ── Recipe lookup ───────────────────────────────────────────────────
recipe_line() { printf '%s\n' "$RECIPES" | awk -F'\t' -v k="$1" '$1==k {print; exit}'; }

load_recipe() {
    local line; line="$(recipe_line "$1")"
    [[ -z "$line" ]] && { fail "unknown target: $1"; echo; cmd_list; exit 1; }
    IFS=$'\t' read -r R_KEY R_PDB R_PRE R_SPLIT R_METALS R_CENTER R_BOX R_REDOCK R_NOTE <<<"$line"
    if [[ -n "$BOX_OVERRIDE" ]]; then R_BOX="$BOX_OVERRIDE"; fi
    if [[ -n "$CENTER_OVERRIDE" ]]; then R_CENTER="$CENTER_OVERRIDE"; fi
    # An explicit success: a function whose last statement is `[[ ... ]] && x` returns 1
    # when the test is false, and under `set -e` that aborts the caller.
    return 0
}

# Residue names that are cofactors, not ligands to be docked: a box may be centred on
# one, but it has to stay in the receptor.
is_metal() {
    case "$1" in
        ZN|MG|MN|FE|CU|NI|CO|CA) return 0 ;;
        *) return 1 ;;
    esac
}

cmd_list() {
    step "Per-target docking recipes (Step 10)"
    printf "%-15s %-42s %-6s %-9s %s\n" "target" "receptor" "box" "centre" "note"
    printf '%s\n' "$RECIPES" | awk -F'\t' '{
        centre = ($6 == "ligand" ? "measured" : ($6 == "user" ? "you pass" : ($6 == "interface" ? "computed" : "fixed")))
        printf "%-15s %-42s %-6s %-9s %s\n", $1, $2, $7 " A", centre, $9
    }'
    echo
    echo "Targets marked 'you pass' need CENTER_OVERRIDE=\"X Y Z\": they are predicted models"
    echo "with no ligand to measure a site from. Find a pocket in PyMOL first, set the centre."
    echo "Targets whose receptor is a MODEL_* file need run_protein_modeling.sh STEP=trim (and"
    echo "STEP=zn-transfer for BoNT/E and /F) to have been run."
    echo
    echo "Known snag: bontb_cat (1EPW) has a truncated TYR830 that meeko cannot template"
    echo "after rebuild. It is far from the catalytic Zn, so it is safe to drop — but the"
    echo "bundled dock_prep.py does not expose meeko's --allow_bad_res, so this receptor"
    echo "needs a manual meeko run. Every other crystal target prepares cleanly."
}

# ── Environment ─────────────────────────────────────────────────────
DOCK_TOOLS=(vina pymol python mk_prepare_receptor.py mk_prepare_ligand.py mk_export.py)

tools_present() {
    local dir="$1" tool
    for tool in "${DOCK_TOOLS[@]}"; do
        if [[ -n "$dir" ]]; then [[ -x "$dir/$tool" ]] || return 1
        else command -v "$tool" &>/dev/null || return 1; fi
    done
}

resolve_env() {
    if [[ -z "$ENV_WANTED" ]] && tools_present ""; then
        info "using the docking tools already on PATH ($(command -v vina))"
        return
    fi
    local candidates=() name prefix
    [[ -n "$ENV_WANTED" ]] && candidates+=("$ENV_WANTED")
    candidates+=(dock-workshop protein_modeling protein_model)
    for name in "${candidates[@]}"; do
        prefix="$(conda info --envs 2>/dev/null | awk -v n="$name" '$1==n {print $NF}')" || true
        [[ -z "$prefix" ]] && continue
        if tools_present "$prefix/bin"; then
            # Prepend to PATH, not just call by absolute path: dock_prep.py shells out to
            # the meeko scripts, so they have to be findable by the child process too.
            export PATH="$prefix/bin:$PATH"
            info "using conda environment '$name' ($prefix)"
            return
        fi
        warn "environment '$name' is missing one of: ${DOCK_TOOLS[*]}"
    done
    fail "no environment with the full docking toolkit found."
    fail "Create one:  conda create -n dock-workshop -c conda-forge python=3.11 \\"
    fail "                 pymol-open-source vina meeko openbabel pdbfixer"
    exit 1
}

preflight() {
    [[ -f "$DATASETS/dock_prep.py" ]] || { fail "dock_prep.py not found in $DATASETS"; exit 1; }
    resolve_env
    TMP="$(mktemp -d)"
    $DRY_RUN || mkdir -p "$OUT_DIR"
}

dock_prep() {
    if $DRY_RUN; then printf '  [dry-run] python dock_prep.py %s\n' "$*"; return 0; fi
    ( cd "$DATASETS" && python dock_prep.py "$@" )
}

run_pymol() {
    local script="$TMP/$1.py"; shift
    cat > "$script"
    if $DRY_RUN; then printf '  [dry-run] pymol -cq %s\n' "$script"; return 0; fi
    ( cd "$DATASETS" && pymol -cq "$script" )
}

# ── Receptor preparation ────────────────────────────────────────────
# Returns through the globals PREPARED_RECEPTOR (basename of the PDBQT/box) and
# REF_LIGAND (the crystal ligand PDB, empty when there is none).
prepare_target() {
    load_recipe "$1"
    local work="$OUT_DIR/$R_KEY"
    $DRY_RUN || mkdir -p "$work"

    step "$R_KEY — $R_NOTE"

    local source_pdb="$DATASETS/$R_PDB"
    if [[ ! -f "$source_pdb" ]] && ! $DRY_RUN; then
        fail "receptor not found: $R_PDB"
        case "$R_PDB" in
            */MODEL_BoNT_*_with_Zn.pdb) fail "run run_protein_modeling.sh with STEP=zn-transfer" ;;
            */MODEL_*)                    fail "run run_protein_modeling.sh with STEP=fetch then STEP=trim" ;;
        esac
        exit 1
    fi

    # Step 1 — restrict to the chains that make up the receptor, if the recipe says so.
    local split_input="$R_PDB"
    if [[ "$R_PRE" != "-" ]]; then
        OW_IN="$R_PDB" OW_OUT="$work/receptor_subset.pdb" OW_KEEP="$R_PRE" \
        run_pymol subset <<'PY'
import os
from pymol import cmd
cmd.load(os.environ["OW_IN"], "src")
before = cmd.count_atoms("src")
cmd.create("sub", "src and ({})".format(os.environ["OW_KEEP"]))
cmd.remove("sub and solvent")
cmd.save(os.environ["OW_OUT"], "sub")
print("subset: {} -> {} atoms kept by `{}`".format(
    before, cmd.count_atoms("sub"), os.environ["OW_KEEP"]))
PY
        split_input="$work/receptor_subset.pdb"
    fi

    # Step 2 — split off the site residue (ligand or metal) and audit the receptor.
    REF_LIGAND=""
    local receptor_pdb
    if [[ "$R_SPLIT" != "-" ]]; then
        local metal_flag=()
        [[ "$R_METALS" == "yes" ]] && metal_flag=(--keep-metals)
        dock_prep split "$split_input" "$R_SPLIT" --out "$work" "${metal_flag[@]}"
        receptor_pdb="$work/receptor.pdb"
        REF_LIGAND="$work/ligand_ref.pdb"

        # `dock_prep.py split` builds its receptor as "everything kept AND NOT the split
        # residue", so splitting on the catalytic metal strips that metal out of the
        # receptor: --keep-metals only protects the *other* metals. For a metal site the
        # split is therefore used only to measure the box centre, and the receptor is
        # rebuilt from the subset with the metal still in place. A zinc metalloprotease
        # without its zinc is not the enzyme you meant to dock into.
        if is_metal "$R_SPLIT"; then
            OW_IN="$split_input" OW_OUT="$work/receptor_with_metal.pdb" OW_METAL="$R_SPLIT" \
            run_pymol keep_metal <<'PYSCRIPT'
import os
from pymol import cmd
metal = os.environ["OW_METAL"]
cmd.load(os.environ["OW_IN"], "rec")
cmd.remove("rec and solvent")
cmd.remove("rec and hydro")
count = cmd.count_atoms("rec and resn {}".format(metal))
cmd.save(os.environ["OW_OUT"], "rec")
print("receptor keeps its cofactor: {} x {} ({} heavy atoms total)".format(
    count, metal, cmd.count_atoms("rec")))
if count == 0:
    raise SystemExit("ERROR: no {} left in the receptor - check the pre-selection".format(metal))
PYSCRIPT
            receptor_pdb="$work/receptor_with_metal.pdb"
        fi
    else
        OW_IN="$split_input" OW_OUT="$work/receptor.pdb" \
        run_pymol strip <<'PY'
import os
from pymol import cmd
cmd.load(os.environ["OW_IN"], "rec")
cmd.remove("rec and solvent")
cmd.remove("rec and hydro")
cmd.save(os.environ["OW_OUT"], "rec")
print("receptor: {} heavy atoms".format(cmd.count_atoms("rec")))
PY
        receptor_pdb="$work/receptor.pdb"
    fi

    # Step 3 — decide the box centre.
    local centre_args=() cx cy cz
    case "$R_CENTER" in
        ligand)
            centre_args=(--ref-ligand "$REF_LIGAND") ;;
        interface)
            local centre
            centre="$(OW_IN="$split_input" run_pymol interface <<'PY'
import os
from pymol import cmd
cmd.load(os.environ["OW_IN"], "src")
# The contact surface between the effector and its chaperone, not either chain's
# own centre of mass — a box on the whole complex would search mostly solvent.
cmd.select("contacts", "(src and chain A and polymer and not hydro) within 4.5 of "
                       "(src and chain B and polymer and not hydro)")
n = cmd.count_atoms("contacts")
if n == 0:
    raise SystemExit("ERROR: no A/B interface contacts found")
# The centroid of the contact atoms, matching what toxin_load.box_from_selection
# reports as "atom centroid" — not the mass-weighted centre.
coords = cmd.get_coords("contacts")
x, y, z = (float(v) for v in coords.mean(axis=0))
print("interface contact atoms: {}".format(n))
print("CENTRE {:.2f} {:.2f} {:.2f}".format(x, y, z))
PY
)"
            if $DRY_RUN; then
                centre_args=(--center-x 0 --center-y 0 --center-z 0)
            else
                printf '%s\n' "$centre" | grep -v '^CENTRE' || true
                read -r cx cy cz <<<"$(printf '%s\n' "$centre" | awk '/^CENTRE/ {print $2, $3, $4}')"
                [[ -z "${cz:-}" ]] && { fail "could not compute the interface centre"; exit 1; }
                info "interface centre: $cx $cy $cz"
                centre_args=(--center-x "$cx" --center-y "$cy" --center-z "$cz")
            fi ;;
        user)
            fail "$R_KEY is a predicted model with no ligand to measure a site from."
            fail "Locate a pocket, then set in the CONTROL PANEL:  CENTER_OVERRIDE=\"X Y Z\""
            fail "In PyMOL:  run toxin_load.py  then  box_from_selection(\"<your selection>\")"
            exit 1 ;;
        *)
            read -r cx cy cz <<<"$R_CENTER"
            [[ -z "${cz:-}" ]] && { fail "malformed centre for $R_KEY: '$R_CENTER'"; exit 1; }
            centre_args=(--center-x "$cx" --center-y "$cy" --center-z "$cz") ;;
    esac

    # Step 4 — rebuild missing side chains, then write the receptor PDBQT and the box.
    dock_prep receptor "$receptor_pdb" --out "$work/$R_KEY" \
        "${centre_args[@]}" --box-size "$R_BOX"

    if is_metal "$R_SPLIT" && ! $DRY_RUN; then
        local metal_lines
        metal_lines="$(grep -c "$R_SPLIT" "$work/$R_KEY.pdbqt" 2>/dev/null || true)"
        if [[ "${metal_lines:-0}" -gt 0 ]]; then
            info "$R_SPLIT survived into the prepared receptor ($metal_lines line(s))"
        else
            fail "$R_SPLIT is not in $R_KEY.pdbqt — the pocket is missing its cofactor."
            fail "Do not dock into this site until you know why."
            exit 1
        fi
    fi
    info "receptor + box ready: $work/$R_KEY.pdbqt, $R_KEY.box.txt"
    echo "     inspect the box:  pymol -q $work/receptor.pdb $work/$R_KEY.box.pdb"
}

# ── Docking ─────────────────────────────────────────────────────────
run_vina() {
    local receptor="$1" ligand="$2" work="$3" out="$3/poses.pdbqt"
    step "Vina — exhaustiveness $EXHAUSTIVENESS, seed $SEED, $NUM_MODES modes"
    if $DRY_RUN; then
        printf '  [dry-run] vina --receptor %s.pdbqt --ligand %s --config %s.box.txt --out %s\n' \
            "$receptor" "$ligand" "$receptor" "$out"
        return 0
    fi
    vina --receptor "$receptor.pdbqt" \
         --ligand "$ligand" \
         --config "$receptor.box.txt" \
         --exhaustiveness "$EXHAUSTIVENESS" --seed "$SEED" --num_modes "$NUM_MODES" \
         --out "$out" | tee "$work/vina.log"
    mk_export.py "$out" -s "$work/poses.sdf"
    info "poses: $work/poses.pdbqt and $work/poses.sdf"
}

# Redock the crystal ligand and report whether the protocol reproduced a known answer.
redock_and_score() {
    local key="$1"
    local work="$OUT_DIR/$key"
    load_recipe "$key"
    prepare_target "$key"

    if [[ "$R_REDOCK" == "-" ]]; then
        fail "$key has no crystal ligand to redock"; exit 1
    fi

    step "$key — ligand chemistry from the RCSB chemical component dictionary"
    dock_prep ligand "$work/ligand_ref.pdb" "$R_REDOCK" --out "$work/ligand.pdbqt"

    run_vina "$work/$key" "$work/ligand.pdbqt" "$work"

    step "$key — RMSD of the poses against the crystal answer"
    if $DRY_RUN; then
        printf '  [dry-run] python dock_prep.py rmsd %s/poses.sdf %s/ligand_correct.sdf\n' "$work" "$work"
        return 0
    fi
    dock_prep rmsd "$work/poses.sdf" "$work/ligand_correct.sdf" --top 5 | tee "$work/rmsd.log"
    grep -q '^PASS' "$work/rmsd.log"
}

cmd_control() {
    step "Step 9 — validate the protocol on a known answer (mannose into FimH)"
    echo "9AT9 holds mannose in its pocket at 1.34 A. Strip it, dock it back, and measure"
    echo "how close it lands. If this fails, no score from this setup means anything."
    if redock_and_score fimh; then
        $DRY_RUN || { mkdir -p "$OUT_DIR"; date -u +"passed %Y-%m-%dT%H:%M:%SZ" > "$CONTROL_STAMP"; }
        info "PASS — the protocol reproduces the crystal pose. Scores are worth interpreting."
    else
        $DRY_RUN || rm -f "$CONTROL_STAMP"
        fail "FAIL — the top pose is more than 2.0 A from the crystal ligand."
        fail "Do not interpret scores from this setup. Work through Step 11 of the manual"
        fail "before changing anything: box too large, wrong ligand chemistry, a deleted"
        fail "pocket side chain and a missing cofactor all produce exactly this."
        exit 1
    fi
}

cmd_cautionary() {
    step "Step 11 — the cautionary control (doxorubicin into BoNT/B)"
    echo "Doxorubicin is described as a BoNT/B inhibitor, and 1I1E holds its crystal pose,"
    echo "so redocking it looks like a second control. It is expected to FAIL: the ligand"
    echo "sits 81 A from the catalytic Zn in the receptor-binding domain, and a confident"
    echo "score on the wrong pose is the failure mode this exercise exists to show."
    if redock_and_score bontb_rbd; then
        warn "This one PASSED, which the manual does not expect — read $OUT_DIR/bontb_rbd/rmsd.log"
        warn "and check which site the box actually covered before drawing any conclusion."
    else
        info "FAIL, as expected — a good-looking affinity on a pose 4+ A out."
        info "Read $OUT_DIR/bontb_rbd/rmsd.log alongside Step 11."
    fi
}

require_control() {
    $SKIP_CONTROL_CHECK && return 0
    $DRY_RUN && return 0
    [[ -f "$CONTROL_STAMP" ]] && return 0
    fail "the FimH control has not passed in this run folder."
    fail "Set COMMAND=control in the CONTROL PANEL and run this script."
    fail "Or set SKIP_CONTROL_CHECK=true if you validated the protocol elsewhere."
    exit 1
}

cmd_dock() {
    [[ -z "$TARGET" ]] && { fail "dock needs TARGET set — set COMMAND=list to see the keys"; exit 1; }
    load_recipe "$TARGET"

    if [[ -n "$USER_LIGAND" ]]; then
        [[ -f "$USER_LIGAND" ]] || { fail "ligand file not found: $USER_LIGAND"; exit 1; }
        require_control
        prepare_target "$TARGET"
        local work="$OUT_DIR/$TARGET"
        step "$TARGET — preparing $USER_LIGAND"
        if $DRY_RUN; then
            printf '  [dry-run] mk_prepare_ligand.py -i %s -o %s/ligand.pdbqt\n' "$USER_LIGAND" "$work"
        else
            mk_prepare_ligand.py -i "$USER_LIGAND" -o "$work/ligand.pdbqt"
        fi
        run_vina "$work/$TARGET" "$work/ligand.pdbqt" "$work"
        step "Result"
        echo "No RMSD is reported: there is no crystal answer for this ligand. The score is"
        echo "only as trustworthy as the control run — see $work/vina.log."
        return
    fi

    if [[ "$R_REDOCK" != "-" ]]; then
        require_control
        redock_and_score "$TARGET" || warn "the redock did not reproduce the crystal pose — see $OUT_DIR/$TARGET/rmsd.log"
        return
    fi

    # No crystal ligand and none supplied: prepare the site and say so plainly.
    prepare_target "$TARGET"
    step "Result"
    echo "$TARGET has no crystal ligand to redock, so nothing was docked."
    echo "Supply a compound by setting in the CONTROL PANEL:"
    echo "  COMMAND=dock  TARGET=$TARGET  USER_LIGAND=my_compound.sdf"
    echo
    echo "To build a ligand from SMILES first:"
    echo "  python - <<'EOF'"
    echo "  from rdkit import Chem"
    echo "  from rdkit.Chem import AllChem"
    echo "  mol = Chem.AddHs(Chem.MolFromSmiles(\"<your SMILES>\"))"
    echo "  AllChem.EmbedMolecule(mol, randomSeed=42)"
    echo "  AllChem.MMFFOptimizeMolecule(mol)"
    echo "  Chem.MolToMolFile(mol, \"my_compound.sdf\")"
    echo "  EOF"
}

cmd_prep_all() {
    # Collect the keys first: the prepared targets run PyMOL and Vina, which must not
    # inherit the recipe list on stdin.
    local keys=() key centre pdb
    while IFS=$'\t' read -r key _ _ _ _ centre _; do
        [[ -z "$key" ]] && continue
        if [[ "$centre" == "user" ]]; then
            warn "skipping $key — a predicted model with no measurable site; set CENTER_OVERRIDE"
            continue
        fi
        keys+=("$key")
    done <<<"$RECIPES"

    local failed=()
    for key in "${keys[@]}"; do
        pdb="$(recipe_line "$key" | cut -f2)"
        if [[ ! -f "$DATASETS/$pdb" ]]; then
            warn "skipping $key — receptor file missing: $pdb"
            continue
        fi
        # Isolate each target in a subshell so one crystal meeko cannot template
        # (e.g. bontb_cat / 1EPW) does not abort the whole batch under `set -e`.
        if ( prepare_target "$key" </dev/null ); then
            :
        else
            failed+=("$key")
            warn "prep failed for $key — see the output above"
        fi
    done
    if [[ ${#failed[@]} -gt 0 ]]; then
        warn "targets needing manual attention: ${failed[*]}"
    fi
}

# ── Dispatch ────────────────────────────────────────────────────────
case "$COMMAND" in
    list)       cmd_list; exit 0 ;;
    control)    preflight; cmd_control ;;
    cautionary) preflight; cmd_cautionary ;;
    prep)       [[ -z "$TARGET" ]] && { fail "prep needs a target"; exit 1; }
                preflight; prepare_target "$TARGET" ;;
    dock)       preflight; cmd_dock ;;
    prep-all)   preflight; cmd_prep_all ;;
esac

step "Done"
echo "run folder : $OUT_DIR"
echo "manual     : docs/Drafts/Protein_Modeling_Workshop_Manual.md (Steps 9-11)"
