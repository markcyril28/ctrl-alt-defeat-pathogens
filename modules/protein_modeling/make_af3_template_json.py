#!/usr/bin/env python3
"""Generate an AlphaFold3 local-dialect JSON with structural templates.

Reads the AF3 server batch JSON program E wrote (metagenomics query
sequences), maps each entry to its matching reference protein from
WORKING_FOLDER/INPUT_DATASETS/from_Database, converts the template PDB to mmCIF, aligns
query→template to build index mappings, and writes an alphafold3-dialect JSON
ready for a local run.

One batch per reference organism, because that is how program E writes them,
so this writes one local JSON and one folder of converted templates per
organism as well:

    E_AlphaFold3_Inputs/alphafold3_local/<organism>/AF3_local_batch_with_templates.json
    E_AlphaFold3_Inputs/alphafold3_local/templates/<organism>/*.cif

Usage:
    python3 modules/protein_modeling/make_af3_template_json.py
"""

import glob
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser, MMCIFIO, Select

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(BASE, "..", ".."))

E_INPUTS = os.path.join(
    ROOT,
    "WORKING_FOLDER/RESULTS/protein_modeling/for_metagenomics_dataset/E_AlphaFold3_Inputs",
)

SERVER_DIR = os.path.join(E_INPUTS, "alphafoldserver")
LOCAL_DIR = os.path.join(E_INPUTS, "alphafold3_local")

# The batch files, one folder per organism. The batch is the input rather than
# the one-job-per-target files beside it because the batch is the thing you
# actually uploaded, so a local run with templates asks about the same set.
BATCH_GLOB = os.path.join(SERVER_DIR, "*", "*_batch_*.json")

OUT_NAME = "AF3_local_batch_with_templates.json"

DB = os.path.join(ROOT, "WORKING_FOLDER/INPUT_DATASETS/from_Database")

TEMPLATE_MAP = {
    "TNT_CpnT": {
        "pdb": os.path.join(
            DB, "1_Mycobacterium_tuberculosis/MODEL_O05442_CpnT_full_AF.pdb"
        ),
        "fasta": os.path.join(
            DB, "1_Mycobacterium_tuberculosis/FASTA_O05442_CpnT.fasta"
        ),
        "chain": "A",
        "label": "CpnT (O05442) AlphaFold model",
    },
    "FimH": {
        "pdb": os.path.join(
            DB,
            "3_Klebsiella_pneumoniae/MODEL_A0A0H3H2I8_FimH_full_AF.pdb",
        ),
        "fasta": os.path.join(
            DB,
            "3_Klebsiella_pneumoniae/"
            "FASTA_A0A0H3H2I8_FimH_Kpneumoniae_ST11.fasta",
        ),
        "chain": "A",
        "label": "FimH (A0A0H3H2I8) AlphaFold model",
    },
    "MrkD": {
        "pdb": os.path.join(
            DB, "3_Klebsiella_pneumoniae/MODEL_P21648_MrkD_AF.pdb"
        ),
        "fasta": os.path.join(
            DB, "3_Klebsiella_pneumoniae/FASTA_P21648_MrkD.fasta"
        ),
        "chain": "A",
        "label": "MrkD (P21648) AlphaFold model",
    },
    "ExoA_toxA": {
        "pdb": os.path.join(
            DB,
            "4_Pseudomonas_aeruginosa/PDB_1IKQ_ExoA_wildtype.pdb",
        ),
        "fasta": os.path.join(
            DB, "4_Pseudomonas_aeruginosa/FASTA_P11439_ExoA.fasta"
        ),
        "chain": "A",
        "label": "ExoA (P11439) PDB 1IKQ wildtype",
    },
    "ExoT": {
        "pdb": os.path.join(
            DB, "4_Pseudomonas_aeruginosa/MODEL_Q9I788_ExoT_AF.pdb"
        ),
        "fasta": os.path.join(
            DB, "4_Pseudomonas_aeruginosa/FASTA_Q9I788_ExoT_PAE.fasta"
        ),
        "chain": "A",
        "label": "ExoT (Q9I788) AlphaFold model",
    },
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}


class SingleChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] == " "


def pdb_chain_sequence(structure, chain_id):
    for chain in structure.get_chains():
        if chain.id == chain_id:
            residues = [r for r in chain.get_residues() if r.id[0] == " "]
            return "".join(
                THREE_TO_ONE.get(r.resname, "X") for r in residues
            )
    return ""


def read_fasta_seq(path):
    for record in SeqIO.parse(path, "fasta"):
        return str(record.seq)
    return ""


def align_and_map(query_seq, template_seq):
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1

    alignments = aligner.align(query_seq, template_seq)
    if not alignments:
        return [], []

    best = alignments[0]
    q_indices = []
    t_indices = []

    aligned = best.aligned
    for (q_start, q_end), (t_start, t_end) in zip(aligned[0], aligned[1]):
        for q_pos, t_pos in zip(range(q_start, q_end), range(t_start, t_end)):
            q_indices.append(q_pos)
            t_indices.append(t_pos)

    return q_indices, t_indices


def convert_pdb_to_cif(pdb_path, cif_path, chain_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("tmpl", pdb_path)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(cif_path, select=SingleChainSelect(chain_id))
    return cif_path


def match_template(job_name):
    for prefix, info in TEMPLATE_MAP.items():
        if job_name.startswith(prefix):
            return prefix, info
    return None, None


def main():
    batches = sorted(glob.glob(BATCH_GLOB))
    if not batches:
        print(f"no server batch files under {SERVER_DIR}/*/", file=sys.stderr)
        print("Run:  bash E_alphafold3_inputs.sh", file=sys.stderr)
        return 1

    for batch_path in batches:
        # The organism is the folder the batch sits in, which is the same
        # folder name program C grouped the sequences into.
        species = os.path.basename(os.path.dirname(batch_path))
        write_batch(species, batch_path)
    return 0


def write_batch(species, batch_path):
    template_out = os.path.join(LOCAL_DIR, "templates", species)
    out_json = os.path.join(LOCAL_DIR, species, OUT_NAME)

    with open(batch_path) as fh:
        batch = json.load(fh)

    os.makedirs(template_out, exist_ok=True)
    os.makedirs(os.path.dirname(out_json), exist_ok=True)

    cif_cache = {}
    struct_seq_cache = {}
    parser = PDBParser(QUIET=True)

    jobs = []
    summary = []

    for entry in batch:
        name = entry["name"]
        query_seq = entry["sequences"][0]["proteinChain"]["sequence"]

        prefix, tmpl_info = match_template(name)
        if tmpl_info is None:
            print(f"WARNING: no template mapping for {name}", file=sys.stderr)
            jobs.append({
                "name": name,
                "modelSeeds": [1],
                "sequences": [
                    {"protein": {"id": "A", "sequence": query_seq}}
                ],
                "dialect": "alphafold3",
                "version": 2,
            })
            continue

        cif_name = os.path.basename(tmpl_info["pdb"]).replace(".pdb", ".cif")
        cif_path = os.path.join(template_out, cif_name)

        if prefix not in cif_cache:
            convert_pdb_to_cif(
                tmpl_info["pdb"], cif_path, tmpl_info["chain"]
            )
            cif_cache[prefix] = cif_path

            structure = parser.get_structure(prefix, tmpl_info["pdb"])
            struct_seq_cache[prefix] = pdb_chain_sequence(
                structure, tmpl_info["chain"]
            )

        cif_rel = os.path.relpath(cif_cache[prefix], os.path.dirname(out_json))
        struct_seq = struct_seq_cache[prefix]

        q_idx, t_idx = align_and_map(query_seq, struct_seq)

        if not q_idx:
            print(
                f"WARNING: no alignment for {name} vs {tmpl_info['label']}",
                file=sys.stderr,
            )

        job = {
            "name": name,
            "modelSeeds": [1],
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": query_seq,
                        "templates": [
                            {
                                "mmcifPath": cif_rel,
                                "queryIndices": q_idx,
                                "templateIndices": t_idx,
                            }
                        ],
                    }
                }
            ],
            "dialect": "alphafold3",
            "version": 2,
        }
        jobs.append(job)

        ref_seq = read_fasta_seq(tmpl_info["fasta"])
        summary.append(
            f"  {name}: {len(query_seq)} aa query, "
            f"{len(struct_seq)} aa template ({tmpl_info['label']}), "
            f"{len(q_idx)} aligned positions"
        )

    with open(out_json, "w") as fh:
        json.dump(jobs, fh, indent=2)
        fh.write("\n")

    print(f"{species}: wrote {len(jobs)} jobs to {out_json}")
    print(f"  templates converted to mmCIF in {template_out}/")
    for line in summary:
        print(line)
    print()


if __name__ == "__main__":
    sys.exit(main())
