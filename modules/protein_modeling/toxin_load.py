# toxin_load.py — helper for the Ologist Workshop PyMOL Toxin Dataset
# Lives in modules/protein_modeling/. Run PyMOL from WORKING_FOLDER/INPUT_DATASETS/from_Database/, because the glob patterns below are relative to it.
# Usage: pymol -c -q -d "cd .../WORKING_FOLDER/INPUT_DATASETS/from_Database; run ../../../modules/protein_modeling/toxin_load.py; load_all_pdbs(); color_by_chain_all(); save toxin_dataset.pse"

from glob import glob
from os.path import basename
from pymol import util

def load_all_pdbs():
    """Load every .pdb file from all species folders into PyMOL."""
    pdb_files = sorted(glob("[1-9]_*/*.pdb"))
    if not pdb_files:
        print("WARNING: no .pdb files found in any species folder")
        return
    for path in pdb_files:
        name = basename(path).replace(".pdb", "")
        cmd.load(path, name)
        print(f"Loaded: {name}")
    print(f"Total loaded: {len(pdb_files)} PDB files")

def color_by_chain_all():
    """Color every loaded object by chain using distinct colors."""
    for obj in cmd.get_object_list():
        util.cbc(obj)

def load_target(pdb_filename, object_name):
    """Load a single PDB file by filename, searching all species folders."""
    matches = glob(f"[1-9]_*/{pdb_filename}")
    if not matches:
        print(f"ERROR: {pdb_filename} not found in any species folder")
        return
    path = matches[0]
    cmd.load(path, object_name)
    print(f"Loaded: {object_name} <- {path}")

def clean_receptor(object_name, outfile):
    """Strip solvent, ions, hydrogens, and cofactors from an object, then save as a clean PDB for docking.

    WARNING: this removes ALL ions, including a catalytic zinc. Do not use it on a
    BoNT light chain without putting the Zn back — see Step 6 of the manual.
    """
    cmd.remove(f"{object_name} and solvent")
    cmd.remove(f"{object_name} and ions")
    cmd.remove(f"{object_name} and hydro")
    cmd.remove(f"{object_name} and resn NAD")
    cmd.remove(f"{object_name} and resn SO4")
    cmd.remove(f"{object_name} and resn MAN")
    cmd.save(outfile, object_name)
    print(f"Cleaned receptor saved: {outfile}")


def box_from_selection(selection, padding=4.0):
    """Report the centre and a suggested size of a docking box around a selection.

    Returns (centre, size) as two 3-tuples. The centre is the midpoint of the
    selection's bounding box, which for an irregular pocket is a better box centre
    than the centre of mass.
    """
    heavy = f"({selection}) and not hydro"
    (min_x, min_y, min_z), (max_x, max_y, max_z) = cmd.get_extent(heavy)
    centre = (round((min_x + max_x) / 2, 2),
              round((min_y + max_y) / 2, 2),
              round((max_z + min_z) / 2, 2))
    size = (round(max_x - min_x + 2 * padding, 1),
            round(max_y - min_y + 2 * padding, 1),
            round(max_z - min_z + 2 * padding, 1))
    coords = cmd.get_coords(heavy)
    centroid = tuple(round(float(v), 2) for v in coords.mean(axis=0))
    print(f"selection      : {selection} ({len(coords)} heavy atoms)")
    print(f"box centre     : {centre[0]}, {centre[1]}, {centre[2]}   (bounding-box midpoint)")
    print(f"atom centroid  : {centroid[0]}, {centroid[1]}, {centroid[2]}   "
          "(matches the Step 10 recipe table)")
    print(f"suggested size : {size[0]} x {size[1]} x {size[2]} A "
          f"(span + {padding} A padding each side)")
    if max(size) > 30:
        print("NOTE: that span is large. If this selection is whole residues, build the box "
              "from\n      the contact atoms instead — drop the `byres` from your selection.")
    return centre, size


def write_box(selection, outfile, size=None, padding=4.0):
    """Write a Vina --config box file for a selection.

    `size` may be a single number for a cube, or None to use the selection's own
    span plus padding.
    """
    centre, suggested = box_from_selection(selection, padding)
    if size is None:
        box = suggested
    elif isinstance(size, (int, float)):
        box = (float(size),) * 3
    else:
        box = tuple(float(s) for s in size)
    with open(outfile, "w") as handle:
        handle.write(f"center_x = {centre[0]}\n")
        handle.write(f"center_y = {centre[1]}\n")
        handle.write(f"center_z = {centre[2]}\n")
        handle.write(f"size_x = {box[0]}\n")
        handle.write(f"size_y = {box[1]}\n")
        handle.write(f"size_z = {box[2]}\n")
    print(f"Vina box file written: {outfile}")
    return outfile
