#!/usr/bin/env python3
"""Render one structure to a PNG, styled by what kind of file it is.

    render_structure.py <structure> <out_png> <width> <height>

Prints one TSV row: image, structure, source, chains, residues, ligands,
metals, style.

A figure of a protein is an argument about what matters in it, so the two
kinds of file in WORKING_FOLDER/INPUT_DATASETS/from_Database are not drawn the same way.

An experimental structure is a measurement. What it has to show is
architecture: how many chains, how they pack, and what crystallised in the
site. A complex is coloured one colour per chain. A single chain is coloured
along its sequence instead, blue at the N terminus to red at the C, because
colouring one chain by chain paints it one flat colour and hides the domain
layout that is the only architecture such a file has. Either way the crystal
ligand is drawn as sticks and a metal as a sphere, because a zinc rendered
as a cartoon is invisible and a zinc is the whole point of a
metalloprotease.

A prediction has no measured ligand and nothing to pack against. Its
B-factor column is not a B-factor: it is per-residue confidence, and that is
the only thing in the file worth looking at. So a MODEL_ file is coloured
by confidence on the AlphaFold palette — dark blue very high, light blue
confident, yellow low, orange very low — which turns the figure into the
one check program F keeps asking for. The confident core and the ragged
low-confidence loops are visible at a glance, and a model that is a good
fold with a bad active-site loop cannot hide behind a good mean.

The scale is detected rather than assumed. AlphaFold writes pLDDT from 0 to
100; SWISS-MODEL usually writes QMEANDisCo from 0 to 1. Thresholds meant for
one read as nonsense on the other, so the bands are rescaled when the values
look like they are on the 0-1 scale.
"""
import os
import sys

from pymol import cmd

METALS = {"ZN", "MG", "MN", "CA", "FE", "CU", "NI", "CO", "NA", "K", "FE2", "FES"}
SOLVENT = {"HOH", "WAT", "DOD"}

# AlphaFold's own bands, as (lower bound on the 0-100 scale, colour).
PLDDT_BANDS = [
    (90, "0x0053D6"),   # very high
    (70, "0x65CBF3"),   # confident
    (50, "0xFFDB13"),   # low
    (0,  "0xFF7D45"),   # very low
]

# Chain colours, in order. Deliberately not PyMOL's util.cbc default, whose
# first chain is green — the same green this file uses for a ligand.
CHAIN_COLOURS = ["skyblue", "salmon", "palegreen", "lightpink",
                 "paleyellow", "lightblue", "wheat", "violetpurple"]


def _is_prediction(stem):
    """True for a model, by the naming convention the dataset and program F set."""
    upper = stem.upper()
    return (upper.startswith("MODEL_") or upper.startswith("AF3")
            or upper.startswith("SWISSMODEL"))


def _residue_names(selection):
    """The distinct residue names in a selection, as a set."""
    names = set()
    cmd.iterate(selection, "names.add(resn)", space={"names": names})
    return names


def _colour_by_confidence(obj):
    """Paint the polymer on the AlphaFold palette, rescaling if needed.

    Returns the label for the style column, so the caller can report which
    scale was detected rather than leaving it implied.
    """
    values = []
    cmd.iterate("{} and polymer and name CA".format(obj),
                "values.append(b)", space={"values": values})
    if not values:
        return "none"

    # QMEANDisCo runs 0 to 1, pLDDT 0 to 100. Nothing in the file says which,
    # so the range decides, and a model with no variation at all is left
    # alone rather than being stretched into bands that mean nothing.
    top = max(values)
    if top <= 1.5:
        factor, scale = 0.01, "qmean_0_1"
    else:
        factor, scale = 1.0, "plddt_0_100"

    # Lowest band first. Each band repaints everything above its own floor,
    # so painting them in descending order would end with "b > 0" covering
    # the whole chain in the very-low colour and every model would look bad.
    for low, colour in reversed(PLDDT_BANDS):
        cmd.color(colour, "{} and polymer and b > {}".format(obj, low * factor))
    return "confidence:" + scale


def _colour_by_chain(obj, chains):
    """One colour per chain, so a complex reads as a complex."""
    for index, chain in enumerate(chains):
        colour = CHAIN_COLOURS[index % len(CHAIN_COLOURS)]
        cmd.color(colour, "{} and polymer and chain {}".format(obj, chain))
    return "chain"


def _colour_n_to_c(obj):
    """Rainbow along the sequence, blue at the N terminus to red at the C.

    Colouring a one-chain file by chain paints it a single flat colour and
    says nothing. These are the structures whose point is domain
    architecture — BoNT/A holotoxin is protease, then translocation, then
    binding, along the chain — and the rainbow is what makes those domains
    separate into blue, green and red instead of one blue blob.
    """
    cmd.spectrum("count", "rainbow", "{} and polymer and name CA".format(obj))
    return "n_to_c"


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)

    structure, out_png = sys.argv[1], sys.argv[2]
    width, height = int(sys.argv[3]), int(sys.argv[4])
    stem = os.path.splitext(os.path.basename(structure))[0]

    cmd.reinitialize()
    cmd.load(structure, "mol")

    # Waters clutter a figure and tell you nothing at this scale. Hydrogens
    # are usually absent from a crystal structure and modelled in a
    # prediction, so dropping them keeps the two kinds of file comparable.
    cmd.remove("mol and solvent")
    cmd.remove("mol and hydro")

    chains = cmd.get_chains("mol and polymer")
    residues = cmd.count_atoms("mol and polymer and name CA")
    het = _residue_names("mol and not polymer")
    metals = sorted(het & METALS)
    ligands = sorted(het - METALS - SOLVENT)

    # ---- geometry: one cartoon, ligands as sticks, metals as spheres ----
    cmd.hide("everything")
    cmd.show("cartoon", "mol and polymer")
    cmd.show("sticks", "mol and organic")
    cmd.show("spheres", "mol and inorganic")

    cmd.set("cartoon_transparency", 0.0)
    cmd.set("sphere_scale", 0.45, "mol and inorganic")
    cmd.set("stick_radius", 0.20)
    cmd.set("cartoon_smooth_loops", 1)

    # ---- colour: by confidence for a model, by chain for a measurement ----
    if _is_prediction(stem):
        style = _colour_by_confidence("mol")
    elif len(chains) > 1:
        style = _colour_by_chain("mol", chains)
    else:
        style = _colour_n_to_c("mol")

    # The ligand is the thing the eye should land on, in both cases, so it
    # gets green carbons against chain colours that avoid green.
    cmd.color("green", "mol and organic and elem C")
    cmd.color("grey70", "mol and inorganic")
    cmd.util.cnc("mol and organic")

    # ---- view ----
    # orient on the polymer, not the whole file: a ligand sitting out at the
    # edge of the asymmetric unit would otherwise swing the camera onto it.
    cmd.orient("mol and polymer")
    # complete=1 fits the drawn cartoon rather than the atom centres, so a
    # wide ribbon does not get clipped at the edge; the small buffer is what
    # keeps the molecule filling the frame instead of floating in white.
    cmd.zoom("mol", buffer=1.5, complete=1)

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ambient", 0.25)
    cmd.set("specular", 0.15)
    cmd.set("ray_shadows", 0)       # shadows read as surface detail that is not there
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_mode", 0)

    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    cmd.png(out_png, width=width, height=height, dpi=150, ray=1)

    if not os.path.exists(out_png):
        raise SystemExit("ERROR: PyMOL wrote no PNG for {}".format(structure))

    source = "prediction" if _is_prediction(stem) else "experimental"
    print("\t".join([
        out_png, stem, source, str(len(chains)), str(residues),
        ",".join(ligands) or "-", ",".join(metals) or "-", style,
    ]))


main()
