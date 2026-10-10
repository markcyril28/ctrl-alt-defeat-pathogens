#!/usr/bin/env python3
"""Print the version of every Python library this pipeline depends on.

    library_versions.py

One line per library, "<name> <version>", and a line saying so for any that
is not installed rather than stopping.

These belong in the version record next to the command-line tools, because
they do the structural work: PyMOL cleans the receptor and measures the
models, RDKit builds the ligand conformers, Meeko assigns the atom types and
charges that Vina scores. A score is not reproducible without them, and
Meeko's mk_prepare_* scripts take no --version flag of their own, so this is
the only place its version is recorded.
"""
import importlib
import sys

LIBRARIES = ("pymol", "rdkit", "meeko", "Bio", "pandas")


def main():
    for name in LIBRARIES:
        try:
            library = importlib.import_module(name)
        except Exception as error:
            print(f"{name} (not available: {error.__class__.__name__})")
            continue

        if name == "pymol":
            # PyMOL keeps its version behind cmd.get_version(), not __version__.
            print(f"pymol {library.cmd.get_version()[0]}")
        else:
            print(f"{name} {getattr(library, '__version__', '(no __version__)')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
