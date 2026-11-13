# ProteinViz

Reproducible PDB/mmCIF visualization using PyMOL: cartoon + surface, optional residue-range colors and labels, an editable PyMOL session, and a JSON provenance record.

[中文说明](README.zh-CN.md)

## Install

Python 3.11 is recommended (the validated Windows environment uses Python 3.11).

```sh
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
python -m pip install -r requirements.txt
```

PyMOL wheels are platform-dependent. If pip cannot resolve PyMOL, use an environment from conda-forge with `pymol-open-source`, `biopython`, `matplotlib`, `pillow` and `numpy`. See [PyMOL upstream](https://github.com/schrodinger/pymol-open-source). `requirements-lock.txt` records the Windows validation environment; it is not a universal cross-platform lock.

## Usage

```sh
python protein_viz.py --list-schemes
python protein_viz.py examples/data/1ubq.cif --inspect
python protein_viz.py examples/data/1ubq.cif --domains @examples/domains_1ubq.json -o outputs/1ubq.png
python protein_viz.py examples/data/1ubq.pdb --scheme rainbow --no-annotate -o outputs/rainbow.png
python protein_viz.py examples/data/1ubq.pdb --domains @examples/domains_1ubq.json --scheme highlight --highlight-domain N-region -o outputs/highlight.png
```

The example splits ubiquitin into **illustrative residue ranges**, not experimentally established domains. See [example provenance](examples/README.md).

By default a structure is rendered in neutral gray with no inferred domains. Supply `--domains` for colored ranges. The legacy AgEcR1 intervals (AB 1–190, C 191–316, D 317–365, E 366–603) are available only through `--preset agecr1`; they must be checked against the actual input protein. ProteinViz visualizes an existing model; it neither predicts folds nor validates biological function, domain assignments or model confidence.

- First coordinate model only; author chain/residue numbering for mmCIF.
- Single protein chain is chosen automatically; multiple chains require `--chain`.
- Missing endpoints, overlapping intervals, invalid settings and mismatched residue identities cause an error. Interior numbering gaps are recorded in JSON; no missing residues are reconstructed.
- Alternate conformers: blank/A atoms are retained; inputs incompatible with this policy are rejected when CA identities differ. Review unusual alternate conformers before use.
- `--width` and `--height` specify the raw 3D render in pixels. Annotated output adds a left gutter of `max(width*0.25, 260)` pixels. `--dpi` defaults to 300.
- Every render produces `.png`, `.raw.png`, `.pse`, and `.json`. Existing outputs require `--overwrite`.
- Camera rotations are relative to PyMOL's initial orientation, not a promise to reproduce an older figure's camera. Final camera parameters are saved.

Labels are placed outside the structure and connected to real CA atoms near each range centroid. Anchors use the saved PyMOL orthoscopic camera. They are visual references, not a claim that the indicated atom lies on the visible surface. Crowded labels should still be reviewed manually.

## Python API

```python
from protein_viz import render
render("examples/data/1ubq.cif", "outputs/api.png",
       domains={"N-region": [1, 35], "C-region": [36, 76]},
       chain="A", width=1600, height=1200, dpi=300)
```

The API raises exceptions instead of exiting the caller. Palette names remain `earth`, `pastel`, `jewel`, `rainbow`, and `highlight`; custom range names are supported. Highlighting requires a supplied range name (legacy `E` is accepted when present).

## Validation

```sh
python -m pip install pytest
python -m pytest -q
```

Tests include real PyMOL renders from both PDB/mmCIF, non-A chains, editable-session reopening, DPI/dimensions, first-model and insertion-code parsing, output protection, invalid parameters and projection calibration against ray-traced markers in portrait/landscape viewports.

MIT License — Copyright (c) 2026 Zikai Zhou.
