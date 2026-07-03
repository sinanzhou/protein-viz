# 🧬 ProteinViz

Publication-quality protein 3D structure visualization with domain annotations.

Renders PDB/mmCIF files into beautiful cartoon + surface overlay figures with labeled domain annotations — ready for journal submission.

![Earth Tones](examples/preview_earth.png)

## ✨ Features

- 🎨 **5 color schemes**: Earth Tones, Soft Pastels, Jewel Tones, Rainbow, Highlight
- 🫧 **Cartoon + Surface overlay** — see both secondary structure and molecular surface
- 🏷️ **Automatic domain annotations** — labeled with leader lines
- 📐 **Publication-ready** — 300 DPI, clean white background
- 🔧 **Supports PDB & mmCIF** — works with AlphaFold2, AlphaFold3, experimental structures

## 🚀 Quick Start

### Install Dependencies

```bash
pip install matplotlib Pillow numpy

# PyMOL (required for 3D rendering)
conda install -c conda-forge pymol-open-source
# or: sudo apt install pymol
```

### Render a Structure

```bash
# Basic usage (default: earth tones scheme)
python protein_viz.py structure.pdb

# Choose a color scheme
python protein_viz.py structure.cif --scheme pastel
python protein_viz.py structure.cif --scheme jewel
python protein_viz.py structure.cif --scheme rainbow

# Adjust surface transparency (0=opaque, 1=invisible)
python protein_viz.py structure.pdb --transparency 0.55

# Custom domain ranges
python protein_viz.py structure.pdb --domains '{"AB":[1,190],"C":[191,316],"D":[317,365],"E":[366,603]}'

# Just the 3D render, no labels
python protein_viz.py structure.pdb --no-annotate
```

### List Available Schemes

```bash
python protein_viz.py --list-schemes
```

## 🎨 Color Schemes

| Scheme | Description |
|--------|-------------|
| `earth` | Morandi earth tones — dusty rose, warm sand, sage, slate blue |
| `pastel` | Soft pastels — lavender, gold, teal, coral |
| `jewel` | Deep jewel tones — violet, amber, teal, crimson |
| `rainbow` | Classic N→C rainbow gradient |
| `highlight` | Gray structure with red-highlighted target domain |

## 📖 For AI Agents

This tool is designed to be called by AI agents (like Hermes Agent) for automated protein figure generation:

```python
# Python API
from protein_viz import render

render(
    input_file="structure.pdb",
    output_file="figure.png",
    scheme="earth",
    transparency=0.55,
    annotate=True
)
```

### Hermes Agent Skill

To install as a Hermes Agent skill, add to `~/.hermes/skills/protein-viz/SKILL.md`:

```yaml
---
name: protein-viz
description: "Render publication-quality protein 3D structure figures with domain annotations."
version: 1.0.0
author: Zikai Zhou
license: MIT
---
```

Then place `protein_viz.py` in the skill directory and invoke via:
```
python ~/.hermes/skills/protein-viz/protein_viz.py input.pdb --scheme earth
```

## 🏗️ Requirements

- Python 3.8+
- PyMOL (open-source)
- matplotlib
- Pillow
- numpy

## 📋 Default Domain Ranges

The default domains are for **Aphis gossypii Ecdysone Receptor-like Protein (AgEcR1)**:

| Domain | Residues | Function |
|--------|----------|----------|
| A/B | 1–190 | Transactivation domain |
| C | 191–316 | DNA-binding domain (DBD) |
| D | 317–365 | Hinge region |
| E | 366–603 | Ligand-binding domain (LBD), pesticide target |

Use `--domains` to customize for your protein.

## 📜 License

MIT License — Copyright (c) 2026 Zikai Zhou

## 🙏 Acknowledgments

Built with [PyMOL](https://pymol.org/) and [Matplotlib](https://matplotlib.org/). Inspired by publication figures in Nature, Science, and Cell structural biology papers.
