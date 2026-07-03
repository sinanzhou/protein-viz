#!/usr/bin/env python3
"""
ProteinViz: Publication-quality protein 3D structure visualization.
Domain-colored cartoon + semi-transparent surface overlay with annotations.

Author: Zikai Zhou
License: MIT
"""

import subprocess
import sys
import os
import argparse
import json
import numpy as np
from pathlib import Path

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
except ImportError:
    print("Error: matplotlib required. Install: pip install matplotlib")
    sys.exit(1)

try:
    from PIL import Image
except ImportError:
    print("Error: Pillow required. Install: pip install Pillow")
    sys.exit(1)


# === Built-in color schemes ===
SCHEMES = {
    "earth": {
        "name": "Earth Tones (Morandi)",
        "colors": {
            "AB": [190, 120, 130],
            "C":  [210, 185, 140],
            "D":  [130, 165, 130],
            "E":  [110, 120, 170],
        },
        "labels": {
            "AB": "#BE7882",
            "C":  "#D2B98C",
            "D":  "#82A582",
            "E":  "#6E78AA",
        }
    },
    "pastel": {
        "name": "Soft Pastels",
        "colors": {
            "AB": [180, 148, 210],
            "C":  [235, 195, 82],
            "D":  [100, 195, 195],
            "E":  [220, 120, 100],
        },
        "labels": {
            "AB": "#B494D2",
            "C":  "#EBC352",
            "D":  "#64C3C3",
            "E":  "#DC7864",
        }
    },
    "jewel": {
        "name": "Deep Jewel Tones",
        "colors": {
            "AB": [120, 80, 170],
            "C":  [210, 160, 50],
            "D":  [50, 150, 170],
            "E":  [180, 50, 60],
        },
        "labels": {
            "AB": "#7850AA",
            "C":  "#D2A032",
            "D":  "#3296AA",
            "E":  "#B4323C",
        }
    },
    "rainbow": {
        "name": "Rainbow (N→C)",
        "colors": None,  # handled specially
        "labels": {
            "AB": "#8B00FF",
            "C":  "#00BFFF",
            "D":  "#32CD32",
            "E":  "#FF4500",
        }
    },
    "highlight": {
        "name": "Gray + Target Highlight",
        "colors": None,  # gray everything, highlight E
        "labels": {
            "AB": "#999999",
            "C":  "#999999",
            "D":  "#999999",
            "E":  "#CC3333",
        }
    },
}

DEFAULT_DOMAINS = {
    "AB": (1, 190),
    "C":  (191, 316),
    "D":  (317, 365),
    "E":  (366, 603),
}


def find_pymol():
    """Find PyMOL executable."""
    for name in ["pymol", "pymol.exe"]:
        result = subprocess.run(["which", name], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
    raise FileNotFoundError("PyMOL not found. Install: conda install -c conda-forge pymol-open-source")


def parse_structure(filepath):
    """Parse PDB or mmCIF file, return chain and residue info."""
    filepath = str(filepath)
    with open(filepath) as f:
        first_line = f.readline()

    is_cif = filepath.endswith('.cif') or filepath.endswith('.mmcif')

    chains = {}
    with open(filepath) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if is_cif:
                parts = line.split()
                if len(parts) < 13:
                    continue
                atom_name = parts[3]
                chain = parts[6]
                try:
                    resi = int(parts[8])
                except (ValueError, IndexError):
                    continue
                x, y, z = float(parts[10]), float(parts[11]), float(parts[12])
            else:
                if len(line) < 54:
                    continue
                atom_name = line[12:16].strip()
                chain = line[21]
                try:
                    resi = int(line[22:26].strip())
                except ValueError:
                    continue
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])

            if atom_name != "CA":
                continue
            if chain not in chains:
                chains[chain] = []
            chains[chain].append((resi, x, y, z))

    return chains


def generate_pymol_script(filepath, output_png, domains, scheme_name, 
                          transparency=0.55, resolution=(3200, 2000)):
    """Generate PyMOL script for rendering."""
    scheme = SCHEMES.get(scheme_name, SCHEMES["earth"])
    w, h = resolution

    lines = [
        f'load {filepath}, prot',
        'hide all',
        'bg_color white',
    ]

    # Domain selections
    for name, (start, end) in domains.items():
        lines.append(f'select {name.lower()}, chain A and resi {start}-{end}')

    lines += [
        'show cartoon, prot',
        'show surface, prot',
        'color white, chain A',
    ]

    if scheme_name == "rainbow":
        lines.append('spectrum rainbow, chain A')
    elif scheme_name == "highlight":
        lines.append('color gray70, chain A')
        e_start, e_end = domains["E"]
        lines.append(f'color red, chain A and resi {e_start}-{e_end}')
    else:
        for name, rgb in scheme["colors"].items():
            sel = name.lower()
            lines.append(f'set_color c_{sel}, [{rgb[0]}, {rgb[1]}, {rgb[2]}]')
            lines.append(f'color c_{sel}, {sel}')

    lines += [
        f'set transparency, {transparency}',
        'set surface_mode, 1',
        'set surface_solvent, off',
        'set cartoon_fancy_helices, on',
        'set cartoon_smooth_loops, on',
        'set cartoon_width, 2.0',
        'set specular, 0.8',
        'set shininess, 60',
        'set ambient, 0.45',
        'set direct, 0.55',
        'set depth_cue, 0',
        'set antialias, 2',
        'center prot',
        'zoom prot',
        'turn z, -30',
        'turn x, 15',
        f'ray {w // 2}, {h // 2}',
        f'png {output_png}, {w}, {h}, dpi=150',
    ]

    return "
".join(lines)


def project_domain_centers(filepath, domains):
    """Project 3D domain centers to 2D pixel coordinates."""
    chains = parse_structure(filepath)
    chain_a = chains.get("A", [])
    if not chain_a:
        # Try first available chain
        chain_a = list(chains.values())[0] if chains else []

    centers_3d = {}
    for name, (start, end) in domains.items():
        pts = [(x, y, z) for resi, x, y, z in chain_a if start <= resi <= end]
        if pts:
            centers_3d[name] = np.mean(pts, axis=0)

    # N-term and C-term
    if chain_a:
        resis = [r for r, *_ in chain_a]
        min_resi = min(resis)
        max_resi = max(resis)
        for resi, x, y, z in chain_a:
            if resi == min_resi:
                centers_3d["NT"] = np.array([x, y, z])
            if resi == max_resi:
                centers_3d["CT"] = np.array([x, y, z])

    if not centers_3d:
        return {}

    # Rotation: turn z -30, turn x 15
    def rot_z(a):
        r = np.radians(a)
        return np.array([[np.cos(r), -np.sin(r), 0],
                         [np.sin(r), np.cos(r), 0], [0, 0, 1]])

    def rot_x(a):
        r = np.radians(a)
        return np.array([[1, 0, 0], [0, np.cos(r), -np.sin(r)],
                         [0, np.sin(r), np.cos(r)]])

    R = rot_x(15) @ rot_z(-30)

    keys = list(centers_3d.keys())
    all_proj = np.array([R @ centers_3d[k] for k in keys])

    canvas_w, canvas_h = 3200, 2000
    x_range = all_proj[:, 0].max() - all_proj[:, 0].min()
    y_range = all_proj[:, 1].max() - all_proj[:, 1].min()
    scale = min((canvas_w * 0.65) / max(x_range, 1),
                (canvas_h * 0.65) / max(y_range, 1))
    x_c, y_c = all_proj[:, 0].mean(), all_proj[:, 1].mean()

    result = {}
    for i, k in enumerate(keys):
        px = int((all_proj[i, 0] - x_c) * scale + canvas_w / 2)
        py = int(-(all_proj[i, 1] - y_c) * scale + canvas_h / 2)
        result[k] = (px, py)

    return result


def add_annotations(image_path, output_path, domains, scheme_name, 
                    font_size=28, line_width=1.8):
    """Add domain labels with leader lines using Matplotlib."""
    scheme = SCHEMES.get(scheme_name, SCHEMES["earth"])
    labels_colors = scheme["labels"]

    bg = np.array(Image.open(image_path).convert("RGB"))
    h, w = bg.shape[:2]

    targets = project_domain_centers(image_path.replace(".png", ".png"), domains)
    # Fallback: re-project from the rendered image's metadata
    # For simplicity, use PyMOL-labeled image approach

    fig, ax = plt.subplots(figsize=(32, 20), dpi=150)
    ax.imshow(bg)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis('off')

    # Label positions: left side for most, right side for N/C term
    label_configs = []
    y_positions = {"AB": 0.2, "C": 0.35, "D": 0.5, "E": 0.65}

    for name in ["AB", "C", "D", "E"]:
        lx = w * 0.12
        ly = h * y_positions[name]
        color = labels_colors.get(name, "#333333")
        display = name.replace("AB", "A/B")
        if name == "E":
            display += "
(LBD)"
        label_configs.append((display, lx, ly, targets.get(name), color))

    # N-term and C-term on right side
    if "NT" in targets:
        label_configs.append(("N-term", w * 0.85, h * 0.15, targets["NT"], "#555555"))
    if "CT" in targets:
        label_configs.append(("C-term", w * 0.85, h * 0.8, targets["CT"], "#555555"))

    for text, lx, ly, target, color in label_configs:
        if target is None:
            continue
        tx, ty = target
        ax.plot(tx, ty, 'o', color=color, markersize=6, zorder=4)
        ax.plot([tx, lx], [ty, ly], color=color, linewidth=line_width,
                linestyle='--', alpha=0.7, zorder=3)
        ax.text(lx, ly, text, fontsize=font_size, fontweight='bold', color=color,
                ha='center', va='center',
                path_effects=[pe.withStroke(linewidth=5, foreground='white')], zorder=5)

    fig.savefig(output_path, dpi=150, bbox_inches='tight', pad_inches=0.15,
                facecolor='white')
    plt.close()


def render(input_file, output_file=None, scheme="earth", domains=None,
           transparency=0.55, annotate=True, font_size=28):
    """Main render pipeline."""
    pymol = find_pymol()
    input_path = Path(input_file).resolve()

    if not input_path.exists():
        print(f"Error: {input_path} not found")
        sys.exit(1)

    if domains is None:
        domains = DEFAULT_DOMAINS

    if output_file is None:
        output_file = input_path.parent / f"{input_path.stem}_{scheme}_annotated.png"
    output_path = Path(output_file)

    tmp_dir = Path("/tmp/protein_viz_tmp")
    tmp_dir.mkdir(exist_ok=True)
    raw_png = tmp_dir / "raw_render.png"
    pml_script = tmp_dir / "render.pml"

    # Step 1: PyMOL render
    print(f"[1/3] Rendering with PyMOL (scheme: {scheme})...")
    script = generate_pymol_script(str(input_path), str(raw_png), domains, scheme, transparency)
    pml_script.write_text(script)

    result = subprocess.run([pymol, "-cq", str(pml_script)],
                            capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        print(f"PyMOL error: {result.stderr[-500:]}")
        sys.exit(1)

    # Step 2: Add annotations
    if annotate:
        print("[2/3] Adding annotations...")
        add_annotations(str(raw_png), str(output_path), domains, scheme, font_size)
    else:
        import shutil
        shutil.copy(str(raw_png), str(output_path))

    print(f"[3/3] Done! Output: {output_path}")
    return str(output_path)


def main():
    parser = argparse.ArgumentParser(
        description="ProteinViz: Publication-quality protein structure figures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python protein_viz.py structure.pdb
  python protein_viz.py structure.cif --scheme earth --transparency 0.55
  python protein_viz.py model.cif --scheme pastel --no-annotate
  python protein_viz.py protein.pdb --domains '{"AB": [1,190], "C": [191,316], "D": [317,365], "E": [366,603]}'

Available schemes: earth, pastel, jewel, rainbow, highlight
        """)

    parser.add_argument("input", help="PDB or mmCIF structure file")
    parser.add_argument("-o", "--output", help="Output PNG path")
    parser.add_argument("-s", "--scheme", default="earth",
                        choices=list(SCHEMES.keys()),
                        help="Color scheme (default: earth)")
    parser.add_argument("-t", "--transparency", type=float, default=0.55,
                        help="Surface transparency, 0=opaque 1=invisible (default: 0.55)")
    parser.add_argument("--no-annotate", action="store_true",
                        help="Skip label annotations")
    parser.add_argument("--font-size", type=int, default=28,
                        help="Label font size (default: 28)")
    parser.add_argument("--domains", type=str, default=None,
                        help='Domain ranges as JSON, e.g. '{"AB":[1,190],"C":[191,316]}'')
    parser.add_argument("--list-schemes", action="store_true",
                        help="List available color schemes and exit")

    args = parser.parse_args()

    if args.list_schemes:
        print("Available color schemes:")
        for name, scheme in SCHEMES.items():
            print(f"  {name:12s} — {scheme['name']}")
        return

    domains = DEFAULT_DOMAINS
    if args.domains:
        domains = json.loads(args.domains)

    render(args.input, args.output, args.scheme, domains,
           args.transparency, not args.no_annotate, args.font_size)


if __name__ == "__main__":
    main()
