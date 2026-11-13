#!/usr/bin/env python3
"""Reproducible protein structure figures. MIT License, Zikai Zhou."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

VERSION = "0.2.0"
DEFAULT_DOMAINS = {"AB": (1, 190), "C": (191, 316), "D": (317, 365), "E": (366, 603)}
SCHEMES = {
    "earth": ["#BE7882", "#D2B98C", "#82A582", "#6E78AA"],
    "pastel": ["#B494D2", "#EBC352", "#64C3C3", "#DC7864"],
    "jewel": ["#7850AA", "#D2A032", "#3296AA", "#B4323C"],
    "rainbow": ["#333333"], "highlight": ["#999999"],
}


def parse_structure(filepath):
    """First model, author chain/residue IDs, one selected CA per protein residue."""
    from Bio.PDB import MMCIFParser, PDBParser
    path = Path(filepath)
    if not path.is_file():
        raise FileNotFoundError(f"Structure not found: {path}")
    if path.suffix.lower() in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True)
    elif path.suffix.lower() in {".pdb", ".ent"}:
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError("Input must be .pdb, .ent, .cif or .mmcif")
    structure = parser.get_structure("input", str(path))
    model = next(structure.get_models(), None)
    chains = {}
    if model is not None:
        for chain in model:
            residues = []
            for residue in chain:
                if residue.id[0] == " " and "CA" in residue:
                    xyz = [float(v) for v in residue["CA"].coord]
                    if not all(math.isfinite(v) for v in xyz):
                        raise ValueError("Non-finite CA coordinates")
                    residues.append({"number": residue.id[1], "id": str(residue.id[1]) + residue.id[2].strip(), "xyz": xyz})
            if residues:
                chains[chain.id.strip()] = residues
    if not chains:
        raise ValueError("No protein CA residues in first model")
    return chains


def validate_domains(domains, residues):
    if not isinstance(domains, dict):
        raise ValueError("Domains must map names to [start, end]")
    validated, previous = {}, []
    present = {r["number"] for r in residues}
    for name, bounds in domains.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Domain names must be nonempty strings")
        if (not isinstance(bounds, (list, tuple)) or len(bounds) != 2
                or any(type(v) is not int for v in bounds)):
            raise ValueError(f"Domain {name!r} needs two integer residue numbers")
        start, end = bounds
        if start > end:
            raise ValueError(f"Domain {name!r}: start exceeds end")
        if start not in present or end not in present:
            raise ValueError(f"Domain {name!r}: endpoint absent from selected chain")
        if any(start <= b and end >= a for a, b in previous):
            raise ValueError(f"Domain {name!r} overlaps another domain")
        previous.append((start, end))
        validated[name] = (start, end)
    return validated


def domain_colors(domains, scheme, highlight_domain=None):
    palette = SCHEMES[scheme]
    return {name: ("#CC3333" if name == highlight_domain else "#999999")
            if scheme == "highlight" else palette[i % len(palette)]
            for i, name in enumerate(domains)}


def project_point(xyz, view, width, height, field_of_view):
    """Use the actual orthoscopic PyMOL camera."""
    import numpy as np
    # get_view stores this matrix in column-major order.
    camera = np.asarray(view[:9]).reshape(3, 3).T @ (np.asarray(xyz) - np.asarray(view[12:15])) + np.asarray(view[9:12])
    # PyMOL's stored field of view is vertical, including portrait viewports.
    scale = height / (2 * abs(view[11]) * math.tan(math.radians(field_of_view / 2)))
    return [float(width / 2 + camera[0] * scale), float(height / 2 - camera[1] * scale)]


def add_annotations(raw_image, output_path, domains, colors, targets, dpi, font_size):
    """Labels occupy a separate gutter; keep exact raster dimensions."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image
    with Image.open(raw_image) as image:
        bg = image.convert("RGB")
    w, h = bg.size
    gutter = max(int(w * .25), 260)
    total = w + gutter
    fig = plt.figure(figsize=(total / dpi, h / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(bg, extent=(gutter, total, h, 0))
    ax.set(xlim=(0, total), ylim=(h, 0))
    ax.axis("off")
    names = sorted(domains, key=lambda name: targets[name][1])
    for i, name in enumerate(names):
        y = h * (i + 1) / (len(names) + 1)
        tx, ty = targets[name]
        start, end = domains[name]
        ax.annotate(f"{name}\n{start}\u2013{end}", xy=(gutter + tx, ty),
                    xytext=(gutter * .43, y), ha="center", va="center",
                    fontsize=font_size, color=colors[name],
                    arrowprops={"arrowstyle": "-", "color": colors[name], "lw": 1.2})
    try:
        fig.savefig(output_path, dpi=dpi, facecolor="white")
    finally:
        plt.close(fig)


def render(input_file, output_file=None, scheme="earth", domains=None,
           transparency=.55, annotate=True, font_size=14, chain=None,
           width=1600, height=1200, dpi=300, rotation=(-30., 15.),
           highlight_domain=None, overwrite=False):
    """Save PNG, raw PNG, PSE, JSON. None domains means a neutral whole chain.

    The API raises exceptions; it never exits the caller. AgEcR1 ranges are opt-in.
    """
    from PIL import Image
    if scheme not in SCHEMES:
        raise ValueError(f"Unknown scheme: {scheme}")
    if not math.isfinite(transparency) or not 0 <= transparency <= 1:
        raise ValueError("Transparency must be between 0 and 1")
    if any(type(v) is not int or not 64 <= v <= 10000 for v in (width, height)):
        raise ValueError("Width and height must be integers between 64 and 10000")
    if type(dpi) is not int or not 1 <= dpi <= 2400:
        raise ValueError("DPI must be an integer between 1 and 2400")
    if not math.isfinite(font_size) or not 1 <= font_size <= 100:
        raise ValueError("Font size must be between 1 and 100 points")
    if len(rotation) != 2 or not all(math.isfinite(v) for v in rotation):
        raise ValueError("Rotation requires finite Z and X angles")
    input_path = Path(input_file).resolve()
    chains = parse_structure(input_path)
    if chain is None:
        if len(chains) != 1:
            raise ValueError(f"Multiple chains {list(chains)}; select one with --chain")
        chain = next(iter(chains))
    if chain not in chains:
        raise ValueError(f"Chain {chain!r} absent; available chains: {list(chains)}")
    residues = chains[chain]
    domains = validate_domains({} if domains is None else domains, residues)
    if scheme == "highlight":
        if highlight_domain is None and "E" in domains:
            highlight_domain = "E"
        if highlight_domain not in domains:
            raise ValueError("Highlight scheme requires --highlight-domain naming a supplied domain")
    colors = domain_colors(domains, scheme, highlight_domain)
    output_path = Path(output_file or input_path.with_name(f"{input_path.stem}_{scheme}.png")).resolve()
    if output_path.suffix.lower() != ".png":
        raise ValueError("Output must have a .png extension")
    paths = {"image": output_path, "raw": output_path.with_name(output_path.stem + ".raw.png"),
             "session": output_path.with_suffix(".pse"), "manifest": output_path.with_suffix(".json")}
    if input_path in paths.values():
        raise ValueError("Output would overwrite input")
    if not overwrite and any(p.exists() for p in paths.values()):
        raise FileExistsError("Output or sidecar exists; use another name or --overwrite")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import pymol2
    except ImportError as exc:
        raise RuntimeError("PyMOL module missing; install requirements.txt in a supported environment") from exc
    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd
        cmd.set("cif_use_auth", 1)
        cmd.load(str(input_path), "source")
        atoms = cmd.get_model("source", state=1).atom
        indices = [a.index for a in atoms if a.chain.strip() == chain and not a.hetatm]
        if not indices:
            raise ValueError("PyMOL could not select author chain in first model")
        cmd.select_list("chosen", "source", indices, mode="index")
        cmd.create("pv_protein", "chosen", 1, 1)
        cmd.delete("source")
        cmd.delete("chosen")
        cmd.remove("pv_protein and not (alt ''+A)")
        cmd.hide("everything", "all")
        cmd.bg_color("white")
        cmd.color("gray80", "pv_protein")
        cmd.show("cartoon", "pv_protein")
        cmd.show("surface", "pv_protein")
        cmd.dss("pv_protein")
        atoms = cmd.get_model("pv_protein", state=1).atom
        ca = {a.resi: list(a.coord) for a in atoms if a.name == "CA"}
        if set(ca) != {r["id"] for r in residues}:
            raise ValueError("Parser/PyMOL residue identities differ; inspect alternate conformers or numbering")
        anchors = {}
        for i, (name, (start, end)) in enumerate(domains.items()):
            selection = f"domain_{i}"
            residue_ids = {r["id"] for r in residues if start <= r["number"] <= end}
            cmd.select_list(selection, "pv_protein", [a.index for a in atoms if a.resi in residue_ids], mode="index")
            if scheme != "rainbow":
                color = colors[name]
                cmd.set_color(f"pv_color_{i}", [int(color[j:j+2], 16) / 255 for j in (1, 3, 5)])
                cmd.color(f"pv_color_{i}", selection)
            points = [ca[r["id"]] for r in residues if start <= r["number"] <= end]
            mean = [sum(p[j] for p in points) / len(points) for j in range(3)]
            anchors[name] = min(points, key=lambda p: sum((p[j] - mean[j]) ** 2 for j in range(3)))
        if scheme == "rainbow":
            cmd.spectrum("resi", "rainbow", "pv_protein")
        for setting, value in {"transparency": transparency, "surface_solvent": 0,
                               "cartoon_fancy_helices": 1, "depth_cue": 0, "antialias": 2,
                               "orthoscopic": 1, "ray_opaque_background": 1}.items():
            cmd.set(setting, value)
        cmd.viewport(width, height)
        cmd.orient("pv_protein")
        cmd.turn("z", rotation[0])
        cmd.turn("x", rotation[1])
        cmd.zoom("pv_protein", buffer=5, complete=1)
        view, fov = list(cmd.get_view()), float(cmd.get("field_of_view"))
        targets = {name: project_point(xyz, view, width, height, fov) for name, xyz in anchors.items()}
        cmd.png(str(paths["raw"]), width=width, height=height, dpi=dpi, ray=1, quiet=1)
        cmd.save(str(paths["session"]))
        pymol_version = cmd.get_version()[0]
    if annotate and domains:
        add_annotations(paths["raw"], output_path, domains, colors, targets, dpi, font_size)
    else:
        with Image.open(paths["raw"]) as image:
            image.save(output_path, dpi=(dpi, dpi))
    with Image.open(output_path) as image:
        dimensions = list(image.size)
    present = {r["number"] for r in residues}
    manifest = {
        "protein_viz_version": VERSION, "pymol_version": pymol_version,
        "input_name": input_path.name, "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "model": "first", "chain": chain, "residue_count": len(residues),
        "residue_ids": [r["id"] for r in residues], "domains": domains,
        "missing_residue_numbers_in_domains": {name: [n for n in range(a, b+1) if n not in present] for name, (a, b) in domains.items()},
        "scheme": scheme, "domain_colors": colors, "highlight_domain": highlight_domain,
        "surface_transparency": transparency, "rotation_z_x": list(rotation),
        "orthoscopic": True, "view": view, "field_of_view": fov,
        "anchor_xyz": anchors, "anchor_pixels_in_raw_image": targets,
        "raw_dimensions": [width, height], "output_dimensions": dimensions,
        "dpi": dpi, "annotate": bool(annotate and domains), "font_size": font_size,
        "outputs": {key: path.name for key, path in paths.items()},
        "note": "Visualization only; no folding, domain inference or confidence validation.",
    }
    paths["manifest"].write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return str(output_path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", help="PDB/mmCIF input")
    parser.add_argument("-o", "--output", help="Output PNG; sidecars use the same stem")
    parser.add_argument("-s", "--scheme", choices=SCHEMES, default="earth")
    parser.add_argument("--domains", help="JSON object, or @path/to/domains.json")
    parser.add_argument("--preset", choices=["agecr1"], help="Explicit AgEcR1 ranges, not inferred")
    parser.add_argument("--chain", help="Author chain ID (required for multiple chains)")
    parser.add_argument("-t", "--transparency", type=float, default=.55)
    parser.add_argument("--no-annotate", action="store_true")
    parser.add_argument("--font-size", type=float, default=14)
    parser.add_argument("--width", type=int, default=1600, help="Raw width; annotation adds a gutter")
    parser.add_argument("--height", type=int, default=1200)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--rotation", nargs=2, type=float, default=(-30., 15.), metavar=("Z", "X"))
    parser.add_argument("--highlight-domain")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--inspect", action="store_true", help="Report first-model CA residues without PyMOL")
    parser.add_argument("--list-schemes", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)
    if args.list_schemes:
        print("\n".join(SCHEMES))
        return 0
    if not args.input:
        parser.error("input required unless --list-schemes or --version is used")
    try:
        if args.inspect:
            print(json.dumps(parse_structure(args.input), indent=2))
            return 0
        if args.domains and args.preset:
            raise ValueError("Choose --domains or --preset, not both")
        domains = DEFAULT_DOMAINS if args.preset else None
        if args.domains:
            source = args.domains
            if source.startswith("@"):
                source = Path(source[1:]).read_text(encoding="utf-8-sig")
            domains = json.loads(source)
            if not isinstance(domains, dict):
                raise ValueError("Domains must be a JSON object")
        print(render(args.input, args.output, args.scheme, domains, args.transparency,
                     not args.no_annotate, args.font_size, args.chain, args.width,
                     args.height, args.dpi, args.rotation, args.highlight_domain, args.overwrite))
        return 0
    except (ValueError, OSError, RuntimeError, ImportError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
