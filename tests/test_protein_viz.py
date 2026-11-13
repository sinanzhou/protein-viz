import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from PIL import Image
import pytest

import protein_viz as pv

DATA = Path(__file__).resolve().parents[1] / "examples" / "data"


def test_cli_without_structure():
    result = subprocess.run([sys.executable, str(DATA.parents[1] / "protein_viz.py"),
                             "--list-schemes"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "earth" in result.stdout


def test_pdb_and_mmcif_agree():
    pdb, cif = [pv.parse_structure(DATA / f"1ubq.{ext}") for ext in ("pdb", "cif")]
    assert len(pdb["A"]) == 76
    assert [r["id"] for r in pdb["A"]] == [r["id"] for r in cif["A"]]
    np.testing.assert_allclose([r["xyz"] for r in pdb["A"]], [r["xyz"] for r in cif["A"]])


@pytest.mark.parametrize("domains", [[], {"bad": [1, 999]}, {"bad": [10, 1]},
    {"x": [1, 20], "y": [20, 30]}, {"bad": [True, 20]}, {"bad": [1.5, 3]}, {"": [1, 2]}])
def test_invalid_domains(domains):
    with pytest.raises(ValueError):
        pv.validate_domains(domains, pv.parse_structure(DATA / "1ubq.pdb")["A"])


@pytest.mark.parametrize("kwargs", [{"chain": "B"}, {"transparency": -1},
    {"transparency": float("nan")}, {"scheme": "oops"}, {"width": 0},
    {"dpi": 0}, {"font_size": 0}, {"rotation": [float("inf"), 0]},
    {"scheme": "highlight", "domains": {"x": [1, 20]}}])
def test_bad_parameters_fail_before_render(tmp_path, kwargs):
    with pytest.raises(ValueError):
        pv.render(DATA / "1ubq.pdb", tmp_path / "bad.png", **kwargs)
    assert not list(tmp_path.iterdir())


def test_first_model_and_insertion_codes(tmp_path):
    # Explicit synthetic fixture: different chains/models and an insertion code.
    atom = "ATOM  {serial:5d}  CA  ALA {chain}{resi:4d}{icode}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
    text = "MODEL        1\n"
    text += atom.format(serial=1, chain="B", resi=10, icode=" ", x=0, y=0, z=0)
    text += atom.format(serial=2, chain="B", resi=10, icode="A", x=3, y=0, z=0)
    text += "ENDMDL\nMODEL        2\n"
    text += atom.format(serial=3, chain="Z", resi=99, icode=" ", x=99, y=0, z=0)
    text += "ENDMDL\nEND\n"
    path = tmp_path / "models.pdb"
    path.write_text(text)
    result = pv.parse_structure(path)
    assert list(result) == ["B"]
    assert [r["id"] for r in result["B"]] == ["10", "10A"]


def test_multichain_requires_choice(tmp_path):
    source = (DATA / "1ubq.pdb").read_text()
    atoms = [line for line in source.splitlines() if line.startswith("ATOM")]
    path = tmp_path / "two.pdb"
    path.write_text("\n".join(atoms + [line[:21] + "B" + line[22:] for line in atoms]) + "\nEND\n")
    with pytest.raises(ValueError, match="Multiple chains"):
        pv.render(path, tmp_path / "out.png")


def test_actual_cif_render_and_session(tmp_path):
    import pymol2
    path = tmp_path / "含空格 demo.png"
    domains = {"region one": [1, 35], "region / two": [36, 76]}
    pv.render(DATA / "1ubq.cif", path, domains=domains, width=480, height=360, dpi=300, font_size=8)
    manifest = json.loads(path.with_suffix(".json").read_text(encoding="utf-8"))
    assert manifest["residue_count"] == 76
    assert manifest["input_sha256"] == hashlib.sha256((DATA / "1ubq.cif").read_bytes()).hexdigest()
    assert manifest["output_dimensions"] == [740, 360]
    assert manifest["missing_residue_numbers_in_domains"] == {name: [] for name in domains}
    with Image.open(path) as image:
        assert image.size == (740, 360)
        assert abs(image.info["dpi"][0] - 300) < .1
    with Image.open(path.with_name(path.stem + ".raw.png")) as image:
        pixels = np.asarray(image.convert("RGB"))
        assert np.all(pixels[[0, -1]] > 245) and np.all(pixels[:, [0, -1]] > 245)
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(str(path.with_suffix(".pse")))
        assert len(pymol.cmd.get_model("pv_protein and name CA").atom) == 76
        np.testing.assert_allclose(pymol.cmd.get_view(), manifest["view"], atol=1e-5)
    with pytest.raises(FileExistsError):
        pv.render(DATA / "1ubq.cif", path)


def test_actual_non_a_chain_and_highlight(tmp_path):
    source = (DATA / "1ubq.pdb").read_text()
    path = tmp_path / "chain_b.pdb"
    path.write_text("\n".join(line[:21] + "B" + line[22:] if line.startswith("ATOM") else line
                                for line in source.splitlines()))
    output = tmp_path / "b.png"
    pv.render(path, output, scheme="highlight", domains={"target": [1, 35]},
              highlight_domain="target", annotate=False, width=320, height=240)
    manifest = json.loads(output.with_suffix(".json").read_text(encoding="utf-8"))
    assert manifest["chain"] == "B" and manifest["domain_colors"]["target"] == "#CC3333"
    assert manifest["output_dimensions"] == [320, 240]


@pytest.mark.parametrize("size", [(480, 320), (320, 480)])
def test_projection_matches_real_raytraced_marker(tmp_path, size):
    import pymol2
    width, height = size
    xyz = [4., -3., 2.]
    output = tmp_path / "marker.png"
    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd
        cmd.pseudoatom("marker", pos=xyz, color="red", vdw=.5)
        cmd.show("spheres", "marker")
        cmd.bg_color("white")
        cmd.set("orthoscopic", 1)
        cmd.set("ray_opaque_background", 1)
        cmd.viewport(width, height)
        cmd.turn("z", -30)
        cmd.turn("x", 15)
        cmd.zoom("marker", buffer=10, complete=1)
        # Move the marker away from the camera origin to test more than centering.
        cmd.translate([2, 1, 0], "marker", camera=0)
        actual_xyz = cmd.get_model("marker").atom[0].coord
        predicted = pv.project_point(actual_xyz, cmd.get_view(), width, height, float(cmd.get("field_of_view")))
        cmd.png(str(output), width=width, height=height, ray=1)
    pixels = np.asarray(Image.open(output).convert("RGB"))
    yy, xx = np.where((pixels[:, :, 0] > 100) & (pixels[:, :, 1] < 80) & (pixels[:, :, 2] < 80))
    assert len(xx) > 5
    np.testing.assert_allclose(predicted, [xx.mean(), yy.mean()], atol=2)
