# Example provenance

- `data/1ubq.pdb`: https://files.rcsb.org/download/1UBQ.pdb
- `data/1ubq.cif`: https://files.rcsb.org/download/1UBQ.cif
- Entry: https://www.rcsb.org/structure/1UBQ
- Retrieved: 2026-11-13.
- These files contain a public experimental ubiquitin structure, not user research data.
- RCSB PDB data usage policies: https://www.rcsb.org/pages/policies
- `domains_1ubq.json` defines illustrative ranges 1–35 and 36–76 to exercise software. They are not curated biological domains.
- Generate the demonstration with `python protein_viz.py examples/data/1ubq.cif --domains @examples/domains_1ubq.json -o outputs/1ubq.png`.
- PNG/PSE/JSON outputs are generated locally in the ignored `outputs/` folder.
