# Changes

## 0.2.0

- Fix the syntax errors that prevented the July 2026 script from starting.
- Parse PDB/mmCIF with Biopython, using first-model author chain and residue numbering.
- Select a single protein chain explicitly when multiple chains are present.
- Make the legacy AgEcR1 ranges an explicit preset instead of applying them to every input.
- Support custom range names and validate interval endpoints, overlaps and numeric settings.
- Anchor range labels on real CA coordinates using the actual PyMOL orthoscopic camera.
- Place labels in a separate gutter and fit the model after rotation with an additional surface buffer.
- Use the PyMOL Python API with an independent instance per call instead of fixed temporary files and a Unix-only executable lookup.
- Export a raw PNG, annotated PNG, editable PSE and input/parameter provenance JSON.
- Add a public ubiquitin example with source attribution, a validated environment snapshot and 24 regression checks, including real ray-traced renders and camera calibration.

Compatibility: `render` retains its first seven positional arguments, but no longer assigns AgEcR1 domains by default. Raw image dimensions are configurable; annotation adds a gutter. Font size now uses Matplotlib points with a default of 14. The old approximate projection and text-script-generation helpers have been removed. The PyMOL Python module is required for rendering.
