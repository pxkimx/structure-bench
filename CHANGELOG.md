# Changelog

## 0.1.0 — first release (2026-09-23)
- Look up a protein by gene symbol + organism (human, mouse, rat, zebrafish, fly, worm, yeast, E. coli or a taxon ID)
  or by UniProt accession; deep links `#gene=LMNA&species=human` and `#uniprot=P02545` (used by RNAseq Bench).
- UniProt panel (function, location, disease, grouped features), AlphaFold confidence track with domains and
  variants, PAE heatmap with domain-pair summary, 3D viewer (3Dmol.js 2.5.5, vendored), PDB coverage with one-click load.
- Variant mapping (R482W / p.Arg482Trp / G608G …): reference check with signal-peptide / isoform suggestions,
  pLDDT or B-factor, relative SASA (Shrake–Rupley / Tien 2013), interface and ligand contacts, neighbours,
  secondary structure (DSSP or a labelled phi/psi approximation), UniProt annotation, AlphaMissense.
- Structure comparison (CE + per-residue CA deviation with pLDDT overlay), uploads (.pdb/.cif/.bcif, offline).
- Downloads: B-factor-replaced structures, PyMOL/ChimeraX scripts, PDF report. Claude assistant with 10 tools.
- Offline: every download is cached; the LMNA example is bundled; `python -m server.selftest`.
