# Changelog

## 0.1.3 — 2026-09-23
- **The PDF report button could fail outright.** The page is read from disk on every request, but the server's code is
  loaded once when it starts — so updating or rebuilding the app while it is running leaves a new page talking to an
  old server. The new page asks for the report at an address the old server does not have, and the button only said
  it could not be written. It now notices that case, writes the report through the older route instead (everything but
  the 3D pictures, which the new route carries), and says plainly: quit and reopen Structure Bench, then make it again.
  Any other failure now ends with the same hint.
- Both report routes are covered by a test, so the one the fallback needs cannot be removed by accident.

## 0.1.2 — 2026-09-23
**The PDF report now holds everything the page shows, and explains each part.** Before, it had the tiles,
the findings, the pLDDT track, a variant table and the comparison text. It left out the UniProt entry, the PAE,
the 3D model and the experimental structures. Every figure and table now carries the same two texts as on the
page: *What this shows* (how to read it) and *In your data* (what it says for this protein, computed from the
numbers).

The report now has:
- **Findings:** the page's findings, first.
- **Protein:** the UniProt entry (function, location, entry version, isoform models), disease involvement with
  OMIM numbers, and every annotated feature (domains and regions, sites, PTMs, natural variants).
- **Predicted structure:**
  - the AlphaFold model's details
  - a pLDDT band chart
  - the confidence/domain/variant track
  - the **PAE heatmap** with domain boundaries
  - **the 3D model**
  - a table of every structure and chain in the analysis.
- **Experimental structures:** PDB coverage along the sequence, and all entries in a table.
- **Variants:**
  - a lollipop chart (AlphaMissense score or pLDDT, coloured by burial)
  - **the variants on the 3D structure**
  - the summary table
  - every variant in detail (change, contacts, UniProt annotation, readings).
- **Comparisons:** the RMSDs, the per-residue deviation plot with pLDDT, **the superposition in 3D**, and the
  deviating stretches.

**3D pictures.** The page's viewer takes them when you press *PDF report*. It draws the model as a cartoon
coloured by pLDDT, laid along its longest axis and seen from two sides, plus the mapped variants and each
superposition. The server trims the white margins so each picture fills its space. When a report is made
without the page (the assistant, or a direct download), the server draws the CA trace from two directions
instead. Depth is shown by line width, a superposition is cropped to the compared region, and the report says
the picture is a stand-in.

The existing titin regression test caught a problem during this work. A reading longer than a page, placed in a
one-cell table, stopped the whole PDF. The reading box is now a shaded paragraph, which can break across pages.

## 0.1.1 — stress-test fixes (2026-09-23)
Found by running ~60 real look-ups three times, 1,000+ random variants, ~30 PDB entries checked against SIFTS,
comparisons, junk and huge uploads, concurrent and offline use, and the whole page at phone width in dark mode.
- **Gene symbols that are also someone else's synonym.** UniProt's exact gene search also matches synonyms, so `ALB`
  offered Fas-binding factor 1 first (ALB is one of its synonyms). The entry whose primary gene name is the query now
  wins, with a note naming the other; when several entries share the primary name (CDKN2A), synonym-only matches are
  listed last.
- **Fusion partners mapped onto the protein.** In β2-adrenergic receptor structures with T4 lysozyme (2RH1, 3SN6), up
  to 22 lysozyme residues were given receptor numbers, so a variant could land on the wrong protein. And a residue
  next to missing density could be placed three residues off (EGFR 2GS2: E722 read as E749 instead of E746). Mapping
  now uses the file's residue numbering to fill gaps and anchors the ends of every mapped stretch. Checked against
  SIFTS (PDBe's residue mapping) on about 60 chains: no residue gets a different position.
- **Comparisons dragged by one floppy end.** Superposing on every matched pair let a displaced tail pull the whole
  fit: EGFR's AlphaFold model vs 1M17 read as "246 of 312 residues deviate by more than 3 Å" (RMSD 11.7 Å) although
  the kinase domains agree within 1 Å. The fit now uses the largest common core and reports both the core RMSD and
  the all-pairs RMSD, so the curve shows where the structures really differ.
- **Wrong chain accepted.** Mapping HBB variants or comparing on 4HHB's α chains (a different protein) was allowed and
  gave meaningless positions; it is now refused with the chains that are HBB. Bad chain names in a comparison, the
  residue table or a residue query gave a server error (500) — now a plain message.
- **Structures that differ from UniProt** (an engineered G12C, pig insulin, an isoform) are now named per chain, and a
  variant mapped on such a residue says the structure carries a different amino acid there.
- **Identical repeats** (polyubiquitin UBB: 1UBQ fits three places equally) now get a note instead of silently picking one.
- **Gzipped uploads** did not show in the 3D viewer (it was handed the compressed file), and a second upload with
  the same name replaced the first one that comparisons still pointed to. Both fixed. The 200 MB limit now applies to
  the unpacked size: a 50 MB .gz of 3J3Q unpacked to 242 MB, took 4 GB of memory and a minute.
- **Empty or broken files**: a PDB file with only a header gave a server error; a truncated mmCIF gave "list index out
  of range". Both now explain themselves.
- **Eight identical look-ups at once** (double-click, several tabs, the assistant) built eight copies of the same
  protein; now one is built and the others wait for it.
- **Refresh without internet** silently showed the old cached files and said "Use Refresh to fetch the current
  versions". It now says the refresh failed and the copies may be out of date.
- **Proteins without a usual AlphaFold entry**: the SARS-CoV-2 spike's only model (a ColabFold model with an opaque id)
  was ignored as an "isoform"; it is now used and labelled as not from the DeepMind pipeline. For proteins with no
  model at all (PKHD1, titin) the page no longer says "the AlphaFold model is the only 3D view", and the reason given
  for titin (over 2,700 residues) is now accurate for human proteins too.
- **Isoform accessions** such as P02545-2 were quietly shown as the canonical protein; the page now says so and where
  the isoform model is.
- **The PDF report failed for titin**: its PDB reading listed ~280 domains, too tall for one page. Readings now
  describe at most 12 domains one by one and summarise the rest, and long findings break across pages.
- **Report file names** are now `YYYY-MM-DD_StructureBench-<gene>_report.pdf` (the day the analysis was created; for
  uploads a short name of the file), both from the server and on every report download link.
- Small things: "position 0" is no longer called "beyond the end"; the assistant's tools answer missing or mistyped
  arguments with a sentence instead of `KeyError`; the empty 3D viewer explains why when there is no model.
- Saved look-ups from 0.1.0 are rebuilt on the next look-up (schema 2), so they pick up the corrected mappings.

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
