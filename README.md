# Structure Bench

A local app that takes a biologist from **gene → protein → structure**: type a gene, and Structure Bench finds
the protein in UniProt, fetches its AlphaFold model and every experimental structure, shows how confident each
part is, places your variants on it, and says in plain language what each panel means for *this* protein —
with a PDF report, PyMOL/ChimeraX scripts and a built-in Claude assistant. Sister app of RNAseq Bench and
MassSpec Bench; built on Biopython 1.88.

## Start

**Mac app** — unzip, run `xattr -cr` on the app once (see README-FIRST.txt), double-click. The first launch
installs packages into `~/Library/Application Support/StructureBench` (1–2 min). It serves
http://localhost:8767, alongside RNAseq Bench (8765) and MassSpec Bench (8766).

**From source**
```bash
./start.sh                    # creates .venv on first run, opens http://localhost:8767
./dev.sh                      # auto-reload while editing
.venv/bin/python -m server.selftest   # offline check of every analysis and assistant tool (~10 s)
.venv/bin/python -m pytest tests -q   # offline unit tests on the bundled LMNA example
```
To use the assistant, open **Settings & API key** and paste an Anthropic key (stored in
`~/.structure-bench/config.json`; if none is saved, RNAseq Bench's key or `ANTHROPIC_API_KEY` is used).

**Deep links** — `http://localhost:8767/#gene=LMNA&species=human` (species as a word or an NCBI taxon ID) and
`#uniprot=P02545` run the lookup straight away; `#item=<id>` reopens a saved one.

## What it does

**Look up** — gene symbol + organism (human, mouse, rat, zebrafish, fly, worm, yeast, E. coli, or a taxon ID)
or a UniProt accession; several reviewed hits → you choose.
- *UniProt*: accession, name, gene, organism, length, function, location, disease; features grouped into
  domains/regions, sites, PTMs and natural variants (with their disease notes).
- *AlphaFold*: canonical model (isoform models listed and loadable); per-residue pLDDT from the CA B-factor on a
  sequence track in the official colours, with UniProt domains and variant ticks on the same axis.
- *PAE*: heatmap with domain boundaries and hover readout, and a reading of whether the domains' relative
  positions are confident.
- *3D viewer* (3Dmol.js, vendored): colour by pLDDT/B-factor, chain, domain, relative SASA or AlphaMissense;
  cartoon/stick/surface; click a residue for its details and neighbours; PNG export.
- *Experimental structures*: every PDB entry UniProt lists, as coverage bars; **Load** downloads the mmCIF
  (cached), maps its chains to UniProt numbering, and makes it available for mapping and comparison.

**Map variants** — `R482W`, `p.R482W`, `p.Arg482Trp`, `Arg482Trp`, `G608G`, `p.Gly608=`, `R482*`. Each is checked
against the reference (with a suggestion when the numbering looks like the mature chain, skips Met1, or matches
an isoform), then described on the chosen structure: pLDDT or B-factor, relative SASA (buried < 20 %, exposed
> 50 %; alone and in the complex → interfaces with the partner named), residues within 5 Å (ligands, other
chains), secondary structure, UniProt domain and annotation, AlphaMissense, and a heuristic reading. Downloads:
the structure with the B-factor column replaced (variant flag or AlphaMissense), PyMOL `.pml`, ChimeraX `.cxc`.

**Compare structures** — AlphaFold vs a PDB entry, two PDB entries, or your own file: CE RMSD, per-residue CA
deviation after superposing the largest common core of the UniProt-matched pairs (so a flexible tail or a hinged
domain does not drag the fit), plotted with pLDDT, and a reading of where deviation
coincides with low confidence (expected) or not (a real difference). Superposed file downloadable; both shown
in the viewer.

**Upload** — `.pdb`, `.cif`, `.bcif` (optionally `.gz`), analysed fully offline with the file's own numbering. Up to 200 MB
once unpacked.

## Internet and honesty
Every download has a 20 s timeout and is cached in `StructureBench/cache`, so anything looked up once works
offline. Failures become a *NOTE* on the page ("UniProt could not be reached — showing the structure only"), never
a crash. The LMNA example is bundled (`examples/`). Things the app cannot do properly are said on the page:
without the DSSP program, secondary structure is a labelled phi/psi approximation; relative SASA on a single-chain
model ignores partners (coiled coils!); AlphaMissense is a prediction; an interface in a crystal may be a
crystal contact.

## Layout
```
server/app.py        HTTP API                    server/core.py        lookups, variants, comparison, exports
server/uniprot.py    UniProt search + entry      server/structure.py   Bio.PDB: parse, map, SASA, contacts, CE
server/alphafold.py  AFDB model, PAE, AM         server/variants.py    notation, reference check, heuristics
server/net.py        timeouts, cache, offline    server/report.py      PDF (reportlab + matplotlib)
server/agent.py      Claude assistant            server/selftest.py    offline self-test
web/index.html       the whole UI                web/vendor/           3Dmol.js 2.5.5 (BSD-3)
```
Saved analyses: `StructureBench/work/<id>/` (request.json, result.json, structures/, variants.json, report.pdf).

## Data and citation
UniProt and AlphaFold DB data are CC BY 4.0; PDB entries CC0. The Methods section of every result and the PDF list
the citations (Jumper 2021, Varadi 2024, Cheng 2023, UniProt 2025, Shrake & Rupley 1973, Tien 2013, Cock 2009,
Shindyalov & Bourne 1998). Created by Paul H. Kim, Ph.D.
