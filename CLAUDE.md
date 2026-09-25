# Structure Bench — notes for Claude Code

Local app: gene → protein (UniProt) → structure (AlphaFold DB, PDB), variant mapping and structure comparison.
FastAPI + Biopython 1.88 + numpy, single-file vanilla-JS UI, 3Dmol.js vendored, in-app Claude assistant, PDF report.
Sister of RNAseq Bench (:8765) and MassSpec Bench (:8766); this app is :8767 (Clone Bench :8768).
Owner: Paul (bench biologist, lamin/progeria background). Explain the *why* in plain language in every `how` /
`yours` text, and generate `yours` from the actual numbers.

## Layout / contracts
- `server/core.py` orchestrates *items* (one lookup or upload = `SB_HOME/work/<id>/`): `lookup`, `build_protein`,
  `load_pdb`, `load_isoform`, `upload`, `map_variants`, `residue_environment`, `compare_structures`,
  `export_structure`, `viewer_script`, `item_summary`, `methods`/`citations`. All take/return JSON-able dicts; the
  API, the assistant and the self-test call the same functions. No FastAPI imports outside `app.py`.
- `result.json` drives the UI, PDF and assistant. `structures[sid]` (sid = `AF`, `PDB-1IFR`, `UP-…`, isoform entry
  ids) → file, main_chain, chains (mapped to UniProt or not), bkind (`plddt`/`bfactor`). Per-residue tables live in
  `structures/<sid>_<chain>.residues.json` (key "A:482", pos = UniProt position, b, sasa, rsa, rsa_alone, ss).
- Bump `core.SCHEMA` when result.json changes shape: a saved lookup is reused only when its schema matches.
- `server/net.py`: every download → `fetch()` (20 s timeout, cache-first in `SB_HOME/cache`, falls back to
  `examples/*/manifest.json` seeds). `SB_OFFLINE=1` forces the offline path (tests, self-test). Biopython's own
  fetchers get the timeout by replacing `urlopen` in `Bio.UniProt` / `Bio.PDB.alphafold_db`.

## Gotchas (verified, not assumed)
- **Never `Entity.copy()` a structure you will transform or edit**: a copied `DisorderedAtom` keeps pointing at the
  original's selected altloc, so `transform`/`set_bfactor` silently miss atoms with alternate locations (showed up
  as 55 Å "deviations" on 7Z21). Re-parse the file instead (`st.load`) — see `superpose_pairs`, `export_structure`.
  Copies are fine for read-only use (SASA sub-models).
- Biopython 1.88 `BinaryCIFParser`: crashes on all-null string columns (AFDB .bcif `pdbx_PDB_ins_code`) → patched
  decoder in `structure._patch_binary_cif`; and it returns coordinates/B-factors as strings for some encoders →
  converted to numbers after parsing. Residue numbers can be `numpy.int32`: API responses go through `app.J()`.
- `PairwiseAligner`: Biopython ≥ 1.86 names end gaps `end_insertion_score` / `end_deletion_score`.
  `map_sequence(query, ref, nums=author numbers)`: unresolved stretches (numbering gaps ≤ 300) are filled with X
  placeholders before aligning, so a residue next to a gap lands at its own number (2GS2 E722 = EGFR E746, not E749);
  chain residues the alignment inserts count as mismatches in the 15-event ≥ 60 % identity window, and every kept run
  is trimmed until its end is identical, 4 of the next 5 are, and their spacing matches the numbering — or the end is
  numbered exactly as UniProt (4OBE's KRAS4B tail). That keeps fusion partners (2RH1/3SN6 T4 lysozyme), His tags and
  cloning scars (1IFR's GSH) off the protein. Checked once against SIFTS on ~60 chains (by hand, not an automated test).
- Chains that are not this protein (4HHB's α chains in an HBB item) are refused for mapping and comparison
  (`core._check_chain`) — they have no UniProt positions. Chain differences from UniProt (6OIM G12C, pig insulin) are
  listed per chain and flagged on the variant row; identical repeats (polyubiquitin, 1UBQ) get a note.
- Comparison superposes on the largest common core (`structure.superpose_pairs`: seeded on every 20-pair stretch,
  refined on pairs within 3 Å). A plain all-pairs fit put EGFR AF vs 1M17 at 12 Å because of a C-terminal tail.
- The DSSP binary is not installed here: `secondary_structure` uses DSSP only when `mkdssp` is on PATH, else a
  phi/psi approximation — always labelled as such in the UI, flags and methods.
- Several copies of the protein in one entry (1TUP): the main chain is the copy with most contacts to non-copies
  (DNA/partner/ligand). Users can map on any chain.
- 3Dmol: `requestAnimationFrame` does not fire in hidden tabs — cards draw with `setTimeout(…, 0)`.
- The UI's `.grid` must be `minmax(0,1fr)` or wide tables push the page past 375 px.
- Gzipped uploads: `structure.text_copy` writes `x.unzipped.cif` next to `x.cif.gz` (3Dmol.js and MMCIF2Dict cannot read
  gzip); `view_file` points at it. The 200 MB limit applies to the unpacked size (`core.MAX_STRUCTURE_BYTES`).
- `net.fetch` returns source `stale` when a refresh could not reach the network and fell back to the cache — the
  lookup must say so (it used to claim 'loaded from the cache, use Refresh').
- Report downloads are named `YYYY-MM-DD_StructureBench-<gene|file>_report.pdf` (`core.report_filename`, the item's
  creation date); the server sends it in Content-Disposition and the UI sets it as the `download` attribute.

## Run / test
- `./dev.sh` (auto-reload) · `.venv/bin/python -m pytest tests -q` (needs `pymupdf` for the report test — dev only, not in requirements) · `.venv/bin/python -m server.selftest`
- `./build_mac.sh` refreshes the app bundle (copies server, web, examples) and zips it; `macos/make_icon.py` redraws the icon.
- Bump `VERSION` + `CHANGELOG.md` for every change that ships; the launcher replaces an older running server.

- The report has two routes: `POST /api/items/<id>/report` (the page, with 3D snapshots) and `GET .../report.pdf`
  (assistant, direct download, and the page's fallback). Keep both: `web/` is served from disk while `server/` is the
  code loaded at startup, so rebuilding under a running app gives a new page on an old server — the page detects the
  missing route, falls back, and tells the user to restart.
