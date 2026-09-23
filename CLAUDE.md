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
  `map_sequence` drops aligned pairs whose best 15-pair window is < 60 % identical and trims mismatching run ends
  (fusion partners like 6YJD's Gp7, His tags, cloning scars).
- The DSSP binary is not installed here: `secondary_structure` uses DSSP only when `mkdssp` is on PATH, else a
  phi/psi approximation — always labelled as such in the UI, flags and methods.
- Several copies of the protein in one entry (1TUP): the main chain is the copy with most contacts to non-copies
  (DNA/partner/ligand). Users can map on any chain.
- 3Dmol: `requestAnimationFrame` does not fire in hidden tabs — cards draw with `setTimeout(…, 0)`.
- The UI's `.grid` must be `minmax(0,1fr)` or wide tables push the page past 375 px.

## Run / test
- `./dev.sh` (auto-reload) · `.venv/bin/python -m pytest tests -q` · `.venv/bin/python -m server.selftest`
- `./build_mac.sh` refreshes the app bundle (copies server, web, examples) and zips it; `macos/make_icon.py` redraws the icon.
- Bump `VERSION` + `CHANGELOG.md` for every change that ships; the launcher replaces an older running server.
