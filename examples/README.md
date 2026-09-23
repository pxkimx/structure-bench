# Bundled example — LMNA (human prelamin-A/C, UniProt P02545)

Used by the offline self-test, `tests/`, and the "LMNA example" button, and as the fallback when the network is
unavailable (see `lmna/manifest.json`: each key is the cache path the app would download to, each value the
bundled file). Example variants: E145K, R482W, R527P, G608G.

| File | What | Source | Licence / attribution |
|---|---|---|---|
| `lmna/uniprot_P02545.json` | UniProtKB entry P02545 (JSON) | https://rest.uniprot.org/uniprotkb/P02545.json | CC BY 4.0 — The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res* 53:D609–D617 (2025). |
| `lmna/uniprot_search_LMNA_human.json` | result of `gene_exact:LMNA AND organism_id:9606 AND reviewed:true` | UniProt REST search via `Bio.UniProt.search` | CC BY 4.0 — UniProt Consortium |
| `lmna/alphafold_predictions_P02545.json` | AlphaFold DB prediction records (canonical + isoforms) | https://alphafold.com/api/prediction/P02545 | CC BY 4.0 — AlphaFold DB (Google DeepMind / EMBL-EBI) |
| `lmna/AF-P02545-F1-model_v6.cif` | AlphaFold model, database version 6 | https://alphafold.ebi.ac.uk/files/AF-P02545-F1-model_v6.cif | CC BY 4.0 — Jumper J et al. *Nature* 596:583–589 (2021); Varadi M et al. *Nucleic Acids Res* 52:D368–D375 (2024) |
| `lmna/AF-P02545-F1-predicted_aligned_error_v6.json` | PAE matrix for the model | AlphaFold DB (paeDocUrl) | CC BY 4.0 — as above |
| `lmna/AF-P02545-F1-aa-substitutions.csv` | AlphaMissense scores for every missense substitution | AlphaFold DB (amAnnotationsUrl) | CC BY 4.0 — Cheng J et al. *Science* 381:eadg7492 (2023); as distributed by AlphaFold DB |

Downloaded 2026-09-23. These are copies for offline use; the live app fetches current versions (and caches them)
when the internet is available.
