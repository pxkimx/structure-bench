Structure Bench — how to start
================================

1. Unzip, then move "Structure Bench.app" anywhere (Applications or Desktop).

2. First open only — macOS Gatekeeper will say the app is "damaged" or "cannot be verified"
   because it is not signed with a paid Apple developer certificate. Fix it once:
      open Terminal, type   xattr -cr    (with a space after it), drag the app into the
      Terminal window so its path is filled in, press Return.
   Then double-click the app. (Alternative: right-click the app -> Open -> Open.)

3. The first launch installs its packages into
   ~/Library/Application Support/StructureBench   (1–2 minutes, needs internet).
   Later launches take a few seconds. Your browser opens http://localhost:8767 and the
   server keeps running quietly in the background. To stop it: "Quit Structure Bench" in the
   sidebar, or double-click "Stop Structure Bench.command".
   It uses its own port (8767), so it runs side by side with RNAseq Bench (8765) and
   MassSpec Bench (8766).

4. Looking up a protein needs internet the first time (UniProt, AlphaFold DB, RCSB PDB).
   Everything fetched is cached in ~/Library/Application Support/StructureBench/cache, so a
   protein you have looked up once works offline afterwards. The LMNA example is bundled
   and always works offline.

If the app will not open at all: double-click "Start Structure Bench.command" instead —
it runs the same thing inside a Terminal window so you can see every message.
The log is at ~/Library/Application Support/StructureBench/structure-bench.log

"Run self-test.command" checks the installation offline in about 20 seconds.
Needs Python 3.10–3.13 (python.org recommended).
