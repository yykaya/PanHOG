# PanHOG v0.3.0 — Local Upgrade Guide

This guide explains how to move your older local PanHOG to the upgraded
**v0.3.0** code, set it up on your Mac, and then keep improving it in a
**local** working environment (with access to your real data, MAFFT,
PAML `codeml`, and conda).

---

## 1. What changed in v0.3.0

- **Ka/Ks (dN/dS) fixed.** The old built-in Nei–Gojobori routine counted
  synonymous/non-synonymous *sites* with arbitrary constants and only on
  differing codons, so the numbers were not interpretable. It is now delegated
  to validated engines and undefined values are reported as `NaN`.
- **New dN/dS engines** via `--kaks-method`: `biopython` (models NG86/LWL85/
  YN00/ML, set with `--kaks-model`), **`codeml` (PAML pairwise)**, and
  `kakscalculator`. Writes a per-HOG summary *and* a per-pair table.
- **Modern BioPython fix.** On BioPython ≥ 1.85 a removed import was silently
  disabling *all* BioPython features; fixed.
- **Phylogeny-aware analysis** (`--species-tree`): deterministic node naming +
  annotated Newick tree, species/tip validation, a Faith's PD fraction per HOG,
  and **Dollo-parsimony gain/loss** with per-branch counts.
- **Packaging:** ships `panhog_dnds` and `panhog_phylo`, unified conda recipe,
  `pyproject.toml`, `.gitignore`, and a `pytest` suite (20 tests).

Full details are in `README.md` (Changelog).

---

## 2. Folder layout on your Mac

```
/Users/yasin/Desktop/Update_PanHOG/
├── 19Nov25/                        # your older PanHOG (may hold local-only data/configs)
├── 19Nov25_improved/               # v0.3.0 (this branch) — the upgraded code
└── dN-dS-Snakemake-pipeline-main/  # the PAML reference pipeline
```

---

## 3. Get the improved code into your workflow

You already cloned `19Nov25_improved`. Pick one path:

**Option A — Use the clone as your new working copy (recommended).**
Just work in `19Nov25_improved`. Copy over any local-only inputs (your real
`N0.tsv`, peptide/CDS FASTA, configs, notes) that live only in `19Nov25`.
Prompt 1 below does this diff for you safely.

**Option B — Pull the branch into an existing git repo.**
If `19Nov25` is already a git clone of the same GitHub repo, fetch and check out
the v0.3.0 branch you pushed:
```bash
cd /Users/yasin/Desktop/Update_PanHOG/19Nov25
git fetch origin <your-v0.3.0-branch>
git checkout <your-v0.3.0-branch>
```

**Option C — Copy files (last resort, no git).**
Copy `PanHOG.py`, `panhog_dnds.py`, `panhog_phylo.py`, `setup.py`,
`pyproject.toml`, `conda-recipe/`, `environment.yml`, `tests/` from
`19Nov25_improved` into your working folder, overwriting the old files.

> To refresh the clone later: `cd 19Nov25_improved && git pull`.

---

## 4. One-time setup (do this once)

```bash
cd /Users/yasin/Desktop/Update_PanHOG/19Nov25_improved

# Full toolchain (mafft, muscle, pal2nal, paml, blast, scipy, ...)
conda env create -f environment.yml
conda activate panhog

# Editable install with test/optional deps (scipy enables YN00/ML dN/dS models)
pip install -e ".[dev]"

# Verify: all 20 tests should pass
pytest -q
```

---

## 5. Running the new features

**dN/dS with the PAML `codeml` engine (the pipeline's approach):**
```bash
panhog --hog <N0.tsv> --fasta <peptides_dir> --cds <cds_dir> --pan \
  --kaks --kaks-type core --kaks-method codeml \
  -o results_codeml/ -p cml_
```

**dN/dS with the corrected built-in engine (no external tools):**
```bash
panhog --hog <N0.tsv> --fasta <peptides_dir> --cds <cds_dir> --pan \
  --kaks --kaks-type core --kaks-method biopython --kaks-model NG86 \
  -o results_bp/ -p bp_
# Outputs: bp_kaks_results_core.tsv (per-HOG) + bp_kaks_pairwise_core.tsv (per-pair)
```

**Phylogeny-aware analysis (LCA + PD fraction + Dollo gain/loss):**
```bash
panhog --hog <N0.tsv> --fasta <peptides_dir> --pan \
  --species-tree <species_tree.nwk> -o results_phylo/ -p phy_
# Outputs: phy_hog_lca_analysis.tsv, phy_hog_gainloss.tsv,
#          phy_phylo_node_summary.tsv, phy_species_tree_annotated.nwk
```

---

## 6. Recommended local follow-up tasks

Run these from inside `19Nov25_improved`, one at a time.
Replace `<...>` placeholders with your real paths.

### Task 0 — Orient yourself (read first)
```
I'm working on PanHOG, a phylogeny-aware pangenome/HOG toolkit. This folder is
v0.3.0. Read README.md, PanHOG.py, panhog_dnds.py, panhog_phylo.py and the
tests/ folder, then give me a short summary of the architecture and the dN/dS
and phylogeny engines. Don't change anything yet.
```

### Prompt 1 — Bring my older local files across, safely
```
Two sibling folders exist: ../19Nov25 (my older PanHOG, which may contain
local-only data, configs or notes that were never pushed to GitHub) and this
folder (upgraded v0.3.0). Compare them. List every file that exists only in
../19Nov25 or differs, and classify each as (a) my local data/config/notes
worth keeping or (b) old code now superseded by the upgrade. Do NOT overwrite
anything yet — propose a plan to copy only my local-only inputs into this
folder, then wait for my approval.
```

### Prompt 2 — Environment + tests
```
Set up the environment and run the tests. Create and activate the conda env
from environment.yml, then `pip install -e ".[dev]"`, then run `pytest -q`.
Report how many of the 20 tests pass. If conda is slow or fails, offer a
pip-only fallback (numpy pandas biopython scipy pytest) that still runs the
unit tests, and tell me which tests need the external tools.
```

### Task 3 — Validate the PAML `codeml` engine (not yet run locally)
```
Verify the PAML codeml dN/dS engine end-to-end; this code path was written but
never executed against a real codeml binary. Confirm `codeml` is on PATH
(install with `conda install -c bioconda paml` if missing). Then run, on my
data (HOG=<N0.tsv>, peptides=<peptides_dir>, CDS=<cds_dir>):
  panhog --hog <N0.tsv> --fasta <peptides_dir> --cds <cds_dir> --pan \
    --kaks --kaks-type core --kaks-method codeml -o codeml_test/ -p cml_
Then run the same with `--kaks-method biopython --kaks-model NG86` into
bp_test/. Compare per-HOG dN/dS between the two engines, flag HOGs where they
disagree strongly, and confirm codeml ran cleanly. Show me the head of both
kaks_pairwise_core.tsv files.
```

### Prompt 3b — Cross-check my codeml setup against the real pipeline
```
Read the reference pipeline at ../dN-dS-Snakemake-pipeline-main (its Snakefile
and scripts). It computes pairwise dN/dS with PAML. Compare its exact codeml
control-file settings (runmode, seqtype, model, NSsites, CodonFreq, icode,
cleandata, kappa/omega handling) against my dnds_pair_codeml() in
panhog_dnds.py. List any differences that would change the dN/dS numbers, and
whether the pipeline does something important I'm missing (a specific
ortholog-calling or filtering step). Propose the minimal patch to align them
and show me a diff before applying anything.
```

### Prompt 4 — Real full run on my dataset
```
Run the full v0.3.0 pipeline on my real data and sanity-check it. Inputs:
HOG=<N0.tsv>, peptides=<peptides_dir>, CDS=<cds_dir>, species tree=<tree.nwk>.
Run pangenome classification, Ka/Ks on core (biopython NG86), and the
phylogeny-aware analysis. Then summarize: counts of core/shell/private HOGs,
the dN/dS distribution (median and % under purifying selection, dN/dS < 1), and
the Dollo gain/loss totals per major node from phylo_node_summary.tsv. Point out
anything suspicious — e.g. species in N0.tsv missing from the tree, or empty
alignments.
```

### Prompt 5 — Continue the roadmap (Tier 2.1 / 2.2)
```
Extend panhog_phylo.py with the next roadmap items, matching the existing code
style and adding pytest tests in tests/ (use a small hand-worked tree like
tests/test_phylo.py). (1) PD-weighted core/shell/private classification: classify
each HOG by the fraction of total tree branch length its carriers span, with
configurable thresholds, exposed as a new `--pan-weighted` flag. (2) Clade-
conditioned compartments: for each internal node, compute core/shell/private
within that clade's subtree, written to a per-clade table. Run pytest, update
the README, and commit on a new branch (don't touch main).
```

### Prompt 6 — Scope Tier 3 reconciliation (for the meeting with my collaborator)
```
Help me scope gene-tree/species-tree reconciliation (Tier 3). My OrthoFinder run
has per-HOG gene trees in <Gene_Trees_dir>. Draft a design for a new
panhog_reconcile.py (design only, don't implement) that reconciles each HOG gene
tree against the species tree to get per-HOG duplication/loss (and, for
prokaryotes, transfer) counts using an external tool. Evaluate GeneRax, Notung
and ecceTERA, recommend one and say why. Specify inputs, outputs, new CLI flags,
and a test strategy. Write it as a plan I can review with my phylogenetics
collaborator before we build.
```

### Prompt 7 — Build and smoke-test the conda package
```
Build the conda package and smoke-test it. Run `conda build conda-recipe/`
(use mamba if available), fix any recipe errors, install the built package into
a fresh env, and confirm `panhog --help`, `pangenehog --help`, and
`python -c "import panhog_dnds, panhog_phylo"` all work. Report the built
artifact path and any changes you made to the recipe.
```

---

## 7. Tips for local development

- Work **inside** `19Nov25_improved` so your tooling sees this code, the tests,
  and (via `../`) both your older folder and the reference pipeline.
- Work on a branch for new features (Prompts 5–6) and review diffs before
  committing; keep `main` clean.
- Tasks 3 and 3b are the highest priority — they close the one gap not yet
  tested locally: that PAML `codeml` runs correctly and matches your
  reference pipeline.
- Re-run `pytest` after any change; the suite is fast (< 1s).
