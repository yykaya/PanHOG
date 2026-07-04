# The `phylogeny_weighted/` outputs — what each file means

`--species-tree` (+ `--pan-weighted`, `--gene-trees`) writes everything phylogeny-
related into `results/compartments/phylogeny_weighted/`. This is the phylogeny-aware
layer on top of the frequency-based `--pan` classification: it asks *where on the
species tree* each HOG sits, *how much of the tree* it spans, and *whether its own
gene tree is coherent*.

## Files

### `<prefix>hog_lca_analysis.tsv` — where each HOG originates
One row per HOG mapped onto the species tree.

| Column | Meaning |
|---|---|
| `HOG` | orthogroup id |
| `Num_Species` | number of carrier accessions |
| `LCA_Node` | the **lowest common ancestor** node (postorder-named `N0…Nk`) of all carriers — the deepest point on the species tree from which the HOG could descend |
| `PD_Fraction` | **Faith's phylogenetic diversity** the carriers span, as a fraction of the tree's total branch length (0–1). 1.0 = carriers reach every corner of the tree; small = confined to a shallow cluster |
| `Species_List` | the carrier accessions |

### `<prefix>hog_pd_weighted_class.tsv` — PD-weighted compartment
Re-classifies each HOG by how much of the **tree** it spans, not by a raw count.

| Column | Meaning |
|---|---|
| `HOG`, `Num_Species` | as above |
| `PD_Fraction` | Faith's PD fraction (0–1) |
| `Weighted_Class` | `core` if PD ≥ `--pan-weighted-core` (default 0.9), `private` if ≤ `--pan-weighted-private` (0.1), else `shell` |

Why it matters: a HOG present in 11/12 accessions but missing one *closely related*
tip still spans ~all the branch length → it is biologically core, though a flat
count calls it shell.

### `<prefix>classification_confidence.tsv` — frequency vs phylogeny (**confidence**)
Directly compares the two views per HOG.

| Column | Meaning |
|---|---|
| `Flat_Class` | frequency-based (core = all accessions, private = 1, else shell) |
| `Weighted_Class` | PD-weighted (above) |
| `Agreement` | `Same` or `Reclassified` |
| `Change` | e.g. `shell->core` when the two disagree |

**On this dataset (11,031 HOGs): 8,096 (73.4 %) agree; 2,935 (26.6 %) are
reclassified — almost all `shell → core` (2,886).** Those are HOGs the flat count
demotes to shell for missing an accession or two, but which span ≥90 % of the tree
and are phylogenetically core. Treat the **73 % that agree as high-confidence
compartment calls**; scrutinise the reclassified 27 % (mostly near-core HOGs).

### `<prefix>clade_compartments.tsv` — compartments *within* each clade
For every internal node, core/shell/private tallied among that clade's own tips
(`Clade_Core/Shell/Private`, `Num_Tips`, `HOGs_Present`) — lets you see, e.g., a
gene that is core within one lineage but absent from another.

### `<prefix>phylo_node_summary.tsv` — Dollo gain/loss per branch
`HOGs_LCA` (origins), `Gains`, `Losses` per node. Dollo parsimony places each HOG
once at its carriers' LCA and counts losses on branches to carrier-free clades.
(This dataset: 11,031 gains, 25,220 losses across 23 nodes.)

### `<prefix>species_tree_annotated.nwk`
The species tree with deterministic internal-node names (`N0…Nk`) so the `LCA_Node`
/ node-summary labels are identifiable in iTOL/FigTree.

### `<prefix>genetree_validation.tsv` + `gene_trees/`
Per-HOG **gene-tree** validation (see `docs/GENETREE_VALIDATION.md`): protein +
codon ML trees, scored for divergence, support and species-tree concordance, with a
`Status` (Confirmed / Redundant / Conflict / Low-signal). `Conflict` = the HOG's own
gene tree is well-supported but disagrees with the species tree → the members are
probably **not clean orthologs** and the compartment call is suspect. The
`gene_trees/` folder holds each `.raxml.support` tree plus example plots
(`genetree_core_example.png`, `genetree_shell_example.png`) to eyeball clustering.

## How to decide "confident" compartments
1. **Agree in `classification_confidence.tsv`** (73 %) → trust the call.
2. **Reclassified** (27 %) → prefer `Weighted_Class` when accessions are unevenly
   sampled on the tree (mostly near-core shell→core here).
3. **Gene-tree `Conflict`** → flag/split regardless of count.
4. **Private** → BLAST-validate (most are missed orthologs; see the guidance doc).

See also `docs/CLASSIFICATION_GUIDANCE.md`.
