# Deciding pangenome compartments: from proportion to phylogeny

Frequency alone (how many genomes carry a HOG) is a fast first pass but it trusts
the ortholog grouping and ignores biology. PanHOG now gives you three independent
lines of evidence to turn a raw count into a **confident** core / shell / private
call. Use them in this order.

## 1. Proportion (frequency) — the starting point
`--pan`: core = present in all genomes, private = 1 genome, shell = in between.
Fast, but it treats "present in 2 sister accessions" and "present in 2 deeply
divergent lineages" identically, and it inherits every ortholog-grouping error.

## 2. Weight the proportion by phylogeny — `--pan-weighted`
Classify by the **fraction of total tree branch length** a HOG's carriers span
(Faith's PD), not the raw count. A HOG spanning ≥ 90 % of the tree is "really"
core even if one tip is missing; one confined to a shallow cluster is shell/private
regardless of count. In practice this **promotes many near-complete HOGs** — present
in all-but-one or two tips, yet spanning almost the whole tree — from shell to core.
**Use PD-weighting when your accessions are unevenly sampled across the tree.**

## 3. Validate the ortholog set with gene trees — `--gene-trees`
Build a per-HOG **protein** tree (default) and, with `--cds`, a **codon** tree
(better resolved at shallow divergence). The status tells you whether to trust the
compartment:

| Status | What it means for the call |
|---|---|
| **Confirmed** | codon tree resolved + concordant with species tree → the ortholog set is coherent; keep the call. |
| **Conflict** | codon tree well-supported but conflicts with the species tree → hidden paralogy / mis-grouping → the count is unreliable; **split the HOG** and re-count. |
| **Redundant** | members near-identical even synonymously → duplicate/over-merge; keep but flag. |
| **Low-signal** | unresolved (normal for short, conserved genes at shallow divergence) → fall back to proportion + Ka/Ks; use `--supermatrix` concatenation for a species-level tree. |

## 4. Validate private genes with BLAST — (runs inside `--gene-trees`)
Private is the **most error-prone** compartment: a "private" call is often just an
ortholog OrthoFinder failed to group. Each private gene is BLASTed against every
other accession; a strong full-length hit elsewhere (≥ 50 % id, ≥ 50 % cov) means
it is **not really private** and should be merged and reclassified. In practice a
large share of "private" calls turn out to have a near-identical, full-length homolog
in another accession, so treat private counts as an **upper bound** until BLAST-checked.

## 5. Cross-check with selection — `--kaks-compartments`
Selection pressure is an **independent** axis that should track the compartment:
core genes sit under the strongest purifying selection (lowest Ka/Ks), dispensable
genes are the most relaxed (highest Ka/Ks). If your "core" shows high Ka/Ks, or
"dispensable" shows low Ka/Ks, the classification (or the ortholog grouping) is
suspect. Use per-HOG dN/dS as the **ratio of means** (the default) — never the mean
of per-pair ratios.

## A decision workflow

```
        frequency (--pan)
              │
     PD-weight (--pan-weighted)   ← uneven tree sampling? re-rank by branch length
              │
   gene-tree status (--gene-trees, +--cds)
        ├─ Confirmed  → keep
        ├─ Conflict   → split HOG, re-count  (compartment may change)
        ├─ Redundant  → keep + flag
        └─ Low-signal → rely on proportion + Ka/Ks
              │
   private?  → BLAST: homolog elsewhere → reclassify (not private)
              │
   Ka/Ks sanity (--kaks-compartments): core low, dispensable high?
        └─ pattern violated → revisit grouping/classification
              │
        final compartment
```

## Recommended thresholds (tune per dataset / divergence)
- Gene tree: bootstrap support ≥ 70 to trust topology; normalised RF ≤ 0.5 = concordant.
- PD-weighted: core ≥ 0.9, private ≤ 0.1 of total branch length.
- Private BLAST: identity ≥ 50 %, coverage ≥ 50 % ⇒ not private.
- Selection: dN/dS < 1 = purifying (expected for core); dN/dS > 1 = positive-selection candidate.

## Caveats
- At shallow pangenome divergence a single short gene rarely resolves — many HOGs
  are `Low-signal`; that is honest, not a failure. Concatenate single-copy core
  genes (`--supermatrix`) for a resolved species tree.
- Single-copy **private** genes have no pair, so their dN/dS is undefined; multi-copy
  private gives **paralog** dN/dS (a younger comparison than ortholog dN/dS) — do not
  pool it naively with core/shell ortholog dN/dS.

## Two data artefacts that distort the compartments (and the occupancy histogram)
The occupancy histogram (`occupancy_histogram.png`) can look "wrong" for reasons that
are in the **input**, not the plot:

1. **A distant outgroup deflates the core.** "core = present in ALL accessions" forces
   a gene to also exist in the outgroup, so any family the outgroup happens to lack
   drops out of core. Including a divergent outgroup can therefore shrink the core
   fraction substantially and lower frequency↔PD agreement. **Restrict the pangenome
   to the in-group with `--clade`** (or drop the outgroup with `--outgroup`) — otherwise
   "core" really means "conserved all the way to the outgroup". Removing the outgroup
   typically lifts both the core fraction and the agreement.
2. **Annotation twins inflate a mid-occupancy bar and suppress private.** If two
   accessions were annotated by **lift-over from the same reference** they share gene
   models and stay near-identical in copy number across most HOGs. That pair then
   dominates the "present-in-exactly-2" families — a spike at occupancy 2 that is an
   annotation-method artefact, not biology — and it steals genes that would otherwise
   be private. PanHOG prints a `[WARNING]` when one accession's mean copy number ≫ the
   others (see `docs/HEATMAP_NORMALISATION.md`). Fix upstream with consistent
   annotation, or treat such twins as a single accession.

The plotting itself is faithful: it simply counts, for each occupancy level, how many
gene families sit there, coloured by compartment.
