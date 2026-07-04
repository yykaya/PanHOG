# Copy-number heatmap (`--genevar`) — why it differs from the summary, and normalisation

## The heatmap and the summary show *different quantities*
- `tf_summary_stats_plot.png` (`--summary`) = **total gene counts per compartment
  per accession** (and Gene/HOG ratios). One number per accession per compartment.
- `tf_genevar_heatmap.png` (`--genevar`) = **per-HOG copy number** (`dGeneNumbers`),
  one cell per HOG × accession, and it is **filtered to HOGs whose copy number varies
  by ≥ 2 across accessions** (the invariant ones are dropped). So it deliberately
  shows only the *variable / multi-copy* families — a different, non-comparable view.

They are not supposed to match: one is a per-compartment total, the other is a
per-family copy-number map of the most variable HOGs.

## Why Col-0 and Tanz-1 look "higher" — annotation bias, not biology
Measured on the input `N0.tsv`:

| | mean copies/HOG | % multi-copy | total genes |
|---|---|---|---|
| **Col-0** | 1.14 | 17.8 % | 12,564 |
| **Tanz-1** | 1.11 | 17.7 % | 12,278 |
| other 10 accessions | ~0.64 | ~1.3 % | ~7,000 |

And **Col-0 vs Tanz-1 have identical copy counts in 91.3 % of HOGs (r = 0.78)**, while
Col-0 vs A_lyrata r = 0.07. Tanz-1's annotation is a **reference lift-over of Col-0**
(both use TAIR `AT#G` IDs), so the two carry the same, richer/more-fragmented gene
models — ~2× the genes and ~14× the multi-copy rate of the de-novo-annotated
accessions. The heatmap faithfully shows this; it is an **annotation-method
artefact**, not true copy-number variation. PanHOG now prints a `[WARNING]` when one
accession's mean copy number exceeds 1.5× the median, exactly to flag this.

## Normalisation options (`plot_genevar_heatmap`)
- **raw** (default): shows absolute copy number — dominated by the abundant accessions.
- **`--genevar` + log**: `log2(count+1)` — compresses the scale.
- **`--zscore`**: per-HOG z-score (`log2` then row-standardised) — removes each HOG's
  absolute level and shows **relative** variation across accessions. This is the right
  choice when accessions were annotated differently, because it de-emphasises the
  systematic Col-0/Tanz-1 offset.

**Recommendation:** for cross-accession comparison use `--zscore`; and treat
copy-number differences between differently-annotated accessions (lift-over vs
de-novo) as suspect. No plot normalisation fixes a genuine 2× gene-count difference
that comes from the annotation pipeline — the fix is consistent annotation upstream.
