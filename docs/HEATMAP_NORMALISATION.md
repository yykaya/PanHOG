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

## Why some accessions look "higher" — annotation bias, not biology
When one or two accessions show a systematically higher mean copy number, more
multi-copy families and roughly double the gene count of the rest, the usual cause is
**annotation method, not biology**. An accession annotated by **reference lift-over**
inherits the reference's richer / more-fragmented gene models (and its gene IDs),
whereas **de-novo-annotated** accessions carry leaner models — so the lift-over
accession has more genes and a much higher multi-copy rate. Two lift-over accessions
built from the *same* reference will also correlate almost perfectly in copy number,
while correlating poorly with a divergent genome. The heatmap faithfully shows this;
it is an **annotation-method artefact**, not true copy-number variation. PanHOG prints
a `[WARNING]` when one accession's mean copy number exceeds 1.5× the median, exactly
to flag this.

## Normalisation options (`plot_genevar_heatmap`)
- **raw** (default): shows absolute copy number — dominated by the abundant accessions.
- **`--genevar` + log**: `log2(count+1)` — compresses the scale.
- **`--zscore`**: per-HOG z-score (`log2` then row-standardised) — removes each HOG's
  absolute level and shows **relative** variation across accessions. This is the right
  choice when accessions were annotated differently, because it de-emphasises the
  systematic per-accession offset introduced by lift-over annotation.

**Recommendation:** for cross-accession comparison use `--zscore`; and treat
copy-number differences between differently-annotated accessions (lift-over vs
de-novo) as suspect. No plot normalisation fixes a genuine 2× gene-count difference
that comes from the annotation pipeline — the fix is consistent annotation upstream.
