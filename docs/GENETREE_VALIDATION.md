# Gene-tree validation of pangenome compartments (`--gene-trees`)

## Motivation

PanHOG's compartments are defined by **presence/absence frequency**: 12/12 = core,
1/12 = private, 2–11/12 = shell. That count trusts the upstream ortholog grouping
completely — it never asks whether the genes a HOG lumps together are a coherent,
divergent ortholog set. `--gene-trees` adds that check and assigns each candidate
HOG a phylogenetically-informed **status**, so the compartment becomes a supported
statement rather than a bare count.

## Method

**Non-private HOGs — two trees, compared.** For each HOG:

1. Extract member sequences with **accession-labelled tips** (each accession's
   FASTA is read separately, so shared gene IDs never collide).
2. **Protein tree**: MAFFT protein alignment → RAxML-NG (`LG+G`).
3. **Codon tree**: the same protein alignment back-translated to a **codon (CDS)
   alignment** → RAxML-NG (`GTR+G`).
4. Score each tree (divergence, bootstrap support, normalised Robinson-Foulds vs
   the species tree), measure protein-vs-codon agreement, and assign a status.

The **codon tree is the primary arbiter**: it has ~3× the sites and captures
*synonymous* variation, so at the shallow divergence typical of a pangenome it is
better resolved than the protein tree (a gene can be identical in protein yet vary
synonymously). Status:

| Status | Meaning |
|---|---|
| **Confirmed** | codon tree well resolved (support ≥ 70) **and** concordant with the species tree → ortholog set is coherent; the compartment call holds. |
| **Redundant** | members near-identical even in synonymous sites → possible duplicate / over-merge. |
| **Conflict** | codon tree well supported but **conflicts** with the species tree → hidden paralogy / mis-grouping → compartment call suspect; candidate for HOG splitting. |
| **Low-signal** | neither tree resolved → cannot confirm (expected for conserved short genes over closely-related accessions). |

**Private HOGs — homology search.** A single gene has no tree, so the private gene
is **BLASTed** (blastp) against every other accession's proteome. A strong hit
elsewhere (identity ≥ 50 %, coverage ≥ 50 %) means the gene is **not truly
private** — its ortholog exists but was not grouped into the HOG — so the private
call is likely an ortholog-grouping artefact.

## Result on the toy data (12 *Arabidopsis* Chr1 accessions)

```
panhog --hog N0.tsv --fasta peptide/ --cds cds/ --pan \
       --species-tree species_tree.nwk --gene-trees
```

**Protein vs codon (non-private):**

| HOG | comp | prot support | codon support | codon nRF | prot~codon | Status |
|---|---|---|---|---|---|---|
| HOG0000060 | core | 39 | **52** | 1.0 | 0.44 | Low-signal |
| HOG0000261 | core | 1.5 | **14** | 1.0 | 1.0 | Low-signal |
| HOG0000026 | shell | 35 | 24 | 0.0 | 0.0 | Low-signal |
| HOG0000074 | shell | 75 | 75 | 1.0 | **0.0** | **Conflict** |

**Private (BLAST):**

| HOG | gene | best hit elsewhere | Status |
|---|---|---|---|
| HOG0000863 | Col-0 `AT1G06190.2` | **Nemrut-1, 99.6 % id, 95 % cov** | **REVIEW: likely missed ortholog** |

### What this tells us

- **The codon tree adds signal.** Support rises vs the protein tree (HOG0000060
  39→52), and HOG0000261 — *identical* at the protein level — shows tiny
  synonymous variation only the codon tree sees.
- **HOG0000074 is a real red flag.** Its protein and codon trees **agree with each
  other** (prot~codon nRF = 0) but **both conflict with the species tree at 75 %
  support** → the six "present" genes form a strongly-supported topology that is
  *not* the species topology: hidden paralogy / mis-grouping. Its "shell" call
  should be reviewed (candidate for splitting).
- **HOG0000863 is not really private.** Its Col-0 gene has a 99.6 %-identical
  homolog in Nemrut-1 that was not grouped into the HOG → the "private" label is an
  ortholog-calling miss, and the gene should be **reclassified** (merged with the
  Nemrut-1 ortholog).
- **Single-gene resolution is limited at this depth.** Most core/shell genes are
  `Low-signal` even with codon data — these 12 accessions are very closely related,
  so a single short gene rarely recovers the species tree. The tool correctly says
  "unresolved" rather than a false "Confirmed"; concatenation (`--supermatrix`) is
  the route to a resolved species-level tree.

## Reclassification

The status feeds compartment refinement:
- `Conflict` HOGs → flag for splitting; after correction the compartment may change.
- Private genes with a strong homolog elsewhere → **reclassify** (they are not
  private; merge with the homolog's HOG).
- `Redundant` core HOGs → check for assembly/annotation duplicates.

## How to run

```bash
conda activate panhog   # has mafft, raxml-ng, blast, ete3
panhog --hog N0.tsv --fasta peptide/ --cds cds/ --pan \
       --species-tree species_tree.nwk \
       --gene-trees --gene-trees-per-class 5      # or --gene-tree-hogs H1,H2,...
# -> results/compartments/<prefix>genetree_validation.tsv  + gene_trees/*.raxml.support
```

Omit `--cds` for a protein-only run; provide it to add the codon tree and the
protein-vs-codon comparison. Private validation needs BLAST+ (`makeblastdb`/`blastp`).

## Limitations / next steps

- **Multi-copy HOGs**: species-tree concordance is skipped (tips ≠ accessions); a
  duplication-aware reconciliation check fits here.
- **Thresholds** (support ≥ 70, nRF ≤ 0.5, identity/coverage ≥ 50 %) are defaults —
  tune per dataset / divergence.
- **Feedback loop**: automatically split `Conflict` HOGs and re-assign private genes
  with an ortholog elsewhere, then recompute compartments.
