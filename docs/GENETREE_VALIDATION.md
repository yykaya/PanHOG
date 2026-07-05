# Gene-tree validation of pangenome compartments (`--gene-trees`)

## Motivation

PanHOG's compartments are defined by **presence/absence frequency**: present in all
accessions = core, in exactly one = private, in between = shell. That count trusts the
upstream ortholog grouping
completely — it never asks whether the genes a HOG lumps together are a coherent,
divergent ortholog set. `--gene-trees` adds that check and assigns each candidate
HOG a phylogenetically-informed **status**, so the compartment becomes a supported
statement rather than a bare count.

## Method

**Non-private HOGs — two trees, compared.** For each HOG:

1. Extract member sequences with **accession-labelled tips** (each accession's
   FASTA is read separately, so shared gene IDs never collide).
2. **Protein tree** (default; peptides are always available): MAFFT protein
   alignment → RAxML-NG (`LG+G`).
3. **Codon tree** (optional, built only when `--cds` is given): the same protein
   alignment back-translated to a **codon (CDS) alignment** → RAxML-NG (`GTR+G`).
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

## How to read the output

`genetree_validation.tsv` has one row per candidate HOG. The columns to read together:

| Column | Meaning |
|---|---|
| `prot_support` / `codon_support` | mean bootstrap support of the protein / codon tree (higher = better resolved) |
| `codon_nRF` | normalised Robinson-Foulds distance of the codon tree vs the species tree (0 = identical topology, 1 = maximally different) |
| `prot~codon` | protein-vs-codon topology distance (0 = the two gene trees agree) |
| `Status` | Confirmed / Redundant / Conflict / Low-signal |

Interpreting them:

- **The codon tree usually adds signal.** It has ~3× the sites and captures synonymous
  variation, so a HOG that is *identical* at the protein level can still be resolved by
  the codon tree — bootstrap support typically rises relative to the protein tree. Use
  the codon tree as the arbiter when `--cds` is provided.
- **Protein and codon agree with each other but both conflict with the species tree**
  (high support, high `codon_nRF`, `prot~codon` ≈ 0) is the key red flag → `Conflict`:
  the members form a strongly-supported topology that is *not* the species topology
  (hidden paralogy / mis-grouping). That compartment call should be reviewed — the HOG
  is a candidate for splitting.
- **Private BLAST:** a "private" gene with a strong full-length hit in another
  accession (e.g. ≥ 95 % identity, near-complete coverage) is flagged
  **REVIEW: likely missed ortholog** — it is not truly private and should be merged
  with that accession's ortholog.
- **Single-gene resolution is limited at shallow divergence.** When accessions are
  closely related, a single short gene rarely recovers the species tree, so many
  core/shell HOGs are honestly `Low-signal` rather than a false `Confirmed`. That is
  the expected, correct behaviour; concatenation (`--supermatrix`) is the route to a
  resolved species-level tree.

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
