#!/usr/bin/env python3
"""
panhog_genetrees.py
===================
Phylogenetic validation of pangenome compartment calls.

PanHOG's core / shell / private classification is based purely on the
presence/absence *frequency* of a HOG across genomes (12/12 = core, 1/12 =
private, in between = shell). That is fast, but it trusts the upstream ortholog
grouping completely. This module double-checks each candidate HOG with gene
trees and, for private genes, a homology search, then assigns a phylogenetically
informed *status* to the frequency-based call.

For non-private HOGs it builds two ML trees per HOG and compares them:

  * protein tree : MAFFT protein alignment -> RAxML-NG (LG+G).
  * codon  tree  : same protein alignment back-translated to a codon (CDS)
                   alignment -> RAxML-NG (GTR+G).

The codon tree has ~3x the sites and captures *synonymous* variation, so at the
shallow divergence typical of a pangenome it is far better resolved than the
protein tree — it is the primary arbiter of whether the grouped genes are a
coherent, divergent ortholog set that is topologically consistent with the
species tree. Each HOG gets a status:

  Confirmed  – codon tree well resolved and concordant with the species tree.
  Redundant  – members near-identical even in synonymous sites (possible
               duplicate / over-merge).
  Conflict   – codon tree well supported but conflicts with the species tree
               (hidden paralogy / mis-grouping) -> the compartment call is
               suspect and the HOG may need splitting.
  Low-signal – neither tree resolved (expected for very conserved short genes).

For private HOGs (one accession, no tree) it BLASTs the gene against every other
accession's proteome: a strong hit elsewhere means the "private" call is likely a
missed ortholog, not a genuinely accession-specific gene.

Tip labels use the *accession* name (each accession's FASTA is read separately),
so identical gene IDs shared across accessions never collide.

External tools (optional; absence is reported, never faked): MAFFT, RAxML-NG,
and (for private validation) BLAST+ (makeblastdb / blastp). Tree analysis uses
ete3 when available.
"""

from __future__ import annotations

import csv
import math
import os
import shutil
import statistics
import subprocess
from itertools import combinations

try:
    from ete3 import Tree as EteTree
    HAS_ETE3 = True
except Exception:  # pragma: no cover - ete3 optional
    HAS_ETE3 = False


# ---------------------------------------------------------------------------
# HOG occupancy + candidate selection
# ---------------------------------------------------------------------------

def parse_hog_occupancy(hogsfile):
    """
    Parse an OrthoFinder ``N0.tsv`` into (species_list, occupancy).

    ``occupancy`` is ``{hog_id: {species: [gene_ids]}}`` holding only the
    species that carry the HOG.
    """
    with open(hogsfile) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        species = header[3:]
        occ = {}
        for row in reader:
            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))
            genes = {}
            for i, sp in enumerate(species):
                cell = row[3 + i].strip()
                if cell:
                    genes[sp] = [g.strip() for g in cell.split(",") if g.strip()]
            if genes:
                occ[row[0]] = genes
    return species, occ


def compartment_of(occ_hog, n_species):
    """Frequency-based compartment for one HOG's occupancy dict."""
    npres = len(occ_hog)
    if npres >= n_species:
        return "core"
    if npres == 1:
        return "private"
    return "shell"


def select_candidates(occ, n_species, per_class=3, single_copy_only=True):
    """Pick a few example HOGs per compartment for validation."""
    buckets = {"core": [], "shell": [], "private": []}
    for hog, genes in occ.items():
        if single_copy_only and any(len(v) != 1 for v in genes.values()):
            continue
        c = compartment_of(genes, n_species)
        if len(buckets[c]) < per_class:
            buckets[c].append(hog)
        if all(len(buckets[k]) >= per_class for k in buckets):
            break
    return buckets


# ---------------------------------------------------------------------------
# Per-accession sequence store (collision-safe, protein or CDS)
# ---------------------------------------------------------------------------

class SeqStore:
    """
    Lazily read per-accession FASTA files, keyed by a *normalised* accession
    name so the same store maps a HOG's accession (an N0 column such as
    ``Col-0_Chr1.pep``) to either the peptide file ``Col-0_Chr1.pep.fa`` or the
    CDS file ``Col-0_Chr1.cds.fa``. Reading each accession separately means
    identical gene IDs in different accessions never overwrite each other.
    """

    _EXTS = (".fasta", ".faa", ".fna", ".fa")
    _TAGS = (".pep", ".cds", ".prot", ".aa", ".nt")

    @classmethod
    def _norm(cls, key):
        for ext in cls._EXTS:
            if key.lower().endswith(ext):
                key = key[: -len(ext)]
                break
        for tag in cls._TAGS:
            if key.lower().endswith(tag):
                key = key[: -len(tag)]
                break
        return key

    def __init__(self, fasta_dir):
        self.fasta_dir = fasta_dir
        self._file_of = {}
        for fn in os.listdir(fasta_dir):
            self._file_of[self._norm(fn)] = os.path.join(fasta_dir, fn)
        self._cache = {}

    def _load(self, accession):
        acc = self._norm(accession)
        if acc in self._cache:
            return self._cache[acc]
        path = self._file_of.get(acc)
        seqs = {}
        if path and os.path.exists(path):
            cur, buf = None, []
            with open(path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if cur:
                            seqs[cur] = "".join(buf)
                        cur = line[1:].split("|")[0].split()[0]
                        buf = []
                    else:
                        buf.append(line.strip())
            if cur:
                seqs[cur] = "".join(buf)
        self._cache[acc] = seqs
        return seqs

    def get(self, accession, gene):
        return self._load(accession).get(gene)


# ---------------------------------------------------------------------------
# FASTA -> alignment -> back-translation -> tree
# ---------------------------------------------------------------------------

def _safe_tip(name):
    for ch in "():,; '\t":
        name = name.replace(ch, "_")
    return name


def hog_tips(occ_hog):
    """Yield (tip_label, accession, gene) for a HOG (single- or multi-copy)."""
    single_copy = all(len(v) == 1 for v in occ_hog.values())
    for acc, genes in occ_hog.items():
        for g in genes:
            tip = _safe_tip(acc if single_copy else f"{acc}__{g}")
            yield tip, acc, g


def write_protein_fasta(occ_hog, store, path):
    """Write accession-labelled protein sequences. Returns {tip: (acc, gene)}."""
    tip_map = {}
    with open(path, "w") as out:
        for tip, acc, gene in hog_tips(occ_hog):
            seq = store.get(acc, gene)
            if seq:
                out.write(f">{tip}\n{seq}\n")
                tip_map[tip] = (acc, gene)
    return tip_map


def read_fasta(path):
    """Parse a FASTA/alignment into an ordered list of (id, seq)."""
    recs, cur, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    recs.append((cur, "".join(buf)))
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        recs.append((cur, "".join(buf)))
    return recs


def backtranslate(prot_aln_recs, tip_map, cds_store):
    """
    Turn a protein alignment into a codon alignment using each tip's CDS.

    Every aligned residue maps to its 3 CDS nucleotides; gaps become ``---``.
    Returns list of (tip, codon_seq), skipping tips with no CDS. Length checks
    are lenient (missing bases -> gap) so partial/odd CDS never crash the run.
    """
    out = []
    for tip, aap in prot_aln_recs:
        acc, gene = tip_map.get(tip, (None, None))
        cds = cds_store.get(acc, gene) if acc else None
        if not cds:
            continue
        codons, i = [], 0
        for aa in aap:
            if aa == "-":
                codons.append("---")
            else:
                codon = cds[i:i + 3]
                codons.append(codon if len(codon) == 3 else "---")
                i += 3
        out.append((tip, "".join(codons)))
    return out


def align_mafft(in_fasta, out_aln, mafft_path="mafft"):
    if shutil.which(mafft_path) is None:
        raise FileNotFoundError(f"MAFFT not found (looked for '{mafft_path}').")
    with open(out_aln, "w") as out:
        r = subprocess.run([mafft_path, "--auto", "--quiet", in_fasta],
                           stdout=out, stderr=subprocess.DEVNULL)
    return r.returncode == 0 and os.path.exists(out_aln) and os.path.getsize(out_aln) > 0


def build_tree_raxml(aln, prefix, raxml_path="raxml-ng", model="LG+G",
                     bs_trees=100, seed=42, threads=1):
    """
    Run RAxML-NG (``--all``: ML search + bootstrap + support).
    Returns (tree_path, status). tree_path is the support tree or None.
    """
    if shutil.which(raxml_path) is None:
        raise FileNotFoundError(
            f"RAxML-NG not found (looked for '{raxml_path}'). "
            f"Install e.g. `conda install -c bioconda raxml-ng`."
        )
    cmd = [raxml_path, "--all", "--msa", aln, "--model", model,
           "--bs-trees", str(bs_trees), "--prefix", prefix,
           "--threads", str(threads), "--seed", str(seed), "--force", "--redo"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    support = f"{prefix}.raxml.support"
    best = f"{prefix}.raxml.bestTree"
    if os.path.exists(support):
        return support, "ok"
    if os.path.exists(best):
        return best, "ok (no support)"
    tail = (r.stderr or r.stdout or "").strip().splitlines()
    return None, "raxml failed: " + (tail[-1] if tail else "unknown")


# ---------------------------------------------------------------------------
# Tree QC metrics + concordance
# ---------------------------------------------------------------------------

def tree_metrics(tree_path):
    """n_tips, total tree length, mean pairwise patristic distance, mean support."""
    if not HAS_ETE3 or not tree_path:
        return {}
    t = EteTree(tree_path)
    leaves = t.get_leaves()
    total_len = sum(nd.dist for nd in t.traverse())
    dists = [a.get_distance(b) for a, b in combinations(leaves, 2)]
    sups = [nd.support for nd in t.traverse()
            if not nd.is_leaf() and nd.support is not None]
    return {
        "n_tips": len(leaves),
        "tree_len": round(total_len, 4),
        "mean_patristic": round(statistics.mean(dists), 4) if dists else 0.0,
        "mean_support": round(statistics.mean(sups), 1) if sups else float("nan"),
    }


def _nrf(tree_a_path, tree_b_path, fmt_b=0):
    if not HAS_ETE3 or not tree_a_path or not tree_b_path:
        return None
    a = EteTree(tree_a_path)
    try:
        b = EteTree(tree_b_path, format=fmt_b)
    except Exception:
        b = EteTree(tree_b_path)
    common = set(a.get_leaf_names()) & set(b.get_leaf_names())
    if len(common) < 4:
        return None
    a.prune(common, preserve_branch_length=True)
    b.prune(common, preserve_branch_length=True)
    rf = a.robinson_foulds(b, unrooted_trees=True)
    return round(rf[0] / rf[1], 3) if rf[1] else 0.0


def species_concordance(gene_tree_path, species_tree_path):
    """Normalised Robinson-Foulds vs the species tree (0 = identical topology)."""
    return _nrf(gene_tree_path, species_tree_path, fmt_b=1)


def tree_agreement(tree_a_path, tree_b_path):
    """Normalised RF between two gene trees (e.g. protein vs codon)."""
    return _nrf(tree_a_path, tree_b_path, fmt_b=0)


def make_verdict(metrics, nrf, min_divergence=0.005, max_nrf=0.5, min_support=70.0):
    """OK / REVIEW verdict for a single tree (kept for one-tree use + tests)."""
    n = metrics.get("n_tips", 0)
    div = metrics.get("mean_patristic", 1.0)
    sup = metrics.get("mean_support")
    has_sup = sup is not None and not (isinstance(sup, float) and math.isnan(sup))
    flags, notes = [], []
    if n >= 4 and div < min_divergence:
        flags.append("near-identical sequences")
    if nrf is not None and nrf > max_nrf:
        if has_sup and sup >= min_support:
            flags.append(f"well-supported species-tree conflict (nRF={nrf}, support={sup})")
        else:
            notes.append(f"low phylogenetic signal (nRF={nrf}, support={sup if has_sup else 'NA'})")
    if flags:
        return "REVIEW: " + "; ".join(flags)
    if notes:
        return "OK — " + "; ".join(notes)
    return "OK"


def compute_status(prot_m, prot_nrf, codon_m, codon_nrf,
                   min_div=0.005, min_support=70.0, max_nrf=0.5):
    """
    Assign a compartment-validation status, using the codon tree as the primary
    arbiter (more sites / synonymous signal) and falling back to the protein
    tree. Returns (status, reason).
    """
    if codon_m:
        m, nrf, arbiter = codon_m, codon_nrf, "codon"
    elif prot_m:
        m, nrf, arbiter = prot_m, prot_nrf, "protein"
    else:
        return "no-tree", "no tree could be built"

    n = m.get("n_tips", 0)
    div = m.get("mean_patristic", 1.0)
    sup = m.get("mean_support")
    has_sup = sup is not None and not (isinstance(sup, float) and math.isnan(sup))
    if n < 4:
        return "Low-signal", "too few tips for a tree"
    if div < min_div:
        return "Redundant", f"near-identical even in the {arbiter} tree (patristic={div})"
    if has_sup and sup >= min_support:
        if nrf is not None and nrf > max_nrf:
            return "Conflict", (f"{arbiter} tree well supported (support={sup}) but conflicts "
                                f"with the species tree (nRF={nrf}) — possible mis-grouping")
        return "Confirmed", (f"{arbiter} tree resolved (support={sup})"
                             + (f" & concordant (nRF={nrf})" if nrf is not None else ""))
    return "Low-signal", f"{arbiter} tree unresolved (support={sup if has_sup else 'NA'})"


# ---------------------------------------------------------------------------
# Private-gene homology validation (BLAST)
# ---------------------------------------------------------------------------

def build_proteome_db(species, store, workdir, makeblastdb_path="makeblastdb"):
    """
    Concatenate every accession's proteome with ``accession|gene`` headers into a
    single BLAST protein database. Returns the db path, or None on failure.
    """
    if shutil.which(makeblastdb_path) is None:
        return None
    combined = os.path.join(workdir, "all_proteomes.faa")
    with open(combined, "w") as out:
        for acc in species:
            for gene, seq in store._load(acc).items():
                if seq:
                    out.write(f">{acc}|{gene}\n{seq}\n")
    db = os.path.join(workdir, "all_proteomes_db")
    r = subprocess.run([makeblastdb_path, "-in", combined, "-dbtype", "prot",
                        "-out", db], capture_output=True, text=True)
    return db if r.returncode == 0 else None


def blast_best_other_hit(gene_seq, self_accession, db, workdir,
                         blastp_path="blastp", evalue="1e-3"):
    """
    BLAST one protein against the all-accession DB and return the best hit from a
    *different* accession as {acc, pident, qcov, evalue}, or None.
    """
    if shutil.which(blastp_path) is None or not db:
        return None
    q = os.path.join(workdir, "query.faa")
    with open(q, "w") as fh:
        fh.write(f">query\n{gene_seq}\n")
    out = os.path.join(workdir, "blast.tsv")
    with open(out, "w") as fh:
        subprocess.run(
            [blastp_path, "-query", q, "-db", db, "-evalue", evalue,
             "-max_target_seqs", "50",
             "-outfmt", "6 sseqid pident length evalue qcovs"],
            stdout=fh, stderr=subprocess.DEVNULL)
    best = None
    with open(out) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            acc = f[0].split("|")[0]
            if acc == SeqStore._norm(self_accession) or acc == self_accession:
                continue
            hit = {"acc": acc, "pident": float(f[1]),
                   "evalue": float(f[3]), "qcov": float(f[4])}
            if best is None or hit["pident"] > best["pident"]:
                best = hit
    return best


def validate_private(hog, occ_hog, species, store, db, workdir,
                     blastp_path="blastp", min_pident=50.0, min_qcov=50.0):
    """BLAST a private gene against all other accessions to see why it is private."""
    acc = next(iter(occ_hog))
    gene = occ_hog[acc][0]
    seq = store.get(acc, gene)
    row = {"HOG": hog, "Compartment": "private", "N_species": 1,
           "N_genes": len(occ_hog[acc]), "Accession": acc, "Gene": gene}
    if not seq or not db:
        row["Status"] = "unchecked"
        row["Reason"] = "no sequence or BLAST DB unavailable"
        return row
    best = blast_best_other_hit(seq, acc, db, workdir, blastp_path)
    if best is None:
        row["Status"] = "Confirmed private"
        row["Reason"] = "no homolog found in any other accession (genuinely accession-specific)"
    else:
        row["BLAST_hit"] = best["acc"]
        row["BLAST_pident"] = round(best["pident"], 1)
        row["BLAST_qcov"] = round(best["qcov"], 1)
        if best["pident"] >= min_pident and best["qcov"] >= min_qcov:
            row["Status"] = "REVIEW: likely missed ortholog"
            row["Reason"] = (f"strong homolog in {best['acc']} "
                             f"(id={best['pident']:.0f}%, cov={best['qcov']:.0f}%) — "
                             f"not truly private; may be an ortholog-grouping miss")
        else:
            row["Status"] = "Confirmed private"
            row["Reason"] = (f"only weak homology elsewhere (best {best['acc']} "
                             f"id={best['pident']:.0f}%, cov={best['qcov']:.0f}%)")
    return row


# ---------------------------------------------------------------------------
# Per-HOG gene-tree validation (protein + codon)
# ---------------------------------------------------------------------------

def validate_hog(hog, occ_hog, n_species, prot_store, cds_store, workdir,
                 species_tree=None, mafft_path="mafft", raxml_path="raxml-ng",
                 prot_model="LG+G", codon_model="GTR+G", bs_trees=100):
    """Build protein and (if CDS available) codon trees for one HOG and score them."""
    row = {"HOG": hog, "Compartment": compartment_of(occ_hog, n_species),
           "N_species": len(occ_hog),
           "N_genes": sum(len(v) for v in occ_hog.values())}

    prot_fa = os.path.join(workdir, f"{hog}.faa")
    tip_map = write_protein_fasta(occ_hog, prot_store, prot_fa)
    if len(tip_map) < 4:
        row["Status"] = "Low-signal"
        row["Reason"] = f"only {len(tip_map)} sequences (too few for a tree)"
        return row

    prot_aln = os.path.join(workdir, f"{hog}.prot.aln")
    if not align_mafft(prot_fa, prot_aln, mafft_path):
        row["Status"] = "no-tree"
        row["Reason"] = "protein alignment failed"
        return row

    # Protein tree
    prot_tree, _ = build_tree_raxml(prot_aln, os.path.join(workdir, f"{hog}.prot"),
                                    raxml_path, prot_model, bs_trees)
    prot_m = tree_metrics(prot_tree)
    prot_nrf = species_concordance(prot_tree, species_tree) if species_tree else None

    # Codon tree (protein alignment back-translated with the CDS)
    codon_tree, codon_m, codon_nrf = None, {}, None
    if cds_store is not None:
        codon_recs = backtranslate(read_fasta(prot_aln), tip_map, cds_store)
        if len(codon_recs) >= 4:
            codon_aln = os.path.join(workdir, f"{hog}.codon.aln")
            with open(codon_aln, "w") as out:
                for tip, seq in codon_recs:
                    out.write(f">{tip}\n{seq}\n")
            codon_tree, _ = build_tree_raxml(
                codon_aln, os.path.join(workdir, f"{hog}.codon"),
                raxml_path, codon_model, bs_trees)
            codon_m = tree_metrics(codon_tree)
            codon_nrf = species_concordance(codon_tree, species_tree) if species_tree else None

    pc_nrf = tree_agreement(prot_tree, codon_tree) if codon_tree else None
    status, reason = compute_status(prot_m, prot_nrf, codon_m, codon_nrf)

    row.update({
        "Prot_tips": prot_m.get("n_tips", ""),
        "Prot_support": prot_m.get("mean_support", ""),
        "Prot_nRF": "" if prot_nrf is None else prot_nrf,
        "Codon_support": codon_m.get("mean_support", ""),
        "Codon_patristic": codon_m.get("mean_patristic", ""),
        "Codon_nRF": "" if codon_nrf is None else codon_nrf,
        "Prot_vs_Codon_nRF": "" if pc_nrf is None else pc_nrf,
        "Status": status,
        "Reason": reason,
    })
    row["_tree"] = codon_tree or prot_tree   # representative tree (for plotting)
    return row


# ---------------------------------------------------------------------------
# Tree plotting
# ---------------------------------------------------------------------------

def plot_gene_tree(tree_path, out_png, title=""):
    """
    Render a gene tree (Newick, with bootstrap support) to PNG using Biopython +
    matplotlib (headless). Returns the path, or None if it could not be drawn.
    """
    if not tree_path or not os.path.exists(tree_path):
        return None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from Bio import Phylo
    except Exception:
        return None
    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception:
        return None
    n_tips = tree.count_terminals()
    fig = plt.figure(figsize=(7, max(2.5, 0.4 * n_tips)))
    ax = fig.add_subplot(1, 1, 1)
    try:
        Phylo.draw(tree, do_show=False, axes=ax,
                   show_confidence=True, branch_labels=None)
    except Exception:
        plt.close(fig)
        return None
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

_COLS = ["HOG", "Compartment", "N_species", "N_genes",
         "Prot_tips", "Prot_support", "Prot_nRF",
         "Codon_support", "Codon_patristic", "Codon_nRF", "Prot_vs_Codon_nRF",
         "Accession", "Gene", "BLAST_hit", "BLAST_pident", "BLAST_qcov",
         "Status", "Reason"]


def run(hogsfile, fasta_dir, outdir, prefix="", cds_dir=None, species_tree=None,
        per_class=3, hog_ids=None, mafft_path="mafft", raxml_path="raxml-ng",
        blastp_path="blastp", makeblastdb_path="makeblastdb",
        prot_model="LG+G", codon_model="GTR+G", bs_trees=100, writer=None):
    """
    Validate candidate HOGs and write a QC table.

    Non-private HOGs get a protein tree and (when ``cds_dir`` is given) a codon
    tree, compared and scored. Private HOGs are validated by BLAST against all
    other accessions. Trees go to ``<outdir>/gene_trees/`` and the table to
    ``<outdir>/<prefix>genetree_validation.tsv``.
    """
    species, occ = parse_hog_occupancy(hogsfile)
    n_species = len(species)
    prot_store = SeqStore(fasta_dir)
    cds_store = SeqStore(cds_dir) if cds_dir else None

    if hog_ids:
        targets = [h for h in hog_ids if h in occ]
    else:
        buckets = select_candidates(occ, n_species, per_class)
        targets = [h for c in ("core", "shell", "private") for h in buckets[c]]

    tree_dir = os.path.join(outdir, "gene_trees")
    os.makedirs(tree_dir, exist_ok=True)

    # BLAST DB (built once) if any private HOG is in the target set.
    db = None
    if any(compartment_of(occ[h], n_species) == "private" for h in targets):
        db = build_proteome_db(species, prot_store, tree_dir, makeblastdb_path)
        if db is None:
            print("[WARNING] BLAST+ (makeblastdb/blastp) not available — "
                  "private HOGs cannot be validated.")

    rows = []
    for hog in targets:
        comp = compartment_of(occ[hog], n_species)
        print(f"[INFO] gene-tree validation: {hog} ({comp})")
        if comp == "private":
            rows.append(validate_private(hog, occ[hog], species, prot_store, db,
                                         tree_dir, blastp_path))
        else:
            rows.append(validate_hog(
                hog, occ[hog], n_species, prot_store, cds_store, tree_dir,
                species_tree, mafft_path, raxml_path,
                prot_model, codon_model, bs_trees))

    # Render one representative gene tree per compartment so the clustering can be
    # eyeballed (core / shell; private has no tree).
    plots = {}
    for comp in ("core", "shell"):
        rep = next((r for r in rows
                    if r.get("Compartment") == comp and r.get("_tree")), None)
        if rep:
            png = os.path.join(tree_dir, f"{prefix}genetree_{comp}_example.png")
            if plot_gene_tree(rep["_tree"], png,
                              title=f"{comp.capitalize()} example: {rep['HOG']} "
                                    f"[{rep.get('Status', '')}]"):
                plots[comp] = png
                print(f"[INFO] Saved {comp} example gene-tree plot -> {png}")

    out_tsv = os.path.join(outdir, f"{prefix}genetree_validation.tsv")
    if writer is None:
        os.makedirs(outdir, exist_ok=True)
        with open(out_tsv, "w") as fh:
            fh.write("\t".join(_COLS) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "")) for c in _COLS) + "\n")
    else:
        writer(rows, out_tsv)

    # Flagged HOGs: gene tree conflicts with / can't confirm the ortholog grouping —
    # the ones to be careful about when interpreting the compartment.
    flagged = [r for r in rows if str(r.get("Status", "")).startswith(
        ("REVIEW", "Conflict", "Redundant"))]
    flagged_tsv = os.path.join(outdir, f"{prefix}genetree_flagged.tsv")
    if writer is None:
        with open(flagged_tsv, "w") as fh:
            fh.write("\t".join(_COLS) + "\n")
            for r in flagged:
                fh.write("\t".join(str(r.get(c, "")) for c in _COLS) + "\n")
    else:
        writer(flagged, flagged_tsv)

    print(f"[INFO] Saved gene-tree validation -> {out_tsv} "
          f"({len(rows)} HOGs, {len(flagged)} flagged -> {flagged_tsv})")
    return {"table": out_tsv, "flagged": flagged_tsv, "tree_dir": tree_dir,
            "rows": rows, "plots": plots}
