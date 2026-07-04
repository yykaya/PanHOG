#!/usr/bin/env python3
"""
panhog_kaks_compartments.py
===========================
Compare codon-based dN/dS (Ka/Ks) across pangenome compartments.

Reproduces the classic pangenome selection pattern — core genes under the
strongest purifying selection (lowest Ka/Ks), shell / dispensable genes more
relaxed — by sampling HOGs from each compartment, computing one dN/dS per HOG,
and drawing a box plot of the distributions with a Kruskal-Wallis test.

Per-HOG dN/dS is the **ratio of means** (mean(dN)/mean(dS)); the mean of the
per-pair ratios is unstable (a near-zero dS sends a pair's ratio to infinity)
and is *not* used — see ``panhog_dnds.summarize_hog``.

Notes on compartments:
  * core / shell  — pairwise dN/dS among the orthologous copies (needs ≥2
    accessions present).
  * private       — a single-copy private gene has no pair and is skipped;
    multi-copy private HOGs give *paralog* dN/dS, which is a different (younger)
    comparison than the ortholog dN/dS of core/shell. Interpret accordingly.

Reuses collision-safe per-accession sequence loading + codon back-translation
from ``panhog_genetrees`` and the dN/dS engine from ``panhog_dnds``.
"""

from __future__ import annotations

import math
import os
import statistics

import panhog_dnds
import panhog_genetrees as gt

try:
    from scipy.stats import kruskal, mannwhitneyu
    HAS_SCIPY = True
except Exception:  # pragma: no cover
    HAS_SCIPY = False


def sample_by_compartment(occ, n_species, per_compartment=100,
                          min_seqs=2, max_seqs=12):
    """
    Pick HOGs per compartment with ≥ ``min_seqs`` sequences (needed for a
    pairwise dN/dS) and ≤ ``max_seqs`` (to keep alignments fast). Returns
    ``{compartment: [hog_ids]}``.
    """
    buckets = {"core": [], "shell": [], "private": []}
    for hog, genes in occ.items():
        nseq = sum(len(v) for v in genes.values())
        if nseq < min_seqs or nseq > max_seqs:
            continue
        c = gt.compartment_of(genes, n_species)
        if len(buckets[c]) < per_compartment:
            buckets[c].append(hog)
    return buckets


def hog_dnds(hog, occ_hog, prot_store, cds_store, workdir,
             method="biopython", model="NG86", mafft_path="mafft", reference=None):
    """
    Codon-align a HOG and return its robust per-HOG dN/dS (ratio of means), or None.

    With ``reference`` (an accession name) only HOGs that contain the reference are
    used, and dN/dS is computed **reference-vs-rest** — one value per reference gene
    against its orthologs, matching a "<reference> Ka/Ks" analysis. (A private gene
    unique to the reference then has no ortholog to compare against and is dropped,
    which is the correct behaviour.)
    """
    if reference and reference not in occ_hog:
        return None
    prot_fa = os.path.join(workdir, f"{hog}.faa")
    tip_map = gt.write_protein_fasta(occ_hog, prot_store, prot_fa)
    if len(tip_map) < 2:
        return None
    aln = os.path.join(workdir, f"{hog}.aln")
    if not gt.align_mafft(prot_fa, aln, mafft_path):
        return None
    codon_recs = gt.backtranslate(gt.read_fasta(aln), tip_map, cds_store)
    if len(codon_recs) < 2:
        return None
    species_of = {tip: acc for tip, (acc, _g) in tip_map.items()}
    rows = panhog_dnds.dnds_from_alignment(codon_recs, method=method, model=model,
                                           species_of=species_of, reference=reference)
    v = panhog_dnds.summarize_hog(rows)["dN_dS"]
    return v if isinstance(v, float) and not math.isnan(v) else None


def compartment_stats(data):
    """Medians per compartment + Kruskal-Wallis P across compartments."""
    out = {}
    for c, v in data.items():
        if v:
            out[f"{c}_median"] = round(statistics.median(v), 4)
            out[f"{c}_n"] = len(v)
    groups = [v for v in (data.get("core"), data.get("shell"), data.get("private"))
              if v and len(v) >= 2]
    if HAS_SCIPY and len(groups) >= 2:
        out["kruskal_p"] = float(kruskal(*groups).pvalue)
    return out


def plot_compartment_dnds(data, outpath, pvalue=None, ymax=None):
    """
    Box plot of per-HOG dN/dS by compartment (core red, shell blue, private
    green), styled after published pangenome Ka/Ks figures. Writes PNG/PDF/SVG.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"core": "#d62728", "shell": "#1f77b4", "private": "#2ca02c"}
    labels = {"core": "Core", "shell": "Shell", "private": "Private"}
    order = [c for c in ("core", "shell", "private") if data.get(c)]
    vals = [data[c] for c in order]

    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    bp = ax.boxplot(vals, patch_artist=True, widths=0.6, showfliers=False,
                    medianprops=dict(color="black", linewidth=1.4),
                    whiskerprops=dict(color="black"), capprops=dict(color="black"))
    for patch, c in zip(bp["boxes"], order):
        patch.set_facecolor(colors[c])
        patch.set_edgecolor("black")
        patch.set_alpha(0.85)
    ax.set_xticklabels([f"{labels[c]}\n(n={len(data[c])})" for c in order])
    ax.set_ylabel(r"$K_a/K_s$", fontsize=12)
    ax.set_ylim(0, ymax if ymax else None)
    ax.axhline(1.0, ls="--", lw=0.8, color="grey", alpha=0.7)  # neutral expectation
    if pvalue is not None:
        txt = ("P < 2.2e-16" if pvalue < 2.2e-16
               else f"Kruskal–Wallis\nP = {pvalue:.2g}")
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top", fontsize=9)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    stem = outpath.rsplit(".", 1)[0]
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{stem}.{ext}", dpi=150)
    plt.close(fig)
    return f"{stem}.png"


def run(hogsfile, fasta_dir, cds_dir, outdir, prefix="", per_compartment=100,
        max_seqs=12, method="biopython", model="NG86", mafft_path="mafft",
        reference=None, make_plot=True):
    """
    Sample HOGs per compartment, compute per-HOG dN/dS, write a per-HOG table and
    a compartment box plot. Requires ``cds_dir`` (codon alignments).

    With ``reference`` (an accession) dN/dS is reference-vs-rest — each retained
    HOG contributes the selection on the reference gene against its orthologs
    (a "<reference> Ka/Ks" analysis); otherwise all unordered pairs are used.
    """
    if not cds_dir:
        print("[ERROR] --cds is required for the compartment Ka/Ks comparison. Skipping.")
        return {}
    species, occ = gt.parse_hog_occupancy(hogsfile)
    n_species = len(species)
    prot_store = gt.SeqStore(fasta_dir)
    cds_store = gt.SeqStore(cds_dir)
    buckets = sample_by_compartment(occ, n_species, per_compartment, 2, max_seqs)

    workdir = os.path.join(outdir, "kaks_compartments_tmp")
    os.makedirs(workdir, exist_ok=True)

    data = {"core": [], "shell": [], "private": []}
    rows_out = []
    for comp in ("core", "shell", "private"):
        for hog in buckets[comp]:
            v = hog_dnds(hog, occ[hog], prot_store, cds_store, workdir,
                         method, model, mafft_path, reference=reference)
            if v is not None:
                data[comp].append(v)
                rows_out.append((hog, comp, v))
        print(f"[INFO] {comp}: {len(data[comp])} HOGs with usable dN/dS")

    os.makedirs(outdir, exist_ok=True)
    tsv = os.path.join(outdir, f"{prefix}kaks_by_compartment.tsv")
    with open(tsv, "w") as fh:
        fh.write("HOG\tCompartment\tdN_dS\n")
        for hog, comp, v in rows_out:
            fh.write(f"{hog}\t{comp}\t{v:.4f}\n")

    stats = compartment_stats(data)
    plot_path = None
    if make_plot and any(data.values()):
        plot_path = plot_compartment_dnds(
            data, os.path.join(outdir, f"{prefix}kaks_by_compartment.png"),
            pvalue=stats.get("kruskal_p"))

    print(f"[INFO] Saved compartment Ka/Ks table -> {tsv}")
    for c in ("core", "shell", "private"):
        if data[c]:
            print(f"[INFO]   {c}: median dN/dS = {stats.get(f'{c}_median')} "
                  f"(n={len(data[c])})")
    if "kruskal_p" in stats:
        print(f"[INFO]   Kruskal-Wallis P = {stats['kruskal_p']:.3g}")
    if plot_path:
        print(f"[INFO] Saved compartment Ka/Ks plot  -> {plot_path}")
    return {"table": tsv, "plot": plot_path, "data": data, "stats": stats}
