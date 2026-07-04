#!/usr/bin/env python3
"""
panhog_pansummary.py
====================
Pangenome-composition summary + plots for a core/shell/private classification.

Works for either the frequency-based classification (--pan) or the PD-weighted one
(--pan-weighted), so the two can be produced and compared with identical plots:
  * <prefix>pangenome_summary.tsv    compartment counts + percentages
  * <prefix>occupancy_histogram.tsv  families per occupancy level per compartment
  * <prefix>pangenome_pie.png/svg    pie chart
  * <prefix>occupancy_histogram.png/svg  U-shaped occupancy histogram (stacked)
  * <prefix>compartment_stacked_bar.png/svg  single stacked bar

Both plotted TSVs are written so users can re-draw the figures themselves.
"""

import os
from collections import Counter

# "Watermelon Sorbet" palette (matches the Ka/Ks compartment plot).
WATERMELON = {"core": "#ef476f", "shell": "#ffd166", "private": "#06d6a0"}
_ORDER = ["core", "shell", "private"]
_LABELS = {"core": "Core", "shell": "Shell", "private": "Private"}


def summarize(class_of, occupancy_of, outdir, prefix, approach="frequency",
              n_accessions=None, colors=None):
    """
    class_of     : {hog: 'core'|'shell'|'private'}
    occupancy_of : {hog: number_of_accessions_carrying_it}
    approach     : label used in titles/filenames context ('frequency' or 'PD-weighted')
    Returns a dict of written paths.
    """
    colors = colors or WATERMELON
    os.makedirs(outdir, exist_ok=True)
    comp_counts = Counter(class_of.values())
    total = sum(comp_counts.values()) or 1

    # --- summary tsv ---
    summ_tsv = os.path.join(outdir, f"{prefix}pangenome_summary.tsv")
    with open(summ_tsv, "w") as fh:
        fh.write("Compartment\tN_Families\tPercent\tApproach\n")
        for c in _ORDER:
            n = comp_counts.get(c, 0)
            fh.write(f"{c}\t{n}\t{100 * n / total:.2f}\t{approach}\n")
        fh.write(f"total\t{total}\t100.00\t{approach}\n")

    # --- occupancy histogram data ---
    if n_accessions is None:
        n_accessions = max(occupancy_of.values()) if occupancy_of else 0
    hist = {k: Counter() for k in range(1, n_accessions + 1)}
    for hog, k in occupancy_of.items():
        if 1 <= k <= n_accessions:
            hist[k][class_of.get(hog, "shell")] += 1
    hist_tsv = os.path.join(outdir, f"{prefix}occupancy_histogram.tsv")
    with open(hist_tsv, "w") as fh:
        fh.write("Num_Accessions\tCore\tShell\tPrivate\tTotal\n")
        for k in range(1, n_accessions + 1):
            c = hist[k]
            fh.write(f"{k}\t{c.get('core', 0)}\t{c.get('shell', 0)}\t"
                     f"{c.get('private', 0)}\t{sum(c.values())}\n")

    out = {"summary": summ_tsv, "histogram_tsv": hist_tsv}

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover
        print("[WARNING] matplotlib unavailable; wrote pangenome summary TSVs only.")
        return out

    def _save(fig, name):
        for ext in ("png", "svg"):
            fig.savefig(os.path.join(outdir, f"{prefix}{name}.{ext}"),
                        dpi=150, bbox_inches="tight")
        plt.close(fig)

    # --- pie ---
    fig, ax = plt.subplots(figsize=(4.6, 4.6))
    vals = [comp_counts.get(c, 0) for c in _ORDER]
    wedges, _, _ = ax.pie(
        vals, colors=[colors[c] for c in _ORDER],
        autopct=lambda p: f"{p:.1f}%", startangle=90,
        wedgeprops=dict(edgecolor="white", linewidth=1.5),
        textprops=dict(color="black", fontsize=10, fontweight="bold"))
    ax.legend(wedges,
              [f"{_LABELS[c]}: {comp_counts.get(c, 0):,} "
               f"({100 * comp_counts.get(c, 0) / total:.1f}%)" for c in _ORDER],
              loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=1,
              fontsize=9, frameon=False)
    ax.set_title(f"Pangenome composition ({approach})")
    _save(fig, "pangenome_pie")

    # --- U-shaped occupancy histogram (stacked by compartment) ---
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ks = list(range(1, n_accessions + 1))
    bottoms = [0] * len(ks)
    for c in _ORDER:
        h = [hist[k].get(c, 0) for k in ks]
        ax.bar(ks, h, bottom=bottoms, color=colors[c], label=_LABELS[c], width=0.85)
        bottoms = [b + x for b, x in zip(bottoms, h)]
    for k in ks:
        tot = sum(hist[k].values())
        if tot and (k == 1 or k == n_accessions):
            ax.text(k, tot, f"{tot:,}", ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("Number of accessions")
    ax.set_ylabel("Number of gene families")
    ax.set_title(f"Gene-family occupancy ({approach})")
    ax.legend(frameon=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    _save(fig, "occupancy_histogram")

    # --- single stacked bar ---
    fig, ax = plt.subplots(figsize=(3.4, 4.6))
    b = 0
    for c in _ORDER:
        n = comp_counts.get(c, 0)
        ax.bar(0, n, bottom=b, color=colors[c], width=0.6, edgecolor="white")
        if n:
            ax.text(0, b + n / 2, f"{_LABELS[c]}\n{n:,}", ha="center",
                    va="center", fontsize=9, fontweight="bold")
        b += n
    ax.set_xticks([])
    ax.set_ylabel("Number of gene families")
    ax.set_title(f"Compartments ({approach})")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    _save(fig, "compartment_stacked_bar")

    print(f"[INFO] Saved pangenome summary ({approach}) -> {summ_tsv} "
          f"+ pie / occupancy histogram / stacked bar (png+svg) + {hist_tsv}")
    return out
