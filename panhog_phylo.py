#!/usr/bin/env python3
"""
panhog_phylo.py
===============
Phylogeny-aware analysis of HOGs on a species tree.

This replaces the earlier single-pass LCA routine (which named internal nodes
with ``Internal_Node_<id(...)>`` — a Python memory address, so results were not
reproducible and per-node counts collapsed together) and adds the features that
make the "phylogeny-aware" label meaningful:

* Deterministic internal-node naming (postorder ``N0..Nk``) and an annotated
  Newick tree, so every node in the output tables is identifiable and can be
  visualised in iTOL / FigTree.
* Species-set validation against the tree tips (with light name normalisation).
* HOG -> LCA mapping with a phylogenetic-diversity (Faith's PD) fraction.
* Dollo-parsimony gain/loss reconstruction: each HOG is gained once at the LCA
  of the species that carry it and lost on the branches leading to maximal
  clades that lack it, giving per-branch gene-family gain/loss counts.

Only Biopython (already a PanHOG dependency) is used — no ete3/dendropy needed.
"""

from __future__ import annotations

import os
from collections import Counter

try:
    from Bio import Phylo
    HAS_BIOPYTHON = True
except ImportError:  # pragma: no cover
    HAS_BIOPYTHON = False


# ---------------------------------------------------------------------------
# Tree loading / preparation
# ---------------------------------------------------------------------------

def load_species_tree(path):
    """Read a Newick species tree. Raises on failure (caller handles)."""
    if not HAS_BIOPYTHON:
        raise RuntimeError("Biopython is required for phylogenetic analysis.")
    return Phylo.read(path, "newick")


def tip_names(tree):
    return [t.name for t in tree.get_terminals()]


def name_internal_nodes(tree, prefix="N"):
    """
    Assign deterministic names to unnamed internal nodes in postorder
    (``N0``, ``N1``, ...). Existing internal names are preserved. Returns the
    tree for chaining.
    """
    i = 0
    for clade in tree.find_clades(order="postorder"):
        if not clade.is_terminal() and not clade.name:
            clade.name = f"{prefix}{i}"
            i += 1
    return tree


def tips_below(tree):
    """Map id(clade) -> frozenset of descendant tip names (one postorder pass)."""
    below = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            below[id(clade)] = frozenset([clade.name])
        else:
            acc = set()
            for child in clade.clades:
                acc |= below[id(child)]
            below[id(clade)] = frozenset(acc)
    return below


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _normalize(name):
    """Lowercase and strip common FASTA suffixes for fuzzy tip matching."""
    n = name.strip()
    for suf in (".pep.fa", ".pep", ".faa", ".fasta", ".fa", ".fna", ".cds"):
        if n.lower().endswith(suf):
            n = n[: -len(suf)]
            break
    return n.lower()


def validate_tree_species(tree, species):
    """
    Compare tree tips with the analysis species set.

    Returns (name_map, missing_in_tree, extra_in_tree) where name_map maps
    analysis-species -> tree-tip name (via exact match first, then a normalised
    match). ``missing_in_tree`` are species with no tip; ``extra_in_tree`` are
    tips with no matching species.
    """
    tips = tip_names(tree)
    tip_set = set(tips)
    norm_tip = {}
    for t in tips:
        norm_tip.setdefault(_normalize(t), t)

    name_map = {}
    missing = []
    for sp in species:
        if sp in tip_set:
            name_map[sp] = sp
        elif _normalize(sp) in norm_tip:
            name_map[sp] = norm_tip[_normalize(sp)]
        else:
            missing.append(sp)

    mapped_tips = set(name_map.values())
    extra = [t for t in tips if t not in mapped_tips]
    return name_map, missing, extra


# ---------------------------------------------------------------------------
# LCA mapping + phylogenetic diversity
# ---------------------------------------------------------------------------

def total_branch_length(tree):
    return sum(c.branch_length for c in tree.find_clades() if c.branch_length)


def _lca(tree, tip_clades):
    if len(tip_clades) == 1:
        return tip_clades[0]
    return tree.common_ancestor(tip_clades)


def pd_fraction(tree, present_tips, below=None, total_bl=None):
    """
    Faith's PD of ``present_tips`` as a fraction of the tree's total branch
    length: the summed length of the minimal subtree spanning those tips over
    the summed length of all branches. Falls back to tip fraction when the tree
    has no branch lengths (a cladogram).
    """
    tips = tip_names(tree)
    n_all = len(tips)
    present = set(present_tips)
    if not present:
        return 0.0
    if total_bl is None:
        total_bl = total_branch_length(tree)
    if not total_bl:  # cladogram: no lengths
        return len(present) / n_all if n_all else 0.0
    if below is None:
        below = tips_below(tree)

    tip_by_name = {t.name: t for t in tree.get_terminals()}
    present_clades = [tip_by_name[s] for s in present if s in tip_by_name]
    if not present_clades:
        return 0.0
    lca = _lca(tree, present_clades)
    # Sum branch lengths of clades in the LCA subtree that lead to a present tip.
    pd = 0.0
    for clade in lca.find_clades():
        if clade is lca:
            continue
        if below[id(clade)] & present and clade.branch_length:
            pd += clade.branch_length
    return pd / total_bl if total_bl else 0.0


def map_hogs_to_lca(tree, presence, below=None):
    """
    For each HOG, find the LCA node of the species that carry it.

    ``presence`` is {hog_id: iterable_of_tip_names}. Returns
    (rows, per_node_counter):
      rows: list of dicts (HOG, Num_Species, LCA_Node, PD_Fraction, Species_List)
      per_node_counter: Counter of HOGs whose LCA is each node.
    """
    if below is None:
        below = tips_below(tree)
    total_bl = total_branch_length(tree)
    tip_by_name = {t.name: t for t in tree.get_terminals()}

    rows = []
    per_node = Counter()
    for hog, present in presence.items():
        present = [s for s in present if s in tip_by_name]
        if not present:
            continue
        lca = _lca(tree, [tip_by_name[s] for s in present])
        node_name = lca.name if lca.name else "root"
        per_node[node_name] += 1
        rows.append({
            "HOG": hog,
            "Num_Species": len(present),
            "LCA_Node": node_name,
            "PD_Fraction": round(pd_fraction(tree, present, below, total_bl), 4),
            "Species_List": ",".join(sorted(present)),
        })
    return rows, per_node


# ---------------------------------------------------------------------------
# Dollo parsimony gain / loss
# ---------------------------------------------------------------------------

def dollo_events(tree, present_tips, below=None):
    """
    Dollo reconstruction for one HOG (gained once, lost many times).

    Returns (gain_node_name, [loss_node_names]). The gene is placed at the LCA
    of the carriers; losses are the maximal clades under that LCA that contain
    no carrier.
    """
    if below is None:
        below = tips_below(tree)
    tip_by_name = {t.name: t for t in tree.get_terminals()}
    present = set(s for s in present_tips if s in tip_by_name)
    if not present:
        return None, []

    gain = _lca(tree, [tip_by_name[s] for s in present])
    losses = []

    def recurse(clade):
        # clade is known to contain >=1 carrier; descend to find lost subclades.
        if clade.is_terminal():
            return
        for child in clade.clades:
            if below[id(child)] & present:
                recurse(child)
            else:
                losses.append(child.name if child.name else "unnamed")

    recurse(gain)
    gain_name = gain.name if gain.name else "root"
    return gain_name, losses


def gain_loss_table(tree, presence, below=None):
    """
    Aggregate Dollo events across all HOGs.

    Returns (hog_rows, node_gain_counter, node_loss_counter):
      hog_rows: list of dicts (HOG, Gain_Node, Num_Losses, Loss_Nodes)
      node_gain_counter / node_loss_counter: per-node event totals.
    """
    if below is None:
        below = tips_below(tree)
    hog_rows = []
    gains = Counter()
    losses = Counter()
    for hog, present in presence.items():
        gain_node, loss_nodes = dollo_events(tree, present, below)
        if gain_node is None:
            continue
        gains[gain_node] += 1
        for ln in loss_nodes:
            losses[ln] += 1
        hog_rows.append({
            "HOG": hog,
            "Gain_Node": gain_node,
            "Num_Losses": len(loss_nodes),
            "Loss_Nodes": ",".join(loss_nodes),
        })
    return hog_rows, gains, losses


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def presence_from_gene_numbers(dGeneNumbers, dSpecies):
    """
    Build {hog_id: set(species_present)} from PanHOG's dGeneNumbers (per-HOG
    count vector in sorted-species-index order) and dSpecies (index -> name).
    """
    order = [dSpecies[i] for i in sorted(dSpecies.keys())]
    presence = {}
    for hog, counts in dGeneNumbers.items():
        present = {order[i] for i, c in enumerate(counts) if i < len(order) and c > 0}
        if present:
            presence[hog] = present
    return presence


def analyze(dGeneNumbers, dSpecies, species_tree_file, outdir, prefix, writer=None):
    """
    Run the full phylogeny-aware analysis and write output tables + an annotated
    tree. ``writer`` is an optional callable(df, path) for TSV output; when None
    a pandas-based writer is used. Returns a dict of output file paths.
    """
    if not HAS_BIOPYTHON:
        print("[ERROR] Biopython is required for phylogenetic analysis. Skipping.")
        return {}

    try:
        tree = load_species_tree(species_tree_file)
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Failed to read species tree {species_tree_file}: {e}")
        return {}

    name_internal_nodes(tree)
    below = tips_below(tree)

    species = [dSpecies[i] for i in sorted(dSpecies.keys())]
    name_map, missing, extra = validate_tree_species(tree, species)
    if missing:
        print(f"[WARNING] {len(missing)} analysis species not found as tree tips: "
              f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}")
    if extra:
        print(f"[WARNING] {len(extra)} tree tips have no matching species: "
              f"{', '.join(extra[:8])}{' ...' if len(extra) > 8 else ''}")

    # Map analysis species names onto the actual tip names before building presence.
    presence = presence_from_gene_numbers(dGeneNumbers, dSpecies)
    presence = {
        hog: {name_map.get(sp, sp) for sp in sps}
        for hog, sps in presence.items()
    }

    print(f"\n[INFO] Phylogeny-aware analysis on {len(species)} species, "
          f"{len(presence)} HOGs...")

    lca_rows, per_node_lca = map_hogs_to_lca(tree, presence, below)
    gl_rows, gains, losses = gain_loss_table(tree, presence, below)

    # Per-node summary combining LCA origins, gains and losses.
    all_nodes = set(per_node_lca) | set(gains) | set(losses)
    node_rows = [{
        "Node": n,
        "HOGs_LCA": per_node_lca.get(n, 0),
        "Gains": gains.get(n, 0),
        "Losses": losses.get(n, 0),
    } for n in sorted(all_nodes)]

    if writer is None:
        import pandas as pd

        def writer(rows, path):
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    os.makedirs(outdir, exist_ok=True)
    paths = {
        "lca": os.path.join(outdir, f"{prefix}hog_lca_analysis.tsv"),
        "gainloss": os.path.join(outdir, f"{prefix}hog_gainloss.tsv"),
        "nodes": os.path.join(outdir, f"{prefix}phylo_node_summary.tsv"),
        "tree": os.path.join(outdir, f"{prefix}species_tree_annotated.nwk"),
    }
    writer(lca_rows, paths["lca"])
    writer(gl_rows, paths["gainloss"])
    writer(node_rows, paths["nodes"])
    Phylo.write(tree, paths["tree"], "newick")

    total_gains = sum(gains.values())
    total_losses = sum(losses.values())
    print(f"[INFO] Saved LCA table       -> {paths['lca']}")
    print(f"[INFO] Saved gain/loss table -> {paths['gainloss']}")
    print(f"[INFO] Saved per-node summary-> {paths['nodes']}")
    print(f"[INFO] Saved annotated tree  -> {paths['tree']}")
    print(f"[INFO] Dollo reconstruction: {total_gains} gains, {total_losses} losses "
          f"across {len(tree.get_nonterminals())} internal nodes.")
    return paths
