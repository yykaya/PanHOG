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
# Tier 2.1: phylogenetic-diversity-weighted core / shell / private
# ---------------------------------------------------------------------------

def pd_weighted_classification(tree, presence, core_threshold=0.9,
                               private_threshold=0.1, below=None, total_bl=None):
    """
    Classify each HOG by the fraction of total tree branch length its carrier
    species span (Faith's PD fraction) instead of by a flat species count.

    A HOG whose carriers span >= ``core_threshold`` of the tree is ``core``, one
    that spans <= ``private_threshold`` is ``private``, and anything between is
    ``shell``. This up-weights HOGs spread across deep, divergent lineages and
    down-weights those confined to a few closely related tips, so a HOG in two
    distant clades is not treated the same as one in two sister tips.

    Returns (rows, counts):
      rows   : list of dicts (HOG, Num_Species, PD_Fraction, Weighted_Class)
      counts : {'core': n, 'shell': n, 'private': n}
    """
    if below is None:
        below = tips_below(tree)
    if total_bl is None:
        total_bl = total_branch_length(tree)
    tip_by_name = {t.name: t for t in tree.get_terminals()}

    rows = []
    counts = {"core": 0, "shell": 0, "private": 0}
    for hog, present in presence.items():
        present = [s for s in present if s in tip_by_name]
        if not present:
            continue
        frac = pd_fraction(tree, present, below, total_bl)
        if frac >= core_threshold:
            cls = "core"
        elif frac <= private_threshold:
            cls = "private"
        else:
            cls = "shell"
        counts[cls] += 1
        rows.append({
            "HOG": hog,
            "Num_Species": len(present),
            "PD_Fraction": round(frac, 4),
            "Weighted_Class": cls,
        })
    return rows, counts


# ---------------------------------------------------------------------------
# Tier 2.2: clade-conditioned compartments
# ---------------------------------------------------------------------------

def clade_compartments(tree, presence, below=None, min_clade_size=2):
    """
    For every internal node (clade), classify each HOG's occupancy *within that
    clade's own tip set* and tally clade-core / clade-shell / clade-private.

    For a clade with ``k`` tips, a HOG carried by ``c`` of those tips is:
      * clade-core    if c == k   (present across the whole clade)
      * clade-private if c == 1   (a single tip within the clade)
      * clade-shell   if 1 < c < k
    HOGs absent from the clade (c == 0) are not counted for that clade.

    Returns rows: list of dicts
      (Node, Num_Tips, Clade_Core, Clade_Shell, Clade_Private, HOGs_Present)
    one per internal node with at least ``min_clade_size`` tips.
    """
    if below is None:
        below = tips_below(tree)
    carriers = {hog: set(sps) for hog, sps in presence.items()}

    rows = []
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            continue
        clade_tips = below[id(clade)]
        k = len(clade_tips)
        if k < min_clade_size:
            continue
        core = shell = private = present_total = 0
        for present in carriers.values():
            c = len(clade_tips & present)
            if c == 0:
                continue
            present_total += 1
            if c == k:
                core += 1
            elif c == 1:
                private += 1
            else:
                shell += 1
        rows.append({
            "Node": clade.name if clade.name else "root",
            "Num_Tips": k,
            "Clade_Core": core,
            "Clade_Shell": shell,
            "Clade_Private": private,
            "HOGs_Present": present_total,
        })
    return rows


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


def analyze(dGeneNumbers, dSpecies, species_tree_file, outdir, prefix, writer=None,
            pan_weighted=False, pan_weighted_core=0.9, pan_weighted_private=0.1):
    """
    Run the full phylogeny-aware analysis and write output tables + an annotated
    tree. ``writer`` is an optional callable(df, path) for TSV output; when None
    a pandas-based writer is used. Returns a dict of output file paths.

    When ``pan_weighted`` is True, also writes a phylogenetic-diversity-weighted
    core/shell/private classification (thresholds ``pan_weighted_core`` /
    ``pan_weighted_private``) and a per-clade compartment table.
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

    if pan_weighted:
        w_rows, w_counts = pd_weighted_classification(
            tree, presence, pan_weighted_core, pan_weighted_private, below)
        c_rows = clade_compartments(tree, presence, below)

        # Confidence: does the frequency-based compartment (core = all tips,
        # private = 1 tip, shell = in between) agree with the PD-weighted class?
        n_tips = len(tree.get_terminals())

        def _flat(k):
            return "core" if k >= n_tips else ("private" if k == 1 else "shell")

        conf_rows, n_same, transitions = [], 0, Counter()
        for w in w_rows:
            fc, wc = _flat(w["Num_Species"]), w["Weighted_Class"]
            agree = fc == wc
            n_same += agree
            if not agree:
                transitions[f"{fc}->{wc}"] += 1
            conf_rows.append({
                "HOG": w["HOG"], "Num_Species": w["Num_Species"],
                "Flat_Class": fc, "PD_Fraction": w["PD_Fraction"],
                "Weighted_Class": wc,
                "Agreement": "Same" if agree else "Reclassified",
                "Change": "" if agree else f"{fc}->{wc}",
            })

        paths["weighted"] = os.path.join(outdir, f"{prefix}hog_pd_weighted_class.tsv")
        paths["clade"] = os.path.join(outdir, f"{prefix}clade_compartments.tsv")
        paths["confidence"] = os.path.join(outdir, f"{prefix}classification_confidence.tsv")
        paths["reclassified"] = os.path.join(outdir, f"{prefix}reclassified_HOGs.tsv")
        writer(w_rows, paths["weighted"])
        writer(c_rows, paths["clade"])
        writer(conf_rows, paths["confidence"])
        # Just the HOGs whose compartment changed under PD-weighting (the ones to
        # look at when the count and the phylogeny disagree).
        writer([r for r in conf_rows if r["Agreement"] == "Reclassified"],
               paths["reclassified"])
        print(f"[INFO] Saved PD-weighted class-> {paths['weighted']} "
              f"(core={w_counts['core']}, shell={w_counts['shell']}, "
              f"private={w_counts['private']})")
        print(f"[INFO] Saved clade compartments-> {paths['clade']} "
              f"({len(c_rows)} clades)")
        pct = 100.0 * n_same / len(conf_rows) if conf_rows else 0.0
        print(f"[INFO] Saved classification confidence-> {paths['confidence']} "
              f"({n_same}/{len(conf_rows)} = {pct:.1f}% agree with frequency; "
              f"{len(conf_rows) - n_same} reclassified; "
              f"top: {dict(transitions.most_common(3))})")

    _write_phylo_readme(outdir, prefix, pan_weighted)
    return paths


def _write_phylo_readme(outdir, prefix, pan_weighted):
    """Write a README explaining every file in the phylogeny_weighted directory."""
    p = prefix
    lines = f"""# phylogeny_weighted/ — how to read these files

Phylogeny-aware layer on top of the frequency-based --pan classification. It asks
WHERE on the species tree each HOG sits and HOW MUCH of the tree it spans.

## {p}hog_lca_analysis.tsv
One row per HOG mapped onto the species tree.
  HOG           orthogroup id
  Num_Species   number of carrier accessions
  LCA_Node      lowest common ancestor node (postorder-named N0..Nk) of all carriers
  PD_Fraction   Faith's phylogenetic diversity the carriers span, as a FRACTION of
                the tree's total branch length. Range 0..1 (see below).
  Species_List  the carrier accessions

## PD_Fraction — what the number means (0 -> 1)
It is NOT the compartment; it is the share of the tree's total branch length that the
carrier tips collectively cover:
  * 1.0  carriers reach every corner of the tree  (e.g. a true core gene)
  * ~0   carriers are one tip, or tips sitting on ZERO-length branches
The PD-weighted class is then thresholded on it:
  PD_Fraction >= 0.9 (--pan-weighted-core)     -> core     (spans ~the whole tree)
  PD_Fraction <= 0.1 (--pan-weighted-private)  -> private  (spans ~nothing)
  in between                                   -> shell
Two HOGs with the SAME Num_Species can get DIFFERENT PD_Fraction: 2 deeply divergent
accessions span more branch length than 2 sister accessions. That is the whole point
of weighting by phylogeny instead of by a raw count.
Edge case: a HOG in 2 accessions whose tips are on zero-length branches has
PD_Fraction = 0.0 and is called 'private' even though Num_Species = 2 (it adds no
phylogenetic breadth).

## {p}hog_pd_weighted_class.tsv
  HOG, Num_Species, PD_Fraction, Weighted_Class (core/shell/private by the thresholds).

## {p}classification_confidence.tsv
Frequency-based vs PD-weighted compartment, per HOG:
  Flat_Class     core = all accessions, private = 1, else shell
  Weighted_Class the PD-weighted class
  Agreement      Same | Reclassified
  Change         e.g. shell->core when they disagree
'Same' rows are HIGH-CONFIDENCE compartment calls; 'Reclassified' rows are where the
count and the phylogeny disagree -- inspect those.

## {p}reclassified_HOGs.tsv
Just the 'Reclassified' rows above (the ones to look at).

## {p}clade_compartments.tsv
For every internal node, core/shell/private tallied WITHIN that clade's own tips
(Num_Tips, Clade_Core/Shell/Private, HOGs_Present).

## {p}phylo_node_summary.tsv
Dollo gain/loss per branch: HOGs_LCA (origins), Gains, Losses per node.

## {p}species_tree_annotated.nwk
Species tree with internal-node names (N0..Nk) matching LCA_Node / node summary.

## {p}genetree_validation.tsv + gene_trees/
Per-HOG gene-tree validation (protein + codon ML trees) with a Status
(Confirmed / Redundant / Conflict / Low-signal). 'Conflict' = the HOG's own gene tree
is well supported but disagrees with the species tree -> members are probably NOT
clean orthologs; be careful using this HOG. gene_trees/ holds the trees + example
plots. See {p}genetree_flagged.tsv for just the flagged HOGs.
"""
    try:
        with open(os.path.join(outdir, "README.md"), "w") as fh:
            fh.write(lines)
    except OSError:
        pass
