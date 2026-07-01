"""
Tests for panhog_phylo — LCA mapping, Dollo gain/loss, PD, validation.

Tree used throughout (branch lengths in parentheses):

        ((A:1,B:1):1, ((C:1,D:1):1, E:2):1);

    root
    |-- (A,B)        internal
    |     |-- A
    |     `-- B
    `-- (CD,E)       internal
          |-- (C,D)  internal
          |     |-- C
          |     `-- D
          `-- E

Run with:  pytest -q tests/test_phylo.py
"""

import os
import sys
from io import StringIO

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ph = pytest.importorskip("panhog_phylo")
if not ph.HAS_BIOPYTHON:
    pytest.skip("Biopython not installed", allow_module_level=True)

from Bio import Phylo

NWK = "((A:1.0,B:1.0):1.0,((C:1.0,D:1.0):1.0,E:2.0):1.0);"


def make_tree():
    tree = Phylo.read(StringIO(NWK), "newick")
    ph.name_internal_nodes(tree)
    return tree


def test_internal_nodes_named_deterministically():
    t1 = make_tree()
    t2 = make_tree()
    n1 = sorted(c.name for c in t1.get_nonterminals())
    n2 = sorted(c.name for c in t2.get_nonterminals())
    assert n1 == n2                       # reproducible
    assert all(name for name in n1)       # no None / empty names
    assert len(n1) == 4                   # root + (A,B) + (CD,E) + (C,D)


def test_lca_mapping_and_pd_fraction():
    tree = make_tree()
    below = ph.tips_below(tree)
    presence = {
        "HOG_all": {"A", "B", "C", "D", "E"},
        "HOG_cd": {"C", "D"},
    }
    rows, per_node = ph.map_hogs_to_lca(tree, presence, below)
    by_hog = {r["HOG"]: r for r in rows}

    # Whole-tree HOG spans all branch length -> PD fraction 1.0
    assert by_hog["HOG_all"]["PD_Fraction"] == pytest.approx(1.0)
    # {C,D}: subtree length 2 (C:1 + D:1) over total 9
    assert by_hog["HOG_cd"]["PD_Fraction"] == pytest.approx(2.0 / 9.0, abs=1e-4)
    assert by_hog["HOG_cd"]["Num_Species"] == 2


def test_dollo_single_gain_no_loss_for_clade():
    tree = make_tree()
    below = ph.tips_below(tree)
    gain, losses = ph.dollo_events(tree, {"C", "D"}, below)
    # gained at the (C,D) ancestor, both descendants present -> no losses
    assert losses == []
    assert gain  # a named internal node


def test_dollo_losses_are_maximal_absent_clades():
    tree = make_tree()
    below = ph.tips_below(tree)
    # Present in A and C only -> gained at root; lost in B, D and E lineages.
    gain, losses = ph.dollo_events(tree, {"A", "C"}, below)
    assert set(losses) == {"B", "D", "E"}
    assert len(losses) == 3


def test_dollo_single_tip_is_terminal_gain():
    tree = make_tree()
    below = ph.tips_below(tree)
    gain, losses = ph.dollo_events(tree, {"A"}, below)
    assert gain == "A"
    assert losses == []


def test_gain_loss_table_aggregates_per_node():
    tree = make_tree()
    below = ph.tips_below(tree)
    presence = {
        "h1": {"A", "C"},        # gain root; losses B, D, E
        "h2": {"C", "D"},        # gain (C,D); no loss
    }
    hog_rows, gains, losses = ph.gain_loss_table(tree, presence, below)
    assert len(hog_rows) == 2
    # Two gains total, three losses total (from h1)
    assert sum(gains.values()) == 2
    assert sum(losses.values()) == 3
    assert losses["B"] == 1 and losses["D"] == 1 and losses["E"] == 1


def test_validate_species_with_normalization():
    tree = make_tree()
    # "A.fa"/"b.faa" should map onto tips A / B; "Z" is missing.
    name_map, missing, extra = ph.validate_tree_species(tree, ["A.fa", "b.faa", "Z"])
    assert name_map["A.fa"] == "A"
    assert name_map["b.faa"] == "B"
    assert "Z" in missing


def test_presence_from_gene_numbers():
    dSpecies = {3: "A", 4: "B", 5: "C"}
    dGeneNumbers = {
        "hogX": [1, 0, 2],   # A and C present
        "hogY": [0, 0, 0],   # absent everywhere -> dropped
    }
    presence = ph.presence_from_gene_numbers(dGeneNumbers, dSpecies)
    assert presence["hogX"] == {"A", "C"}
    assert "hogY" not in presence


def test_analyze_writes_all_outputs(tmp_path):
    tree_file = tmp_path / "sp.nwk"
    tree_file.write_text(NWK)
    dSpecies = {3: "A", 4: "B", 5: "C", 6: "D", 7: "E"}
    dGeneNumbers = {
        "HOG1": [1, 1, 1, 1, 1],   # core
        "HOG2": [1, 0, 1, 0, 0],   # A,C
        "HOG3": [0, 0, 1, 1, 0],   # C,D
    }
    paths = ph.analyze(dGeneNumbers, dSpecies, str(tree_file), str(tmp_path), "t_")
    for key in ("lca", "gainloss", "nodes", "tree"):
        assert os.path.exists(paths[key]), f"missing output: {key}"

    import pandas as pd
    nodes = pd.read_csv(paths["nodes"], sep="\t")
    # total gains across the summary equals number of HOGs with >=1 carrier (3)
    assert int(nodes["Gains"].sum()) == 3
    annotated = tree_file.parent / "t_species_tree_annotated.nwk"
    assert "N" in annotated.read_text()  # internal node names were written


# ---- Tier 2.1: PD-weighted classification ----

def test_pd_weighted_classification():
    tree = make_tree()
    below = ph.tips_below(tree)
    presence = {
        "h_all": {"A", "B", "C", "D", "E"},   # spans the whole tree -> core
        "h_cd": {"C", "D"},                    # PD 2/9 ~ 0.22 -> shell
        "h_a": {"A"},                          # single tip, PD 0 -> private
    }
    rows, counts = ph.pd_weighted_classification(
        tree, presence, core_threshold=0.9, private_threshold=0.1, below=below)
    by = {r["HOG"]: r for r in rows}
    assert by["h_all"]["Weighted_Class"] == "core"
    assert by["h_cd"]["Weighted_Class"] == "shell"
    assert by["h_a"]["Weighted_Class"] == "private"
    assert counts == {"core": 1, "shell": 1, "private": 1}


# ---- Tier 2.2: clade-conditioned compartments ----

def test_clade_compartments_counts_within_each_clade():
    tree = make_tree()
    below = ph.tips_below(tree)
    presence = {
        "h_all": {"A", "B", "C", "D", "E"},
        "h_ab": {"A", "B"},
        "h_a": {"A"},
        "h_cde": {"C", "D", "E"},
    }
    rows = ph.clade_compartments(tree, presence, below)
    by_node = {r["Node"]: r for r in rows}
    # Postorder naming: (A,B)=N0, (C,D)=N1, (CD,E)=N2, root=N3
    ab = by_node["N0"]
    assert ab["Num_Tips"] == 2
    assert (ab["Clade_Core"], ab["Clade_Shell"], ab["Clade_Private"]) == (2, 0, 1)
    assert ab["HOGs_Present"] == 3
    root = by_node["N3"]
    assert root["Num_Tips"] == 5
    assert (root["Clade_Core"], root["Clade_Shell"], root["Clade_Private"]) == (1, 2, 1)
    assert root["HOGs_Present"] == 4


def test_analyze_writes_weighted_outputs(tmp_path):
    tree_file = tmp_path / "sp.nwk"
    tree_file.write_text(NWK)
    dSpecies = {3: "A", 4: "B", 5: "C", 6: "D", 7: "E"}
    dGeneNumbers = {
        "HOG1": [1, 1, 1, 1, 1],   # spans whole tree -> core
        "HOG2": [1, 0, 0, 0, 0],   # single tip -> private
    }
    paths = ph.analyze(dGeneNumbers, dSpecies, str(tree_file), str(tmp_path), "t_",
                       pan_weighted=True)
    assert os.path.exists(paths["weighted"])
    assert os.path.exists(paths["clade"])
    import pandas as pd
    w = pd.read_csv(paths["weighted"], sep="\t")
    cls = dict(zip(w["HOG"], w["Weighted_Class"]))
    assert cls["HOG1"] == "core"
    assert cls["HOG2"] == "private"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
