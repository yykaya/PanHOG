"""
Tests for panhog_genetrees — HOG occupancy parsing, collision-safe sequence
loading, tree QC metrics, species-tree concordance, and verdict logic.

The MAFFT/RAxML-NG steps are not unit-tested here (they need external binaries);
the QC/logic layer is. Run with:  pytest -q tests/test_genetrees.py
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

gt = pytest.importorskip("panhog_genetrees")


# ---- occupancy + compartment ----

def test_parse_occupancy_and_compartment(tmp_path):
    n0 = tmp_path / "N0.tsv"
    n0.write_text(
        "HOG\tOG\tGene Tree Parent Clade\tspA\tspB\tspC\n"
        "H_core\tOG0\tn1\tg1\tg2\tg3\n"          # 3/3 -> core
        "H_shell\tOG1\tn1\tg4\t\tg5\n"           # 2/3 -> shell
        "H_priv\tOG2\tn1\t\tg6\t\n"              # 1/3 -> private
    )
    species, occ = gt.parse_hog_occupancy(str(n0))
    assert species == ["spA", "spB", "spC"]
    assert gt.compartment_of(occ["H_core"], 3) == "core"
    assert gt.compartment_of(occ["H_shell"], 3) == "shell"
    assert gt.compartment_of(occ["H_priv"], 3) == "private"
    assert occ["H_shell"] == {"spA": ["g4"], "spC": ["g5"]}


def test_seqstore_no_cross_accession_collision(tmp_path):
    # Same gene ID 'AT1G1' in two accessions with DIFFERENT sequences.
    (tmp_path / "Col.fa").write_text(">AT1G1\nMAAA\n")
    (tmp_path / "Tanz.fa").write_text(">AT1G1\nMBBB\n")
    store = gt.SeqStore(str(tmp_path))
    assert store.get("Col", "AT1G1") == "MAAA"
    assert store.get("Tanz", "AT1G1") == "MBBB"   # not overwritten


# ---- verdict logic (no external tools) ----

def test_verdict_flags_near_identical():
    m = {"n_tips": 12, "mean_patristic": 0.0, "mean_support": 5.0}
    assert gt.make_verdict(m, nrf=1.0).startswith("REVIEW")
    assert "near-identical" in gt.make_verdict(m, nrf=1.0)


def test_verdict_flags_well_supported_conflict():
    m = {"n_tips": 6, "mean_patristic": 0.12, "mean_support": 85.0}
    v = gt.make_verdict(m, nrf=1.0)
    assert v.startswith("REVIEW") and "conflict" in v


def test_verdict_low_signal_is_not_flagged():
    # High nRF but low support -> not a real conflict, just noise.
    m = {"n_tips": 12, "mean_patristic": 0.02, "mean_support": 40.0}
    v = gt.make_verdict(m, nrf=1.0)
    assert v.startswith("OK") and "low phylogenetic signal" in v


def test_verdict_concordant_is_ok():
    m = {"n_tips": 5, "mean_patristic": 0.15, "mean_support": 70.0}
    assert gt.make_verdict(m, nrf=0.0) == "OK"


# ---- tree metrics + concordance (need ete3) ----

def test_tree_metrics(tmp_path):
    pytest.importorskip("ete3")
    tree = tmp_path / "t.nwk"
    tree.write_text("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")
    m = gt.tree_metrics(str(tree))
    assert m["n_tips"] == 4
    assert m["tree_len"] == pytest.approx(0.6, abs=1e-6)
    assert m["mean_patristic"] > 0


def test_species_concordance_identical_and_conflicting(tmp_path):
    pytest.importorskip("ete3")
    species = tmp_path / "sp.nwk"
    species.write_text("((A:1,B:1):1,(C:1,D:1):1);")
    concordant = tmp_path / "g1.nwk"
    concordant.write_text("((A:1,B:1):1,(C:1,D:1):1);")
    conflicting = tmp_path / "g2.nwk"
    conflicting.write_text("((A:1,C:1):1,(B:1,D:1):1);")
    assert gt.species_concordance(str(concordant), str(species)) == 0.0
    assert gt.species_concordance(str(conflicting), str(species)) == 1.0


# ---- CDS back-translation + pep/cds store normalisation ----

def test_seqstore_normalises_pep_and_cds(tmp_path):
    (tmp_path / "Col-0_Chr1.pep.fa").write_text(">g\nMAA\n")
    store = gt.SeqStore(str(tmp_path))
    # An N0 accession is 'Col-0_Chr1.pep'; a '.cds' lookup normalises to the same key.
    assert store.get("Col-0_Chr1.pep", "g") == "MAA"
    assert store.get("Col-0_Chr1.cds", "g") == "MAA"


def test_backtranslate_maps_residues_and_gaps():
    prot_aln = [("t1", "M-K"), ("t2", "MAK")]
    tip_map = {"t1": ("accA", "g1"), "t2": ("accB", "g2")}

    class FakeCDS:
        def get(self, acc, gene):
            return {"g1": "ATGAAA", "g2": "ATGGCTAAA"}[gene]

    out = dict(gt.backtranslate(prot_aln, tip_map, FakeCDS()))
    assert out["t1"] == "ATG---AAA"   # M, gap, K
    assert out["t2"] == "ATGGCTAAA"   # M A K


# ---- compartment-validation status (codon tree as arbiter) ----

def test_status_confirmed():
    s, _ = gt.compute_status(
        {"n_tips": 12, "mean_patristic": 0.1, "mean_support": 40}, 1.0,
        {"n_tips": 12, "mean_patristic": 0.3, "mean_support": 90}, 0.0)
    assert s == "Confirmed"


def test_status_conflict_when_codon_well_supported_but_discordant():
    s, _ = gt.compute_status({}, None,
                             {"n_tips": 6, "mean_patristic": 0.2, "mean_support": 85}, 1.0)
    assert s == "Conflict"


def test_status_redundant_when_near_identical():
    s, _ = gt.compute_status({}, None,
                             {"n_tips": 12, "mean_patristic": 0.0, "mean_support": 5}, 1.0)
    assert s == "Redundant"


def test_status_low_signal_falls_back_to_protein():
    s, _ = gt.compute_status(
        {"n_tips": 12, "mean_patristic": 0.02, "mean_support": 30}, 1.0, {}, None)
    assert s == "Low-signal"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
