"""
Tests for panhog_kaks_compartments — compartment sampling and statistics.
The alignment/dN/dS steps need external tools and are not unit-tested here.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

kc = pytest.importorskip("panhog_kaks_compartments")


def test_sample_by_compartment_filters_by_seqcount():
    occ = {
        "H_core":  {"A": ["a"], "B": ["b"], "C": ["c"]},   # 3 seqs, all 3 -> core
        "H_single": {"A": ["a"]},                           # 1 seq -> too few, skipped
        "H_privmulti": {"A": ["a1", "a2"]},                 # 2 seqs, 1 accession -> private
        "H_big": {"A": ["x" + str(i) for i in range(20)]},  # 20 seqs -> too many, skipped
    }
    b = kc.sample_by_compartment(occ, n_species=3, per_compartment=10,
                                 min_seqs=2, max_seqs=12)
    assert "H_core" in b["core"]
    assert "H_privmulti" in b["private"]
    assert "H_single" not in b["private"]     # single-copy private has no pair
    assert "H_big" not in b["core"]


def test_compartment_stats_medians_and_kruskal():
    data = {"core": [0.1, 0.15, 0.2],
            "shell": [0.3, 0.4, 0.5],
            "private": [0.5, 0.6, 0.7]}
    s = kc.compartment_stats(data)
    assert s["core_median"] == 0.15
    assert s["shell_median"] == 0.4
    if kc.HAS_SCIPY:
        assert "kruskal_p" in s and 0.0 <= s["kruskal_p"] <= 1.0


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
