"""
Regression test for the per-species indexing in generate_summary_stats:
dGeneNumbers vectors are 0-based but species indices start at 3, so every
species after the first (and especially the last ones) must still get the
correct gene counts.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pytest.importorskip("pandas")
ph = pytest.importorskip("PanHOG")


def test_summary_stats_counts_all_species(tmp_path):
    dSpecies = {3: "spA", 4: "spB", 5: "spC"}   # N0 columns 3,4,5
    dGeneNumbers = {"h1": [1, 2, 3]}            # core; spA=1, spB=2, spC=3
    df = ph.generate_summary_stats(dGeneNumbers, {}, dSpecies, str(tmp_path), "t_")
    row = {r["Species"]: r for _, r in df.iterrows()}
    assert row["spA"]["Core_Gene_Count"] == 1
    assert row["spB"]["Core_Gene_Count"] == 2
    assert row["spC"]["Core_Gene_Count"] == 3   # was 0 before the index fix
    assert row["spC"]["Core_HOG_Count"] == 1


def test_random_hog_matrix_counts_all_species(tmp_path):
    dSpecies = {3: "spA", 4: "spB", 5: "spC"}
    dGeneNumbers = {"h1": [0, 1, 2]}   # spA absent, spB single-copy, spC multi-copy
    ph.generate_random_hog_matrix(dGeneNumbers, dSpecies, str(tmp_path), "t_", n_hogs=1)
    import pandas as pd
    df = pd.read_csv(os.path.join(str(tmp_path), "t_random_hog_matrix.tsv"),
                     sep="\t", index_col=0)  # rows=species, cols=HOGs
    assert int(df.loc["spA", "h1"]) == 0
    assert int(df.loc["spB", "h1"]) == 1
    assert int(df.loc["spC", "h1"]) == 2   # was 0 (Absent) before the index fix


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
