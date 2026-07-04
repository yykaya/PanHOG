"""
Tests for panhog_dnds — the corrected dN/dS (Ka/Ks) engine.

These lock in the behaviour that the previous hand-rolled Nei-Gojobori routine
got wrong: real synonymous/non-synonymous site counting (delegated to
BioPython), and honest NaNs instead of fabricated ratios when a quantity is
undefined.

Run with:  pytest -q   (or)   python tests/test_dnds.py
"""

import math
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pd = pytest.importorskip("panhog_dnds")
if not pd.HAS_BIOPYTHON:
    pytest.skip("Biopython not installed", allow_module_level=True)


def test_identical_sequences_have_no_bogus_ratio():
    s = "ATGAAACGTGATTTTGGGCCCTGG"
    dN, dS, omega = pd.dnds_pair(s, s, "NG86")
    assert dN == 0.0
    assert dS == 0.0
    # dS == 0 -> ratio undefined, must be NaN (never 0, inf, or a made-up value)
    assert math.isnan(omega)


def test_purely_synonymous_changes():
    # Same protein, several silent 3rd-position changes -> dN≈0, dS>0.
    a = "ATGAAACGTGATCTTGGGCCCTGG"
    b = "ATGAAACGTGATCTCGGACCATGG"
    dN, dS, omega = pd.dnds_pair(a, b, "NG86")
    assert dN == pytest.approx(0.0, abs=1e-9)
    assert dS > 0.0
    assert omega == pytest.approx(0.0, abs=1e-9)


def test_purely_nonsynonymous_changes():
    # Amino-acid-changing substitutions -> dN>0. dS may be 0 (ratio NaN).
    a = "ATGAAACGTGATCTTGGGCCCTGG"
    b = "ATGGAAGCTAATCATTGGCCCTGG"
    dN, dS, omega = pd.dnds_pair(a, b, "NG86")
    assert dN > 0.0
    if dS == 0.0:
        assert math.isnan(omega)


def test_clean_pair_strips_gaps_and_stops():
    g1 = "ATG---CGTGATTAA"   # gap codon + trailing stop
    g2 = "ATGAAACGTGATGGG"
    c1, c2 = pd._clean_pair(g1, g2)
    # gap codon and stop codon removed; lengths equal & in-frame
    assert c1 == "ATGCGTGAT"
    assert c2 == "ATGCGTGAT"
    assert len(c1) % 3 == 0 and len(c1) == len(c2)


def test_reference_mode_pairs_only_reference_vs_rest():
    recs = [("g1", "ATGAAA"), ("g2", "ATGAAA"),
            ("g3", "ATGAAA"), ("r1", "ATGAAA")]
    species_of = {"g1": "A", "g2": "B", "g3": "C", "r1": "REF"}
    rows = pd.dnds_from_alignment(
        recs, method="biopython", model="NG86",
        species_of=species_of, reference="REF",
    )
    pairs = {(r["Seq1"], r["Seq2"]) for r in rows}
    assert pairs == {("r1", "g1"), ("r1", "g2"), ("r1", "g3")}


def test_all_pairs_mode_uses_every_combination():
    recs = [("g1", "ATGAAA"), ("g2", "ATGAAA"), ("g3", "ATGAAA")]
    rows = pd.dnds_from_alignment(recs, method="biopython", model="NG86")
    assert len(rows) == 3  # C(3,2)


def test_matches_biopython_reference_values():
    # The engine must agree with BioPython's own cal_dn_ds (it wraps it).
    from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
    a = "ATGAAACGTGATCTTGGGCCCTGGGCTAATCATTTT"
    b = "ATGAAGCGTGATCTCGGACCATGGGCAAACCACTTT"
    exp_dN, exp_dS = cal_dn_ds(CodonSeq(a), CodonSeq(b), method="NG86")
    dN, dS, _ = pd.dnds_pair(a, b, "NG86")
    assert dN == pytest.approx(pd._finite(exp_dN), rel=1e-6, nan_ok=True)
    assert dS == pytest.approx(pd._finite(exp_dS), rel=1e-6, nan_ok=True)


def test_summarize_hog_ignores_nans():
    rows = [
        {"dN": 0.1, "dS": 0.2, "dN_dS": 0.5},
        {"dN": 0.3, "dS": 0.6, "dN_dS": 0.5},
        {"dN": float("nan"), "dS": float("nan"), "dN_dS": float("nan")},
    ]
    s = pd.summarize_hog(rows)
    assert s["dN"] == pytest.approx(0.2)
    assert s["dS"] == pytest.approx(0.4)
    assert s["dN_dS"] == pytest.approx(0.5)
    assert s["n_pairs"] == 2


def test_summarize_hog_uses_ratio_of_means_not_mean_of_ratios():
    # One pair has a near-zero dS, which blows its per-pair ratio up. The robust
    # per-HOG dN_dS (ratio of means) must stay sane; the mean of ratios must not
    # be used as the primary value.
    rows = [
        {"dN": 0.02, "dS": 0.20, "dN_dS": 0.10},
        {"dN": 0.02, "dS": 0.001, "dN_dS": 20.0},
    ]
    s = pd.summarize_hog(rows)
    assert s["dN_dS"] == pytest.approx(0.02 / 0.1005, abs=1e-3)   # ratio of means
    assert s["dN_dS"] < 1                                          # sane
    assert s["dN_dS_mean_of_ratios"] == pytest.approx(10.05, abs=0.1)  # kept, but not used


def test_available_models_reports_scipy_state():
    models = pd.available_models()
    assert "NG86" in models and "LWL85" in models
    if pd.HAS_SCIPY:
        assert "YN00" in models and "ML" in models
    else:
        assert "YN00" not in models


if __name__ == "__main__":
    # Allow running without pytest as a quick smoke test.
    import warnings
    warnings.filterwarnings("ignore")
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
            passed += 1
        except Exception as e:  # noqa: BLE001
            print(f"FAIL  {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{passed}/{len(fns)} tests passed")
    sys.exit(0 if passed == len(fns) else 1)
