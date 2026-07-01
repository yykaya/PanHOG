"""
End-to-end test for run_kaks_pipeline (the rewritten dN/dS orchestration).

MAFFT is not required here: the synthetic proteins are already equal length,
so a stub aligner that passes the input through stands in for the real MSA
step. This exercises the full path HOG selection -> codon back-translation ->
pairwise dN/dS -> output files.

Run with:  pytest -q tests/test_kaks_integration.py
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

PanHOG = pytest.importorskip("PanHOG")
if not PanHOG.HAS_BIOPYTHON or not PanHOG.HAS_DNDS:
    pytest.skip("Biopython / panhog_dnds not available", allow_module_level=True)


# A 20-codon ancestral CDS (no stop codons), then four diverged copies.
ANCESTOR = ("ATG AAA CGT GAT CTT GGG CCC TGG GCT AAT "
            "CAT TTT GAA GTT ACC TCA CTG CGC GGA TAC").replace(" ", "")

# (position, new_base) edits per species -> a controlled mix of changes.
SPECIES_EDITS = {
    "spA": [],                              # identical to ancestor
    "spB": [(5, "G"), (14, "C")],           # a couple of substitutions
    "spC": [(8, "T"), (20, "A"), (33, "G")],
    "spD": [(2, "G"), (44, "T")],
}


def _mutate(seq, edits):
    s = list(seq)
    for pos, base in edits:
        s[pos] = base
    return "".join(s)


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")


@pytest.fixture
def synthetic_project(tmp_path):
    from Bio.Seq import Seq

    prot_dir = tmp_path / "proteins"
    cds_dir = tmp_path / "cds"
    out_dir = tmp_path / "results"
    prot_dir.mkdir(); cds_dir.mkdir(); out_dir.mkdir()

    species = list(SPECIES_EDITS)
    gene_ids = {}
    for sp in species:
        cds = _mutate(ANCESTOR, SPECIES_EDITS[sp])
        gid = f"{sp}_g1"
        gene_ids[sp] = gid
        prot = str(Seq(cds).translate())
        assert "*" not in prot, f"{sp} CDS translated to a stop codon"
        _write_fasta(prot_dir / f"{sp}.fa", [(gid, prot)])
        _write_fasta(cds_dir / f"{sp}.fna", [(gid, cds)])

    # One core HOG containing exactly one gene per species.
    hog_tsv = tmp_path / "N0.tsv"
    with open(hog_tsv, "w") as fh:
        fh.write("HOG\tOG\tGene Tree Parent Clade\t" + "\t".join(species) + "\n")
        fh.write("N0.HOG0000001\tOG0000001\tn0\t" +
                 "\t".join(gene_ids[sp] for sp in species) + "\n")

    return {
        "hog": str(hog_tsv), "prot": str(prot_dir),
        "cds": str(cds_dir), "out": str(out_dir), "species": species,
    }


def _stub_aligner(input_fasta, output_aln, *a, **k):
    """Equal-length proteins are already 'aligned'; just copy them through."""
    with open(input_fasta) as fin, open(output_aln, "w") as fout:
        fout.write(fin.read())
    return True


def test_kaks_pipeline_biopython_end_to_end(synthetic_project, monkeypatch):
    monkeypatch.setattr(PanHOG, "run_alignment", _stub_aligner)

    dGeneNumbers, dHOGs, dSpecies, ddHOGs = PanHOG.parseHOGs(synthetic_project["hog"])
    assert len(dSpecies) == 4
    assert dGeneNumbers["N0.HOG0000001"] == [1, 1, 1, 1]  # core

    PanHOG.run_kaks_pipeline(
        "core", "biopython", synthetic_project["cds"], synthetic_project["prot"],
        dGeneNumbers, ddHOGs, dSpecies, synthetic_project["out"], "test_",
        aligner="mafft", backtrans="naive", reference=None, model="NG86",
    )

    summary = os.path.join(synthetic_project["out"], "test_kaks_results_core.tsv")
    pairwise = os.path.join(synthetic_project["out"], "test_kaks_pairwise_core.tsv")
    assert os.path.exists(summary), "per-HOG summary not written"
    assert os.path.exists(pairwise), "pairwise table not written"

    import pandas as pd
    df_sum = pd.read_csv(summary, sep="\t")
    df_pair = pd.read_csv(pairwise, sep="\t")

    assert list(df_sum["HOG"]) == ["N0.HOG0000001"]
    # 4 sequences -> C(4,2) = 6 pairwise comparisons.
    assert len(df_pair) == 6
    assert int(df_sum.iloc[0]["N_Pairs"]) >= 1
    # dN and dS must be real, finite, non-negative numbers (the old code's bug).
    dN = float(df_sum.iloc[0]["dN"])
    dS = float(df_sum.iloc[0]["dS"])
    assert dN >= 0 and dS >= 0


def test_kaks_pipeline_reference_mode(synthetic_project, monkeypatch):
    monkeypatch.setattr(PanHOG, "run_alignment", _stub_aligner)
    dGeneNumbers, dHOGs, dSpecies, ddHOGs = PanHOG.parseHOGs(synthetic_project["hog"])

    PanHOG.run_kaks_pipeline(
        "core", "biopython", synthetic_project["cds"], synthetic_project["prot"],
        dGeneNumbers, ddHOGs, dSpecies, synthetic_project["out"], "ref_",
        aligner="mafft", backtrans="naive", reference="spA", model="NG86",
    )

    import pandas as pd
    df_pair = pd.read_csv(
        os.path.join(synthetic_project["out"], "ref_kaks_pairwise_core.tsv"), sep="\t")
    # Reference spA vs the other 3 species -> exactly 3 pairs, all involving spA.
    assert len(df_pair) == 3
    involved = set(df_pair["Species1"]) | set(df_pair["Species2"])
    assert "spA" in involved
    for _, row in df_pair.iterrows():
        assert "spA" in (row["Species1"], row["Species2"])


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
