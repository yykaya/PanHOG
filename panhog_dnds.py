#!/usr/bin/env python3
"""
panhog_dnds.py
==============
Correct codon-based dN/dS (Ka/Ks) analysis for PanHOG.

Background
----------
The previous built-in "Nei-Gojobori (1986)" routine in ``PanHOG.py``
(``calculate_ng86``) was mathematically incorrect: it counted synonymous and
non-synonymous *sites* with arbitrary constants (``S_sites += 1``,
``N_sites += 2.5``) and only on codons that already differed, so the pS/pN
denominators were meaningless and the resulting Ka, Ks and Ka/Ks values could
not be interpreted. This module replaces that logic with trustworthy engines.

Engines
-------
* ``biopython``      -- ``Bio.codonalign.codonseq.cal_dn_ds`` (validated).
                        Sub-models via ``model=``: ``NG86`` (default), ``LWL85``,
                        ``YN00`` and ``ML`` (the last two require SciPy).
* ``codeml``         -- PAML ``codeml`` run in pairwise mode (``runmode = -2``),
                        driven through ``Bio.Phylo.PAML.codeml``. This mirrors
                        the PAML-based dN-dS-Snakemake-pipeline workflow.
* ``kakscalculator`` -- external ``KaKs_Calculator`` (AXT input); handled by the
                        caller, kept here only for the model listing.

All engines operate on an already-built *codon alignment* (protein alignment
back-translated to nucleotides). Gap and stop columns are removed pairwise
(complete deletion) before the substitution model is applied, which is the
standard treatment and keeps the library engines well-defined.
"""

from __future__ import annotations

import math
import os
import shutil
import tempfile
from itertools import combinations

try:
    from Bio.Seq import Seq
    from Bio.Data import CodonTable
    HAS_BIOPYTHON = True
except ImportError:  # pragma: no cover - exercised only without Biopython
    HAS_BIOPYTHON = False

try:
    import scipy  # noqa: F401  (presence check only)
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# Models that never need SciPy vs. those that do.
_MODELS_NO_SCIPY = ("NG86", "LWL85")
_MODELS_NEED_SCIPY = ("YN00", "ML")

# BioPython uses -1.0 as an "undefined / saturated" sentinel for dN or dS.
_UNDEFINED = -1.0


def available_models():
    """Return the dN/dS sub-models usable in this environment."""
    models = list(_MODELS_NO_SCIPY)
    if HAS_SCIPY:
        models += list(_MODELS_NEED_SCIPY)
    return models


def _clean_pair(seq1, seq2, codon_table_id=1):
    """
    Return a gap/stop-free, equal-length, in-frame codon pair.

    Columns (codons) are dropped when either sequence has a gap, an ambiguous
    base, an incomplete codon or a stop codon. This is pairwise complete
    deletion, the conventional input treatment for NG86/LWL85/YN00.
    """
    if len(seq1) != len(seq2):
        n = min(len(seq1), len(seq2))
        seq1, seq2 = seq1[:n], seq2[:n]

    table = CodonTable.unambiguous_dna_by_id[codon_table_id]
    stops = set(table.stop_codons)
    valid = set("ACGT")

    keep1, keep2 = [], []
    for i in range(0, len(seq1) - 2, 3):
        c1 = seq1[i:i + 3].upper()
        c2 = seq2[i:i + 3].upper()
        if len(c1) < 3 or len(c2) < 3:
            continue
        if not (set(c1) <= valid and set(c2) <= valid):
            continue
        if c1 in stops or c2 in stops:
            continue
        keep1.append(c1)
        keep2.append(c2)
    return "".join(keep1), "".join(keep2)


def _finite(value):
    """Map BioPython's -1 sentinel (and NaN/inf) to Python float('nan')."""
    if value is None:
        return float("nan")
    try:
        v = float(value)
    except (TypeError, ValueError):
        return float("nan")
    if v == _UNDEFINED or math.isinf(v) or math.isnan(v):
        return float("nan")
    return v


def dnds_pair(seq1, seq2, model="NG86", codon_table_id=1):
    """
    Compute (dN, dS, omega) for one aligned codon pair with BioPython.

    Returns floats; any undefined quantity is ``float('nan')`` (never a bogus
    number or a -1 sentinel). ``omega`` is dN/dS, or NaN when dS is 0/undefined.
    """
    if not HAS_BIOPYTHON:
        raise RuntimeError("Biopython is required for the 'biopython' dN/dS engine.")

    model = (model or "NG86").upper()
    if model in _MODELS_NEED_SCIPY and not HAS_SCIPY:
        raise RuntimeError(
            f"dN/dS model '{model}' requires SciPy. Install scipy or use NG86/LWL85."
        )

    from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds

    c1, c2 = _clean_pair(str(seq1), str(seq2), codon_table_id)
    # Need at least a few shared codons for a meaningful estimate.
    if len(c1) < 6 or len(c1) % 3 != 0:
        return float("nan"), float("nan"), float("nan")

    try:
        dN, dS = cal_dn_ds(CodonSeq(c1), CodonSeq(c2), method=model)
    except (ValueError, ZeroDivisionError, KeyError):
        # e.g. LWL85 "math domain error" under saturation.
        return float("nan"), float("nan"), float("nan")

    dN = _finite(dN)
    dS = _finite(dS)
    if not math.isnan(dS) and dS > 0 and not math.isnan(dN):
        omega = dN / dS
    else:
        omega = float("nan")
    return dN, dS, omega


def dnds_pair_codeml(seq1, seq2, id1="seq1", id2="seq2",
                     codeml_path="codeml", codon_table_id=1, workdir=None):
    """
    Compute (dN, dS, omega) for one aligned codon pair with PAML ``codeml``.

    Uses runmode = -2 (pairwise) via ``Bio.Phylo.PAML.codeml``. Raises
    ``FileNotFoundError`` when the ``codeml`` executable is not on PATH so the
    caller can fall back to the ``biopython`` engine with a clear message.
    """
    if not HAS_BIOPYTHON:
        raise RuntimeError("Biopython is required to drive PAML codeml.")
    if shutil.which(codeml_path) is None:
        raise FileNotFoundError(
            f"PAML 'codeml' not found (looked for '{codeml_path}'). "
            f"Install PAML (e.g. `conda install -c bioconda paml`)."
        )

    from Bio.Phylo.PAML import codeml

    c1, c2 = _clean_pair(str(seq1), str(seq2), codon_table_id)
    if len(c1) < 6 or len(c1) % 3 != 0:
        return float("nan"), float("nan"), float("nan")

    created = False
    if workdir is None:
        workdir = tempfile.mkdtemp(prefix="panhog_codeml_")
        created = True
    try:
        n_codons = len(c1) // 3
        # PHYLIP sequential; codeml is tolerant of long names on one line.
        phy = os.path.join(workdir, "pair.phy")
        with open(phy, "w") as fh:
            fh.write(f" 2 {len(c1)}\n")
            fh.write(f"{id1[:30]:<32}{c1}\n")
            fh.write(f"{id2[:30]:<32}{c2}\n")
        # Minimal 2-taxon tree (not used by runmode -2, but codeml wants a file).
        tree = os.path.join(workdir, "pair.tree")
        with open(tree, "w") as fh:
            fh.write(f"({id1[:30]},{id2[:30]});\n")
        out = os.path.join(workdir, "pair.codeml")

        cml = codeml.Codeml(alignment=phy, tree=tree, out_file=out, working_dir=workdir)
        cml.set_options(
            runmode=-2, seqtype=1, CodonFreq=2, model=0, NSsites=[0],
            icode=codon_table_id - 1, fix_kappa=0, kappa=2,
            fix_omega=0, omega=0.4, cleandata=1, verbose=0,
        )
        results = cml.run(command=codeml_path, verbose=False, parse=True)

        pair = results.get("pairwise", {})
        rec = pair.get(id1[:30], {}).get(id2[:30]) or pair.get(id2[:30], {}).get(id1[:30])
        if not rec:
            return float("nan"), float("nan"), float("nan")
        dN = _finite(rec.get("dN"))
        dS = _finite(rec.get("dS"))
        omega = _finite(rec.get("omega"))
        if math.isnan(omega) and not math.isnan(dS) and dS > 0 and not math.isnan(dN):
            omega = dN / dS
        return dN, dS, omega
    finally:
        if created:
            shutil.rmtree(workdir, ignore_errors=True)


def _select_pairs(ids, species_of=None, reference=None):
    """
    Yield (a, b) id pairs. With a reference species, pair every non-reference
    sequence against each reference sequence; otherwise use all unordered pairs.
    """
    if reference and species_of:
        ref_ids = [i for i in ids if species_of.get(i) == reference]
        qry_ids = [i for i in ids if species_of.get(i) != reference]
        if ref_ids:
            for r in ref_ids:
                for q in qry_ids:
                    yield r, q
            return
    yield from combinations(ids, 2)


def dnds_from_alignment(records, method="biopython", model="NG86",
                        species_of=None, reference=None,
                        codeml_path="codeml", codon_table_id=1, workdir=None):
    """
    Compute pairwise dN/dS across an aligned codon set.

    Parameters
    ----------
    records : list[(id, aligned_codon_seq)]
    method  : 'biopython' or 'codeml'
    model   : sub-model for the biopython engine (NG86/LWL85/YN00/ML)
    species_of : optional {gene_id: species} for reference-based pairing
    reference  : optional species name; restrict pairs to reference-vs-rest

    Returns
    -------
    list[dict] with keys: Seq1, Seq2, Species1, Species2, dN, dS, dN_dS.
    """
    ids = [r[0] for r in records]
    seqof = {rid: rseq for rid, rseq in records}
    rows = []
    for a, b in _select_pairs(ids, species_of, reference):
        if method == "codeml":
            try:
                dN, dS, omega = dnds_pair_codeml(
                    seqof[a], seqof[b], a, b, codeml_path, codon_table_id, workdir
                )
            except FileNotFoundError:
                raise  # let the caller decide to fall back
        else:
            dN, dS, omega = dnds_pair(seqof[a], seqof[b], model, codon_table_id)

        rows.append({
            "Seq1": a,
            "Seq2": b,
            "Species1": (species_of or {}).get(a, ""),
            "Species2": (species_of or {}).get(b, ""),
            "dN": dN,
            "dS": dS,
            "dN_dS": omega,
        })
    return rows


def summarize_hog(rows):
    """
    Collapse per-pair rows for one HOG into a mean summary, ignoring NaNs.
    Returns dict with mean dN, dS, dN/dS and the number of valid pairs.
    """
    def _mean(vals):
        good = [v for v in vals if isinstance(v, float) and not math.isnan(v)]
        return (sum(good) / len(good)) if good else float("nan")

    dN = _mean([r["dN"] for r in rows])
    dS = _mean([r["dS"] for r in rows])
    ratios = [r["dN_dS"] for r in rows]
    good_ratio = [v for v in ratios if isinstance(v, float) and not math.isnan(v)]
    return {
        "dN": dN,
        "dS": dS,
        "dN_dS": (sum(good_ratio) / len(good_ratio)) if good_ratio else (
            dN / dS if not math.isnan(dS) and dS > 0 else float("nan")
        ),
        "n_pairs": len(good_ratio),
    }
