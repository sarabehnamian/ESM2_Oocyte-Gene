"""
Shared configuration and helpers.

Not a pipeline step. Imported by the numbered scripts.

The question: do protein language models work on the genes where IVF
actually fails?

Oocyte / IVF-failure genes are maternal-effect genes expressed only in
the egg and early embryo. TUBB8 is primate-specific. Protein language
models learn from conservation across species, so on genes with little
evolutionary depth there may be little signal to learn from. The control
panel is deeply conserved disease genes where these models are known to
work well.
"""

import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT, "data")
RESULT_DIR = os.path.join(ROOT, "results")
FIG_DIR = os.path.join(ROOT, "figures")
for _d in (DATA_DIR, RESULT_DIR, FIG_DIR):
    os.makedirs(_d, exist_ok=True)


# Oocyte maturation arrest, fertilization failure, early embryonic arrest.
# The phenotype for most of these is only observable in an IVF clinic.
OOCYTE_GENES = [
    "TUBB8", "PATL2", "WEE2", "PADI6", "NLRP5", "TLE6",
    "ZP1", "ZP2", "ZP3", "ZP4", "KHDC3L", "NLRP7",
    "TRIP13", "CDC20", "BTG4", "PANX1", "REC114", "MEI1",
]

# Deeply conserved, heavily studied disease genes. Protein language models
# are known to perform well here, so this is the reference ceiling.
CONTROL_GENES = [
    "BRCA1", "BRCA2", "CFTR", "TP53", "MLH1",
    "PTEN", "LDLR", "MSH2", "ATM", "RB1",
]

ALL_GENES = OOCYTE_GENES + CONTROL_GENES
PANEL = {g: "oocyte" for g in OOCYTE_GENES}
PANEL.update({g: "control" for g in CONTROL_GENES})

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def write_jsonl(path, records):
    with open(path, "w", encoding="utf-8") as fh:
        for r in records:
            fh.write(json.dumps(r, sort_keys=True) + "\n")


def read_jsonl(path):
    if not os.path.exists(path):
        raise SystemExit("missing %s -- run the earlier step first" % path)
    out = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def append_jsonl(path, record):
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(record, sort_keys=True) + "\n")


def auroc(scores, labels):
    """
    Area under the ROC curve by rank sum, with ties handled correctly.
    labels: 1 = pathogenic, 0 = benign.
    Returns None when either class is missing.
    """
    scores = np.asarray(scores, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.int64)
    n_pos = int((labels == 1).sum())
    n_neg = int((labels == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return None

    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=np.float64)
    sorted_scores = scores[order]
    i = 0
    while i < len(scores):
        j = i
        while j + 1 < len(scores) and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1

    return float((ranks[labels == 1].sum() - n_pos * (n_pos + 1) / 2.0)
                 / (n_pos * n_neg))


def bootstrap_auroc(scores, labels, n_boot=2000, seed=0):
    """
    Percentile bootstrap interval, stratified by class so neither class
    can vanish from a resample. Small panels need intervals, not points.
    """
    scores = np.asarray(scores, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.int64)
    point = auroc(scores, labels)
    if point is None:
        return None, None, None

    rng = np.random.default_rng(seed)
    pos = np.flatnonzero(labels == 1)
    neg = np.flatnonzero(labels == 0)
    vals = []
    for _ in range(n_boot):
        idx = np.concatenate([rng.choice(pos, len(pos), replace=True),
                              rng.choice(neg, len(neg), replace=True)])
        a = auroc(scores[idx], labels[idx])
        if a is not None:
            vals.append(a)
    if not vals:
        return point, None, None
    return point, float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def fmt_ci(point, lo, hi):
    if point is None:
        return "n/a"
    if lo is None:
        return "%.3f" % point
    return "%.3f [%.3f, %.3f]" % (point, lo, hi)


# ---------------------------------------------------------------------------
# Helpers for the follow-up analyses, steps 06 onwards. Nothing above this
# line changed, so steps 00-05 behave exactly as before.
# ---------------------------------------------------------------------------

import itertools
from collections import defaultdict

# ClinVar review status to gold stars, as ClinVar defines them.
REVIEW_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, multiple submitters": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no interpretation for the single variant": 0,
}


def review_stars(status):
    """Gold stars for a ClinVar review status string. Unknown strings map
    to None so they are counted, not silently treated as zero."""
    if status is None:
        return None
    return REVIEW_STARS.get(status.strip().lower())


def by_gene(records):
    out = defaultdict(list)
    for r in records:
        out[r["gene"]].append(r)
    return out


def gene_aurocs(records, score_keys, min_class=5):
    """
    Point AUROC per gene for each score key. Returns a list of dicts with
    gene, panel, class counts, one AUROC per key and the inclusion flag.
    A record lacking a score for some key is left out of that key's AUROC
    only, and the count of usable records per key is reported.
    """
    rows = []
    for gene, recs in sorted(by_gene(records).items()):
        labels = [r["label"] for r in recs]
        row = {"gene": gene, "panel": recs[0]["panel"],
               "n_pathogenic": labels.count(1), "n_benign": labels.count(0)}
        row["included"] = (row["n_pathogenic"] >= min_class
                           and row["n_benign"] >= min_class)
        for key in score_keys:
            usable = [r for r in recs if r.get(key) is not None]
            row[key + "_n"] = len(usable)
            row[key + "_auroc"] = auroc([r[key] for r in usable],
                                        [r["label"] for r in usable])
        rows.append(row)
    return rows


def exact_panel_permutation(values, is_oocyte):
    """
    Exact permutation test for a panel gap.

    values      one number per gene, e.g. per-gene AUROC
    is_oocyte   booleans, the observed panel labels

    Every way of choosing which genes carry the oocyte label is
    enumerated, keeping the panel sizes fixed. The gap is control mean
    minus oocyte mean. One-sided p = (#splits with gap >= observed) /
    (#splits). The observed split is one of them, so p >= 1/#splits.

    Returns (observed_gap, n_at_least, n_splits, null_gaps).
    """
    values = np.asarray(values, dtype=np.float64)
    is_oocyte = np.asarray(is_oocyte, dtype=bool)
    n, k = len(values), int(is_oocyte.sum())
    observed = float(values[~is_oocyte].mean() - values[is_oocyte].mean())

    null = []
    for chosen in itertools.combinations(range(n), k):
        mask = np.zeros(n, dtype=bool)
        mask[list(chosen)] = True
        null.append(float(values[~mask].mean() - values[mask].mean()))
    null = np.asarray(null)

    # A small tolerance so a split whose gap equals the observed gap is
    # counted even after floating point rounding.
    n_ge = int((null >= observed - 1e-12).sum())
    return observed, n_ge, len(null), null


def auroc_diff_permutation(scores, labels, groups, n_perm=10000, seed=0):
    """
    Permutation test for the difference in AUROC between two groups of
    variants (group True minus group False). Group membership is shuffled
    within each class, so both groups keep their pathogenic and benign
    counts. Two-sided p with the +1 correction.
    """
    scores = np.asarray(scores, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.int64)
    groups = np.asarray(groups, dtype=bool)

    def diff(g):
        a = auroc(scores[g], labels[g])
        b = auroc(scores[~g], labels[~g])
        if a is None or b is None:
            return None
        return a - b

    observed = diff(groups)
    if observed is None:
        return None, None
    rng = np.random.default_rng(seed)
    pos = np.flatnonzero(labels == 1)
    neg = np.flatnonzero(labels == 0)
    hits, done = 0, 0
    for _ in range(n_perm):
        g = groups.copy()
        g[pos] = rng.permutation(groups[pos])
        g[neg] = rng.permutation(groups[neg])
        d = diff(g)
        if d is None:
            continue
        done += 1
        hits += abs(d) >= abs(observed) - 1e-12
    return observed, (hits + 1) / (done + 1)


def spearman_permutation(x, y, n_perm=100000, seed=0):
    """Spearman rho with a two-sided permutation p-value (+1 corrected).
    Small n, so no reliance on the asymptotic t approximation."""
    def ranks(v):
        v = np.asarray(v, dtype=np.float64)
        order = np.argsort(v, kind="mergesort")
        r = np.empty(len(v))
        sv = v[order]
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and sv[j + 1] == sv[i]:
                j += 1
            r[order[i:j + 1]] = 0.5 * (i + j) + 1.0
            i = j + 1
        return r

    rx, ry = ranks(x), ranks(y)
    rho = float(np.corrcoef(rx, ry)[0, 1])
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_perm):
        if abs(np.corrcoef(rx, rng.permutation(ry))[0, 1]) >= abs(rho) - 1e-12:
            hits += 1
    return rho, (hits + 1) / (n_perm + 1)
