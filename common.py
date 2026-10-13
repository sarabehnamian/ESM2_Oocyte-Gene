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
