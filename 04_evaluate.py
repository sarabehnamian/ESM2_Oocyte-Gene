#!/usr/bin/env python3
"""
04_evaluate.py

AUROC per gene and per panel, with bootstrap confidence intervals.

The comparison: does ESM-2 separate pathogenic from population variants
as well on oocyte / IVF-failure genes as it does on deeply conserved
control genes?

Intervals, not point estimates. Some genes have a few dozen variants, so
a bare AUROC would be meaningless. Genes below the minimum class count
are reported but excluded from the panel summary.

The panel difference is tested by a stratified bootstrap over genes and
variants, which gives an interval on the gap. If that interval crosses
zero, there is no detectable difference, and that is the result.

Output
    results/auroc_by_gene.jsonl
    results/panel_summary.jsonl

Usage
    python 04_evaluate.py
    python 04_evaluate.py --min-class 10
"""

import argparse
import os
from collections import defaultdict

import numpy as np

from common import (RESULT_DIR, auroc, bootstrap_auroc, fmt_ci, read_jsonl,
                    write_jsonl)


def panel_bootstrap(genes_data, n_boot=2000, seed=0):
    """
    Resample genes, then variants within each gene. Two-level, because
    genes are the unit we generalise over, not variants.
    Returns (point, lo, hi) for the mean per-gene AUROC.
    """
    rng = np.random.default_rng(seed)
    names = list(genes_data.keys())
    if not names:
        return None, None, None

    point_vals = [auroc(*genes_data[g]) for g in names]
    point_vals = [v for v in point_vals if v is not None]
    if not point_vals:
        return None, None, None
    point = float(np.mean(point_vals))

    draws = []
    for _ in range(n_boot):
        picked = rng.choice(len(names), len(names), replace=True)
        vals = []
        for i in picked:
            scores, labels = genes_data[names[i]]
            scores = np.asarray(scores)
            labels = np.asarray(labels)
            pos = np.flatnonzero(labels == 1)
            neg = np.flatnonzero(labels == 0)
            idx = np.concatenate([rng.choice(pos, len(pos), replace=True),
                                  rng.choice(neg, len(neg), replace=True)])
            a = auroc(scores[idx], labels[idx])
            if a is not None:
                vals.append(a)
        if vals:
            draws.append(float(np.mean(vals)))

    if not draws:
        return point, None, None
    return point, float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def gap_bootstrap(oo, ct, n_boot=2000, seed=1):
    """Interval on (control mean AUROC) minus (oocyte mean AUROC)."""
    rng = np.random.default_rng(seed)
    oo_names, ct_names = list(oo.keys()), list(ct.keys())
    if not oo_names or not ct_names:
        return None, None, None

    def mean_auroc(names, data, gen):
        vals = []
        for i in gen.choice(len(names), len(names), replace=True):
            scores, labels = data[names[i]]
            scores, labels = np.asarray(scores), np.asarray(labels)
            pos = np.flatnonzero(labels == 1)
            neg = np.flatnonzero(labels == 0)
            idx = np.concatenate([gen.choice(pos, len(pos), replace=True),
                                  gen.choice(neg, len(neg), replace=True)])
            a = auroc(scores[idx], labels[idx])
            if a is not None:
                vals.append(a)
        return float(np.mean(vals)) if vals else None

    base_oo = float(np.mean([auroc(*oo[g]) for g in oo_names]))
    base_ct = float(np.mean([auroc(*ct[g]) for g in ct_names]))
    point = base_ct - base_oo

    draws = []
    for _ in range(n_boot):
        a = mean_auroc(oo_names, oo, rng)
        b = mean_auroc(ct_names, ct, rng)
        if a is not None and b is not None:
            draws.append(b - a)
    if not draws:
        return point, None, None
    return point, float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-class", type=int, default=5,
                    help="minimum variants in EACH class to include a gene")
    ap.add_argument("--n-boot", type=int, default=2000)
    args = ap.parse_args()

    scored = read_jsonl(os.path.join(RESULT_DIR, "scores.jsonl"))
    by_gene = defaultdict(list)
    for r in scored:
        by_gene[r["gene"]].append(r)

    rows = []
    usable = {"esm": {"oocyte": {}, "control": {}},
              "blosum": {"oocyte": {}, "control": {}}}

    print("%-10s %-8s %5s %6s  %-22s %-22s" %
          ("gene", "panel", "path", "benign", "ESM-2 AUROC [95% CI]",
           "BLOSUM62 AUROC"))

    for gene in sorted(by_gene):
        g = by_gene[gene]
        panel = g[0]["panel"]
        labels = [r["label"] for r in g]
        esm = [r["esm_score"] for r in g]
        blo = [r["blosum_score"] for r in g]
        n_pos, n_neg = labels.count(1), labels.count(0)

        e_pt, e_lo, e_hi = bootstrap_auroc(esm, labels, args.n_boot)
        b_pt, b_lo, b_hi = bootstrap_auroc(blo, labels, args.n_boot)

        included = n_pos >= args.min_class and n_neg >= args.min_class
        if included:
            usable["esm"][panel][gene] = (esm, labels)
            usable["blosum"][panel][gene] = (blo, labels)

        rows.append({
            "gene": gene, "panel": panel, "n_pathogenic": n_pos, "n_benign": n_neg,
            "esm_auroc": e_pt, "esm_lo": e_lo, "esm_hi": e_hi,
            "blosum_auroc": b_pt, "blosum_lo": b_lo, "blosum_hi": b_hi,
            "included": included,
        })
        print("%-10s %-8s %5d %6d  %-22s %-22s%s"
              % (gene, panel, n_pos, n_neg, fmt_ci(e_pt, e_lo, e_hi),
                 fmt_ci(b_pt, b_lo, b_hi), "" if included else "   (excluded)"))

    write_jsonl(os.path.join(RESULT_DIR, "auroc_by_gene.jsonl"), rows)

    print("\npanel summary, mean AUROC over included genes")
    summary = []
    for method in ("esm", "blosum"):
        for panel in ("oocyte", "control"):
            data = usable[method][panel]
            pt, lo, hi = panel_bootstrap(data, args.n_boot)
            summary.append({"method": method, "panel": panel,
                            "n_genes": len(data), "auroc": pt,
                            "lo": lo, "hi": hi})
            print("  %-7s %-8s  n_genes=%2d  %s"
                  % (method, panel, len(data), fmt_ci(pt, lo, hi)))

    for method in ("esm", "blosum"):
        pt, lo, hi = gap_bootstrap(usable[method]["oocyte"],
                                   usable[method]["control"], args.n_boot)
        summary.append({"method": method, "panel": "control_minus_oocyte",
                        "auroc": pt, "lo": lo, "hi": hi, "n_genes": None})
        if pt is None:
            print("\n  %s gap: not computable" % method)
            continue
        crosses = (lo is not None and lo <= 0.0 <= hi)
        print("\n  %s gap (control minus oocyte): %s  %s"
              % (method, fmt_ci(pt, lo, hi),
                 "-> interval crosses zero, NO detectable difference"
                 if crosses else "-> interval excludes zero"))

    write_jsonl(os.path.join(RESULT_DIR, "panel_summary.jsonl"), summary)
    print("\nwrote results/auroc_by_gene.jsonl and results/panel_summary.jsonl")
    print("\nRead the gap interval before writing anything. If it crosses zero, "
          "the hypothesis is not supported and the honest result is negative.")


if __name__ == "__main__":
    main()
