#!/usr/bin/env python3
"""
07_review_status.py

Does the panel gap come from ClinVar curation depth rather than biology?

Control gene positives are mostly 2+ star ClinVar records (multiple
submitters agreeing, expert panels, practice guidelines). Oocyte gene
positives are almost all 0-1 star (single submitter or no assertion
criteria). If well-reviewed positives are simply easier to separate,
for example because in-silico and conservation evidence fed into their
classification, then comparing panels confounds biology with review
depth.

This step
  1. reports the share of 2+ star positives per panel,
  2. computes per-gene AUROC (ESM-2, BLOSUM62) with positives split into
     0-1 star and 2+ star, each against all of that gene's negatives,
  3. repeats the panel comparison keeping only 0-1 star positives on
     both sides, with the exact permutation test from step 06 and the
     same minimum class count as step 04.

Negatives come from gnomAD and have no review status, so they are used
unchanged throughout.

Output
    results/review_status_by_gene.jsonl
    results/review_status_summary.jsonl

Usage
    python 07_review_status.py
    python 07_review_status.py --min-class 5
"""

import argparse
import os
from collections import Counter

import numpy as np

from common import (RESULT_DIR, bootstrap_auroc, by_gene,
                    exact_panel_permutation, fmt_ci, gene_aurocs, read_jsonl,
                    review_stars, write_jsonl)


def star_class(r):
    s = review_stars(r.get("review_status"))
    if s is None:
        return None
    return "2plus" if s >= 2 else "0to1"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-class", type=int, default=5)
    ap.add_argument("--n-boot", type=int, default=2000)
    args = ap.parse_args()

    scored = read_jsonl(os.path.join(RESULT_DIR, "scores.jsonl"))
    included = {r["gene"] for r in read_jsonl(os.path.join(RESULT_DIR, "auroc_by_gene.jsonl"))
                if r["included"]}

    unknown = Counter(r.get("review_status") for r in scored
                      if r["label"] == 1 and star_class(r) is None)
    if unknown:
        print("review statuses not recognised, left out of star splits:")
        for k, v in unknown.items():
            print("  %4d  %s" % (v, k))

    summary = []

    # 1. Share of 2+ star positives per panel.
    print("share of positives with 2+ stars")
    for scope, genes in (("included genes", included), ("all genes", None)):
        for panel in ("oocyte", "control"):
            pos = [r for r in scored if r["label"] == 1 and r["panel"] == panel
                   and (genes is None or r["gene"] in genes)]
            n2 = sum(1 for r in pos if star_class(r) == "2plus")
            share = n2 / len(pos) if pos else None
            summary.append({"analysis": "share_2plus", "scope": scope,
                            "panel": panel, "n_positives": len(pos),
                            "n_2plus": n2, "share_2plus": share})
            print("  %-15s %-8s %4d of %4d  (%.1f%%)"
                  % (scope, panel, n2, len(pos), 100 * share if share is not None else 0))

    # 2. Per-gene AUROC by star class.
    print("\nper-gene AUROC, positives split by review status, all negatives kept")
    print("%-8s %-8s %-6s %5s %5s  %-22s %-22s" %
          ("gene", "panel", "stars", "path", "benign", "ESM-2", "BLOSUM62"))
    gene_rows = []
    for gene, recs in sorted(by_gene(scored).items()):
        neg = [r for r in recs if r["label"] == 0]
        for cls in ("0to1", "2plus"):
            pos = [r for r in recs if r["label"] == 1 and star_class(r) == cls]
            row = {"gene": gene, "panel": recs[0]["panel"], "stars": cls,
                   "n_pathogenic": len(pos), "n_benign": len(neg),
                   "gene_included_in_step04": gene in included}
            if pos and neg:
                sub = pos + neg
                lab = [r["label"] for r in sub]
                e = bootstrap_auroc([r["esm_score"] for r in sub], lab, args.n_boot)
                b = bootstrap_auroc([r["blosum_score"] for r in sub], lab, args.n_boot)
            else:
                e = b = (None, None, None)
            row.update({"esm_auroc": e[0], "esm_lo": e[1], "esm_hi": e[2],
                        "blosum_auroc": b[0], "blosum_lo": b[1], "blosum_hi": b[2]})
            gene_rows.append(row)
            if gene in included and pos:
                print("%-8s %-8s %-6s %5d %5d  %-22s %-22s"
                      % (gene, row["panel"], cls, len(pos), len(neg),
                         fmt_ci(*e), fmt_ci(*b)))
    write_jsonl(os.path.join(RESULT_DIR, "review_status_by_gene.jsonl"), gene_rows)

    # Within-gene contrast on control genes that have enough of both
    # star classes: is 2+ star easier than 0-1 star in the same gene?
    print("\ncontrol genes with >= %d positives in both star classes" % args.min_class)
    print("%-8s %8s %8s %8s %8s" % ("gene", "esm 0-1", "esm 2+", "blo 0-1", "blo 2+"))
    contrast = []
    for gene in sorted(included):
        rows = {r["stars"]: r for r in gene_rows if r["gene"] == gene}
        if rows["0to1"]["panel"] != "control":
            continue
        if (rows["0to1"]["n_pathogenic"] < args.min_class
                or rows["2plus"]["n_pathogenic"] < args.min_class):
            continue
        contrast.append((rows["0to1"]["esm_auroc"], rows["2plus"]["esm_auroc"],
                         rows["0to1"]["blosum_auroc"], rows["2plus"]["blosum_auroc"]))
        print("%-8s %8.3f %8.3f %8.3f %8.3f" % ((gene,) + contrast[-1]))
    if contrast:
        c = np.asarray(contrast)
        print("%-8s %8.3f %8.3f %8.3f %8.3f" % (("mean",) + tuple(c.mean(axis=0))))
        summary.append({"analysis": "within_gene_star_contrast",
                        "n_genes": len(contrast),
                        "esm_0to1_mean": float(c[:, 0].mean()),
                        "esm_2plus_mean": float(c[:, 1].mean()),
                        "blosum_0to1_mean": float(c[:, 2].mean()),
                        "blosum_2plus_mean": float(c[:, 3].mean())})

    # 3. Panel comparison with 0-1 star positives only, both panels.
    low = [r for r in scored if r["label"] == 0 or star_class(r) == "0to1"]
    rows = gene_aurocs(low, ["esm_score", "blosum_score"], min_class=args.min_class)
    inc = [r for r in rows if r["included"]]
    is_oo = [r["panel"] == "oocyte" for r in inc]
    print("\npanel comparison, 0-1 star positives only, min %d per class" % args.min_class)
    for r in inc:
        print("  %-8s %-8s %4d/%-4d  esm %.3f  blosum %.3f"
              % (r["gene"], r["panel"], r["n_pathogenic"], r["n_benign"],
                 r["esm_score_auroc"], r["blosum_score_auroc"]))
    quantities = {
        "esm": [r["esm_score_auroc"] for r in inc],
        "blosum": [r["blosum_score_auroc"] for r in inc],
        "esm_adv": [r["esm_score_auroc"] - r["blosum_score_auroc"] for r in inc],
    }
    if 0 < sum(is_oo) < len(inc):
        print("\n%-10s %8s %8s %8s %10s %8s"
              % ("quantity", "oocyte", "control", "gap", "n>=obs", "p"))
        for name, vals in quantities.items():
            obs, n_ge, n_splits, _ = exact_panel_permutation(vals, is_oo)
            v = np.asarray(vals)
            m = np.asarray(is_oo)
            summary.append({
                "analysis": "panel_0to1_star_only", "quantity": name,
                "oocyte_mean": float(v[m].mean()), "control_mean": float(v[~m].mean()),
                "observed_gap": obs, "n_at_least": n_ge, "n_splits": n_splits,
                "p_one_sided": n_ge / n_splits,
                "n_oocyte_genes": int(m.sum()), "n_control_genes": int((~m).sum()),
                "genes": [r["gene"] for r in inc],
                "counts": {r["gene"]: [r["n_pathogenic"], r["n_benign"]] for r in inc},
            })
            print("%-10s %8.3f %8.3f %8.3f %5d/%-5d %8.4f"
                  % (name, v[m].mean(), v[~m].mean(), obs, n_ge, n_splits, n_ge / n_splits))
    else:
        print("one panel has no included genes, no test")

    write_jsonl(os.path.join(RESULT_DIR, "review_status_summary.jsonl"), summary)
    print("\nwrote results/review_status_by_gene.jsonl and "
          "results/review_status_summary.jsonl")


if __name__ == "__main__":
    main()
