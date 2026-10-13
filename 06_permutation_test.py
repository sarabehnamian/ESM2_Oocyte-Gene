#!/usr/bin/env python3
"""
06_permutation_test.py

Exact permutation test for the panel gap. This replaces the bootstrap
interval as the primary panel test; the bootstrap intervals from step 04
stay as descriptive summaries.

Why a permutation test. The bootstrap in step 04 resamples four oocyte
genes, which gives an interval with poor coverage at that size. The
permutation test asks a direct question instead: of all the ways to put
the 4 oocyte labels on the 14 included genes, C(14,4) = 1001, how many
give a gap at least as large as the one observed? Gene AUROCs are held
fixed and only the labels move, so the test treats the gene as the unit,
which is what the paper generalises over.

Three quantities are tested:
    esm         control mean minus oocyte mean of per-gene ESM-2 AUROC
    blosum      the same for BLOSUM62
    esm_adv     the same for (ESM-2 AUROC minus BLOSUM62 AUROC) per gene,
                i.e. does the language model's advantage over the
                substitution matrix shrink on oocyte genes

One-sided p = (#splits with gap >= observed) / #splits. The observed
split is one of the splits, so the smallest possible p is 1/#splits.

Sensitivity to the allele frequency floor. Negatives in the scoring set
carry their gnomAD allele frequency, so a stricter floor can be applied
here by filtering. This is NOT identical to rerunning steps 01-03 with
--min-af 1e-3: step 02 dropped variants that were both ClinVar
pathogenic and in gnomAD above 1e-4. At 1e-3 the ones with AF below 1e-3
would no longer clash and would come back as positives, but they were
never scored. Those variants are listed in the output so the size of the
difference is visible. For an exact rerun, run steps 01-03 with
--min-af 1e-3 into a separate copy and pass its scores with --scores.

Output
    results/permutation_test.jsonl     one row per setting and quantity
    results/permutation_null.jsonl     the full null distribution of gaps
    results/permutation_gene_aurocs.jsonl   per-gene AUROCs used

Usage
    python 06_permutation_test.py
    python 06_permutation_test.py --min-af 1e-4 1e-3
    python 06_permutation_test.py --scores other/scores.jsonl --min-af 1e-3 --tag rerun
"""

import argparse
import os

import numpy as np

from common import (DATA_DIR, RESULT_DIR, exact_panel_permutation,
                    gene_aurocs, read_jsonl, write_jsonl)

BUILD_MIN_AF = 1e-4      # the floor the committed scoring set was built at


def variant_key(r):
    return (r["gene"], r["pos"], r["wt"], r["mt"])


def unscored_positives(min_af):
    """
    Pathogenic variants step 02 dropped because they also appeared in
    gnomAD above the build floor, whose AF is below min_af. A true rerun
    at min_af would keep them as positives. Returns a list of dicts.
    """
    seqs = {r["gene"]: r["sequence"]
            for r in read_jsonl(os.path.join(DATA_DIR, "sequences.jsonl"))}

    def valid(r):
        s = seqs.get(r["gene"])
        return s is not None and 1 <= r["pos"] <= len(s) and s[r["pos"] - 1] == r["wt"]

    path = [r for r in read_jsonl(os.path.join(DATA_DIR, "clinvar_pathogenic.jsonl"))
            if r.get("label") == 1 and valid(r)]
    gn = {variant_key(r): r["af"]
          for r in read_jsonl(os.path.join(DATA_DIR, "gnomad_benign.jsonl")) if valid(r)}
    out = []
    for r in path:
        af = gn.get(variant_key(r))
        if af is not None and af < min_af:
            out.append({"gene": r["gene"], "panel": r["panel"], "pos": r["pos"],
                        "wt": r["wt"], "mt": r["mt"], "gnomad_af": af})
    return out


def apply_min_af(records, min_af):
    """Keep all positives and the negatives at or above min_af."""
    return [r for r in records
            if r["label"] == 1 or (r.get("af") is not None and r["af"] >= min_af)]


def run_setting(records, setting, scorers, baseline, min_class):
    keys = [s + "_score" for s in scorers]
    rows = gene_aurocs(records, keys, min_class=min_class)
    inc = [r for r in rows if r["included"]
           and all(r[k + "_auroc"] is not None for k in keys)]
    for r in rows:
        r["setting"] = setting

    is_oo = [r["panel"] == "oocyte" for r in inc]
    quantities = {s: [r[s + "_score_auroc"] for r in inc] for s in scorers}
    for s in scorers:
        if s != baseline:
            quantities[s + "_adv"] = [r[s + "_score_auroc"] - r[baseline + "_score_auroc"]
                                      for r in inc]

    tests, nulls = [], []
    print("\n== %s ==" % setting)
    print("included genes: %d oocyte, %d control"
          % (sum(is_oo), len(inc) - sum(is_oo)))
    print("  oocyte : %s" % ", ".join(r["gene"] for r in inc if r["panel"] == "oocyte"))
    print("  control: %s" % ", ".join(r["gene"] for r in inc if r["panel"] == "control"))
    if sum(is_oo) == 0 or sum(is_oo) == len(inc):
        print("  one panel is empty, nothing to test")
        return rows, tests, nulls

    print("\n%-12s %8s %8s %8s %10s %8s" %
          ("quantity", "oocyte", "control", "gap", "n>=obs", "p"))
    for name, vals in quantities.items():
        obs, n_ge, n_splits, null = exact_panel_permutation(vals, is_oo)
        v = np.asarray(vals)
        oo_mean = float(v[np.asarray(is_oo)].mean())
        ct_mean = float(v[~np.asarray(is_oo)].mean())
        p = n_ge / n_splits
        tests.append({
            "setting": setting, "quantity": name,
            "oocyte_mean": oo_mean, "control_mean": ct_mean,
            "observed_gap": obs, "n_at_least": n_ge, "n_splits": n_splits,
            "p_one_sided": p, "min_possible_p": 1.0 / n_splits,
            "n_oocyte_genes": int(sum(is_oo)),
            "n_control_genes": int(len(inc) - sum(is_oo)),
            "genes": [r["gene"] for r in inc],
        })
        nulls.append({"setting": setting, "quantity": name,
                      "observed_gap": obs, "null_gaps": [round(x, 6) for x in null]})
        print("%-12s %8.3f %8.3f %8.3f %5d/%-5d %8.4f"
              % (name, oo_mean, ct_mean, obs, n_ge, n_splits, p))
    return rows, tests, nulls


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default=os.path.join(RESULT_DIR, "scores.jsonl"))
    ap.add_argument("--min-af", type=float, nargs="+", default=[1e-4, 1e-3],
                    help="allele frequency floors to test; values above the "
                         "build floor are applied by filtering negatives")
    ap.add_argument("--min-class", type=int, default=5)
    ap.add_argument("--scorers", nargs="+", default=["esm", "blosum"],
                    help="score fields to test, without the _score suffix")
    ap.add_argument("--baseline", default="blosum",
                    help="scorer subtracted to form the advantage quantity")
    ap.add_argument("--tag", default="",
                    help="suffix for output files, e.g. 'am' or 'rerun'")
    args = ap.parse_args()

    scored = read_jsonl(args.scores)
    all_rows, all_tests, all_nulls = [], [], []

    for min_af in args.min_af:
        records = apply_min_af(scored, min_af)
        setting = "min_af=%g" % min_af
        n_neg = sum(1 for r in records if r["label"] == 0)
        print("\n%s: %d variants, %d negatives" % (setting, len(records), n_neg))

        rows, tests, nulls = run_setting(records, setting, args.scorers,
                                         args.baseline, args.min_class)

        if min_af > BUILD_MIN_AF and args.scores == os.path.join(RESULT_DIR, "scores.jsonl"):
            missing = unscored_positives(min_af)
            print("\nnote: %d positives a true rerun at %g would add, never scored:"
                  % (len(missing), min_af))
            for m in missing:
                print("  %-7s %s%d%s  gnomAD AF %.2e" %
                      (m["gene"], m["wt"], m["pos"], m["mt"], m["gnomad_af"]))
            for t in tests:
                t["note"] = ("negatives filtered from the 1e-4 scoring set; "
                             "%d positives a true rerun would add were never "
                             "scored" % len(missing))
                t["unscored_positives"] = missing

        all_rows += rows
        all_tests += tests
        all_nulls += nulls

    suffix = ("_" + args.tag) if args.tag else ""
    write_jsonl(os.path.join(RESULT_DIR, "permutation_test%s.jsonl" % suffix), all_tests)
    write_jsonl(os.path.join(RESULT_DIR, "permutation_null%s.jsonl" % suffix), all_nulls)
    write_jsonl(os.path.join(RESULT_DIR, "permutation_gene_aurocs%s.jsonl" % suffix), all_rows)
    print("\nwrote results/permutation_test%s.jsonl, permutation_null%s.jsonl, "
          "permutation_gene_aurocs%s.jsonl" % (suffix, suffix, suffix))


if __name__ == "__main__":
    main()
