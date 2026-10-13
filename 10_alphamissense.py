#!/usr/bin/env python3
"""
10_alphamissense.py

Add AlphaMissense as a third scorer.

AlphaMissense (Cheng et al. 2023, Science) is built on AlphaFold and
fine-tuned with population frequency labels. It is the obvious
comparison a reviewer will ask for, and it scores every possible human
missense change, so no model needs to be run: the public precomputed
table is enough.

Input. The amino-acid substitution table from the AlphaMissense Zenodo
record (https://zenodo.org/records/8208688):
    AlphaMissense_aa_substitutions.tsv.gz      about 1.2 GB
Columns: uniprot_id, protein_variant (e.g. "R201H"), am_pathogenicity,
am_class. This table is keyed by UniProt accession and protein position,
so it does not depend on the genome build.

Steps
  1. filter the big table to the 28 panel accessions from
     data/sequences.jsonl and write data/alphamissense_panel.tsv.gz
     (small, so it can be committed; later runs read this instead)
  2. merge with data/variants.jsonl on (accession, wt, pos, mt) and
     report the match rate per gene. A low rate means AlphaMissense used
     a different isoform or sequence version for that accession
  3. per-gene AUROC for ESM-2, BLOSUM62 and AlphaMissense, computed on
     the variants AlphaMissense covers so the three are compared on the
     same set, with class-stratified bootstrap intervals
  4. panel means and bootstrap gaps with the step 04 code, and the exact
     permutation tests from step 06, now including AlphaMissense and the
     AlphaMissense-minus-BLOSUM advantage

Caveat to state in the paper. AlphaMissense was trained with human and
primate population variants as its benign class. The negatives here are
gnomAD variants above an allele frequency floor, so some of them were
very likely benign-labelled in its training. This favours AlphaMissense
on every gene with negatives, which is all of them, but not necessarily
equally across panels.

Output
    data/alphamissense_panel.tsv.gz
    results/scores_am.jsonl               scoring set plus am_score
    results/alphamissense_match.jsonl     match rate per gene
    results/auroc_by_gene_am.jsonl
    results/panel_summary_am.jsonl
    results/permutation_test_am.jsonl

Usage
    python 10_alphamissense.py --am-file /path/to/AlphaMissense_aa_substitutions.tsv.gz
    python 10_alphamissense.py         # once data/alphamissense_panel.tsv.gz exists
"""

import argparse
import gzip
import importlib
import os
import re
from collections import Counter

import numpy as np

from common import (DATA_DIR, RESULT_DIR, bootstrap_auroc,
                    exact_panel_permutation, fmt_ci, gene_aurocs, read_jsonl,
                    write_jsonl)

PANEL_FILE = os.path.join(DATA_DIR, "alphamissense_panel.tsv.gz")
VARIANT = re.compile(r"^([A-Z])(\d+)([A-Z])$")

# Step 04's bootstrap code, reused as is rather than copied.
step04 = importlib.import_module("04_evaluate")


def filter_panel(src, accessions):
    """Stream the full table and keep rows for the panel accessions."""
    print("filtering %s to %d accessions" % (src, len(accessions)))
    kept, header = 0, None
    with gzip.open(src, "rt") as fin, gzip.open(PANEL_FILE, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            if header is None:
                header = line
                fout.write(line)
                continue
            if line.split("\t", 1)[0] in accessions:
                fout.write(line)
                kept += 1
    print("wrote %s  (%d rows)" % (PANEL_FILE, kept))


def load_panel(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = {name: i for i, name in enumerate(f)}
                continue
            m = VARIANT.match(f[header["protein_variant"]])
            if not m:
                continue
            key = (f[header["uniprot_id"]], m.group(1), int(m.group(2)), m.group(3))
            out[key] = float(f[header["am_pathogenicity"]])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--am-file", default=None,
                    help="full AlphaMissense_aa_substitutions.tsv.gz; filtered "
                         "to data/alphamissense_panel.tsv.gz on first use")
    ap.add_argument("--min-class", type=int, default=5)
    ap.add_argument("--n-boot", type=int, default=2000)
    args = ap.parse_args()

    seqs = read_jsonl(os.path.join(DATA_DIR, "sequences.jsonl"))
    acc = {r["gene"]: r["accession"] for r in seqs}

    if args.am_file:
        filter_panel(args.am_file, set(acc.values()))
    elif not os.path.exists(PANEL_FILE):
        raise SystemExit("no AlphaMissense data. Download "
                         "AlphaMissense_aa_substitutions.tsv.gz from "
                         "https://zenodo.org/records/8208688 and pass it with --am-file")
    am = load_panel(PANEL_FILE)
    print("loaded %d AlphaMissense scores for the panel" % len(am))

    # Match rate against the scoring set.
    variants = read_jsonl(os.path.join(DATA_DIR, "variants.jsonl"))
    hit, tot = Counter(), Counter()
    for v in variants:
        key = (acc.get(v["gene"]), v["wt"], v["pos"], v["mt"])
        tot[(v["gene"], v["label"])] += 1
        if key in am:
            hit[(v["gene"], v["label"])] += 1

    match_rows = []
    print("\n%-8s %-10s %12s %12s %8s" % ("gene", "accession", "pathogenic", "benign", "rate"))
    for gene in sorted({v["gene"] for v in variants}):
        n = tot[(gene, 1)] + tot[(gene, 0)]
        h = hit[(gene, 1)] + hit[(gene, 0)]
        row = {"gene": gene, "accession": acc.get(gene),
               "matched_pathogenic": hit[(gene, 1)], "total_pathogenic": tot[(gene, 1)],
               "matched_benign": hit[(gene, 0)], "total_benign": tot[(gene, 0)],
               "match_rate": h / n if n else None}
        match_rows.append(row)
        print("%-8s %-10s %5d/%-6d %5d/%-6d %7.1f%%"
              % (gene, acc.get(gene), hit[(gene, 1)], tot[(gene, 1)],
                 hit[(gene, 0)], tot[(gene, 0)], 100 * row["match_rate"]))
    write_jsonl(os.path.join(RESULT_DIR, "alphamissense_match.jsonl"), match_rows)

    # Attach scores. Higher am_pathogenicity already means more pathogenic.
    scored = read_jsonl(os.path.join(RESULT_DIR, "scores.jsonl"))
    for r in scored:
        r["am_score"] = am.get((acc.get(r["gene"]), r["wt"], r["pos"], r["mt"]))
    write_jsonl(os.path.join(RESULT_DIR, "scores_am.jsonl"), scored)

    # Same variant set for all three scorers.
    common_set = [r for r in scored if r["am_score"] is not None]
    keys = ["esm_score", "blosum_score", "am_score"]
    rows = gene_aurocs(common_set, keys, min_class=args.min_class)

    print("\nper-gene AUROC on AlphaMissense-covered variants")
    print("%-8s %-8s %5s %6s  %-22s %-22s %-22s"
          % ("gene", "panel", "path", "benign", "ESM-2", "BLOSUM62", "AlphaMissense"))
    usable = {k: {"oocyte": {}, "control": {}} for k in keys}
    from common import by_gene
    grouped = by_gene(common_set)
    for row in rows:
        recs = grouped[row["gene"]]
        labels = [r["label"] for r in recs]
        for k in keys:
            pt, lo, hi = bootstrap_auroc([r[k] for r in recs], labels, args.n_boot)
            name = k.replace("_score", "")
            row[name + "_auroc"], row[name + "_lo"], row[name + "_hi"] = pt, lo, hi
            if row["included"]:
                usable[k][row["panel"]][row["gene"]] = ([r[k] for r in recs], labels)
        print("%-8s %-8s %5d %6d  %-22s %-22s %-22s%s"
              % (row["gene"], row["panel"], row["n_pathogenic"], row["n_benign"],
                 fmt_ci(row["esm_auroc"], row["esm_lo"], row["esm_hi"]),
                 fmt_ci(row["blosum_auroc"], row["blosum_lo"], row["blosum_hi"]),
                 fmt_ci(row["am_auroc"], row["am_lo"], row["am_hi"]),
                 "" if row["included"] else "   (excluded)"))
    write_jsonl(os.path.join(RESULT_DIR, "auroc_by_gene_am.jsonl"), rows)

    # Panel means and bootstrap gaps, descriptive.
    summary = []
    print("\npanel means (two-level bootstrap, as step 04)")
    for k in keys:
        name = k.replace("_score", "")
        for panel in ("oocyte", "control"):
            pt, lo, hi = step04.panel_bootstrap(usable[k][panel], args.n_boot)
            summary.append({"method": name, "panel": panel,
                            "n_genes": len(usable[k][panel]),
                            "auroc": pt, "lo": lo, "hi": hi})
            print("  %-7s %-8s n_genes=%2d  %s"
                  % (name, panel, len(usable[k][panel]), fmt_ci(pt, lo, hi)))
        pt, lo, hi = step04.gap_bootstrap(usable[k]["oocyte"], usable[k]["control"],
                                          args.n_boot)
        summary.append({"method": name, "panel": "control_minus_oocyte",
                        "auroc": pt, "lo": lo, "hi": hi, "n_genes": None})
        print("  %-7s gap      %s" % (name, fmt_ci(pt, lo, hi)))
    write_jsonl(os.path.join(RESULT_DIR, "panel_summary_am.jsonl"), summary)

    # Exact permutation tests, as step 06.
    inc = [r for r in rows if r["included"]]
    is_oo = [r["panel"] == "oocyte" for r in inc]
    tests = []
    if 0 < sum(is_oo) < len(inc):
        quantities = {
            "esm": [r["esm_auroc"] for r in inc],
            "blosum": [r["blosum_auroc"] for r in inc],
            "am": [r["am_auroc"] for r in inc],
            "esm_adv": [r["esm_auroc"] - r["blosum_auroc"] for r in inc],
            "am_adv": [r["am_auroc"] - r["blosum_auroc"] for r in inc],
            "am_minus_esm": [r["am_auroc"] - r["esm_auroc"] for r in inc],
        }
        print("\nexact permutation test, %d oocyte vs %d control genes"
              % (sum(is_oo), len(inc) - sum(is_oo)))
        print("%-13s %8s %8s %8s %10s %8s"
              % ("quantity", "oocyte", "control", "gap", "n>=obs", "p"))
        m = np.asarray(is_oo)
        for name, vals in quantities.items():
            obs, n_ge, n_splits, _ = exact_panel_permutation(vals, is_oo)
            v = np.asarray(vals)
            tests.append({"setting": "alphamissense_covered", "quantity": name,
                          "oocyte_mean": float(v[m].mean()),
                          "control_mean": float(v[~m].mean()),
                          "observed_gap": obs, "n_at_least": n_ge, "n_splits": n_splits,
                          "p_one_sided": n_ge / n_splits,
                          "genes": [r["gene"] for r in inc]})
            print("%-13s %8.3f %8.3f %8.3f %5d/%-5d %8.4f"
                  % (name, v[m].mean(), v[~m].mean(), obs, n_ge, n_splits, n_ge / n_splits))
    else:
        print("one panel has no included genes, no permutation test")
    write_jsonl(os.path.join(RESULT_DIR, "permutation_test_am.jsonl"), tests)
    print("\nwrote results/*_am.jsonl and results/alphamissense_match.jsonl")


if __name__ == "__main__":
    main()
