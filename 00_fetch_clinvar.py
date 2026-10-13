#!/usr/bin/env python3
"""
00_fetch_clinvar.py

Download the ClinVar release and extract PATHOGENIC missense variants for
the gene panel. These are the positives.

Why not take the benign set from ClinVar too: TUBB8, the flagship oocyte
maturation arrest gene, has 42 pathogenic missense variants and ZERO
benign ones. Several other panel genes are just as lopsided. So the
negatives come from population frequency instead, in step 01. That is
standard practice in the variant effect predictor benchmarking
literature.

Downloads about 440 MB. Cached, so it only happens once.

Output
    data/clinvar_pathogenic.jsonl
    data/counts_clinvar.jsonl

Usage
    python 00_fetch_clinvar.py
    python 00_fetch_clinvar.py --keep-vus     # also record VUS, not scored
"""

import argparse
import collections
import gzip
import os
import re
import urllib.request

from common import AA3TO1, ALL_GENES, DATA_DIR, PANEL, append_jsonl, write_jsonl

URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
LOCAL = os.path.join(DATA_DIR, "variant_summary.txt.gz")

# Matches the protein change inside the ClinVar Name field, e.g.
# NM_177987.3(TUBB8):c.602G>A (p.Arg201His)
PROT = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")


def download():
    if os.path.exists(LOCAL) and os.path.getsize(LOCAL) > 10_000_000:
        print("using cached %s" % LOCAL)
        return
    print("downloading ClinVar, about 440 MB, this takes a few minutes")
    urllib.request.urlretrieve(URL, LOCAL)
    print("saved %s" % LOCAL)


def classify(sig):
    s = sig.lower()
    if "conflicting" in s:
        return "conflicting"
    if s.startswith("pathogenic") or s.startswith("likely pathogenic"):
        return "pathogenic"
    if s.startswith("benign") or s.startswith("likely benign"):
        return "benign"
    if "uncertain" in s:
        return "vus"
    return "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep-vus", action="store_true")
    args = ap.parse_args()

    download()
    wanted = set(ALL_GENES)
    seen = set()
    records = []
    counts = collections.defaultdict(collections.Counter)

    with gzip.open(LOCAL, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            gene = f[ix["GeneSymbol"]]
            if gene not in wanted:
                continue
            if f[ix["Assembly"]] != "GRCh38":
                continue
            if f[ix["Type"]] != "single nucleotide variant":
                continue

            m = PROT.search(f[ix["Name"]])
            if not m:
                continue
            wt3, pos, mt3 = m.group(1), int(m.group(2)), m.group(3)
            if wt3 not in AA3TO1 or mt3 not in AA3TO1:
                continue                      # nonsense, stop-loss, selenocysteine
            wt, mt = AA3TO1[wt3], AA3TO1[mt3]
            if wt == mt:
                continue                      # synonymous at protein level

            label = classify(f[ix["ClinicalSignificance"]])
            counts[gene][label] += 1

            if label == "pathogenic" or (args.keep_vus and label == "vus"):
                key = (gene, pos, wt, mt, label)
                if key in seen:
                    continue                  # same change on several transcripts
                seen.add(key)
                records.append({
                    "gene": gene,
                    "panel": PANEL[gene],
                    "pos": pos,
                    "wt": wt,
                    "mt": mt,
                    "source": "clinvar",
                    "clinvar_class": label,
                    "label": 1 if label == "pathogenic" else None,
                    "review_status": f[ix["ReviewStatus"]],
                    "variation_id": f[ix["VariationID"]],
                })

    out = os.path.join(DATA_DIR, "clinvar_pathogenic.jsonl")
    write_jsonl(out, records)

    counts_path = os.path.join(DATA_DIR, "counts_clinvar.jsonl")
    if os.path.exists(counts_path):
        os.remove(counts_path)

    print("\n%-10s %-8s %10s %8s %6s" % ("gene", "panel", "pathogenic", "benign", "vus"))
    for gene in ALL_GENES:
        k = counts[gene]
        append_jsonl(counts_path, {
            "gene": gene, "panel": PANEL[gene],
            "pathogenic": k["pathogenic"], "benign": k["benign"], "vus": k["vus"],
        })
        print("%-10s %-8s %10d %8d %6d"
              % (gene, PANEL[gene], k["pathogenic"], k["benign"], k["vus"]))

    print("\nwrote %s  (%d pathogenic variants)" % (out, sum(
        1 for r in records if r["clinvar_class"] == "pathogenic")))
    print("note the benign columns -- that is why step 01 exists")


if __name__ == "__main__":
    main()
