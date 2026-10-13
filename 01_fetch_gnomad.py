#!/usr/bin/env python3
"""
01_fetch_gnomad.py

Fetch missense variants from gnomAD and keep those common enough to be
treated as putatively benign. These are the negatives.

The logic: a missense variant carried by many healthy people at
appreciable frequency is unlikely to cause a severe phenotype. This is
the standard negative set in variant effect predictor benchmarking, and
it is the only way to get negatives for genes like TUBB8 where ClinVar
has none.

Read the caveat in the README before writing anything up. "Common in
gnomAD" is not the same as "proven harmless", and for recessive
conditions a pathogenic allele can sit at appreciable frequency.

Output
    data/gnomad_benign.jsonl
    data/counts_gnomad.jsonl

Usage
    python 01_fetch_gnomad.py
    python 01_fetch_gnomad.py --min-af 1e-3     # stricter, fewer negatives
"""

import argparse
import json
import os
import re
import time
import urllib.error
import urllib.request

from common import (AA3TO1, ALL_GENES, DATA_DIR, PANEL, append_jsonl,
                    write_jsonl)

API = "https://gnomad.broadinstitute.org/api"
PROT = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")

QUERY = """
{
  gene(gene_symbol: "%s", reference_genome: GRCh38) {
    gene_id
    variants(dataset: gnomad_r4) {
      consequence
      hgvsp
      genome { af ac an }
      exome  { af ac an }
    }
  }
}
"""


def query_gene(symbol, retries=3):
    body = json.dumps({"query": QUERY % symbol}).encode("utf-8")
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                API, data=body,
                headers={"Content-Type": "application/json",
                         "User-Agent": "variant-benchmark/1.0"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                payload = json.load(resp)
            if "errors" in payload and payload["errors"]:
                print("  API error for %s: %s"
                      % (symbol, str(payload["errors"])[:120]))
                return []
            gene = payload.get("data", {}).get("gene")
            return gene["variants"] if gene else []
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            wait = 5 * (attempt + 1)
            print("  %s failed (%s), retrying in %ds" % (symbol, type(exc).__name__, wait))
            time.sleep(wait)
    print("  giving up on %s" % symbol)
    return []


def best_af(variant):
    """Highest allele frequency across the genome and exome call sets."""
    best = 0.0
    for key in ("genome", "exome"):
        block = variant.get(key)
        if block and block.get("af"):
            best = max(best, float(block["af"]))
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-af", type=float, default=1e-4,
                    help="allele frequency floor for calling a variant benign")
    ap.add_argument("--sleep", type=float, default=1.0,
                    help="pause between genes, be polite to the API")
    args = ap.parse_args()

    records, seen = [], set()
    counts_path = os.path.join(DATA_DIR, "counts_gnomad.jsonl")
    if os.path.exists(counts_path):
        os.remove(counts_path)

    print("querying gnomAD r4 for %d genes, AF floor %.0e\n" % (len(ALL_GENES), args.min_af))
    print("%-10s %-8s %10s %10s" % ("gene", "panel", "missense", "kept"))

    for gene in ALL_GENES:
        variants = query_gene(gene)
        n_missense, kept = 0, 0

        for v in variants:
            if v.get("consequence") != "missense_variant":
                continue
            hgvsp = v.get("hgvsp") or ""
            m = PROT.match(hgvsp)
            if not m:
                continue
            wt3, pos, mt3 = m.group(1), int(m.group(2)), m.group(3)
            if wt3 not in AA3TO1 or mt3 not in AA3TO1:
                continue
            wt, mt = AA3TO1[wt3], AA3TO1[mt3]
            if wt == mt:
                continue

            n_missense += 1
            af = best_af(v)
            if af < args.min_af:
                continue
            key = (gene, pos, wt, mt)
            if key in seen:
                continue
            seen.add(key)
            kept += 1
            records.append({
                "gene": gene, "panel": PANEL[gene],
                "pos": pos, "wt": wt, "mt": mt,
                "source": "gnomad", "af": af, "label": 0,
            })

        append_jsonl(counts_path, {
            "gene": gene, "panel": PANEL[gene],
            "missense": n_missense, "kept": kept, "min_af": args.min_af,
        })
        print("%-10s %-8s %10d %10d" % (gene, PANEL[gene], n_missense, kept))
        time.sleep(args.sleep)

    out = os.path.join(DATA_DIR, "gnomad_benign.jsonl")
    write_jsonl(out, records)
    print("\nwrote %s  (%d putatively benign variants)" % (out, len(records)))


if __name__ == "__main__":
    main()
