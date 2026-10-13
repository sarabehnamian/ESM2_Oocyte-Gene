#!/usr/bin/env python3
"""
02_build_dataset.py

Fetch the canonical human protein sequence for each gene, then merge the
pathogenic and benign sets and VALIDATE every variant against the
sequence.

Validation matters. ClinVar and gnomAD annotate against RefSeq
transcripts; UniProt canonical may be a different isoform, so a variant's
stated wild-type residue can disagree with the sequence position. Any
variant whose wild-type residue does not match is DROPPED, and the drop
rate is reported per gene. A high drop rate means the isoform is wrong
and that gene's numbers cannot be trusted.

Output
    data/sequences.jsonl
    data/variants.jsonl        the scoring set
    data/counts_dataset.jsonl

Usage
    python 02_build_dataset.py
"""

import argparse
import json
import os
import time
import urllib.parse
import urllib.request

from common import (ALL_GENES, DATA_DIR, PANEL, append_jsonl, read_jsonl,
                    write_jsonl)

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


def fetch_sequence(gene, retries=3):
    """Reviewed human entry for this gene symbol, longest if several."""
    params = {
        "query": "gene_exact:%s AND organism_id:9606 AND reviewed:true" % gene,
        "format": "json",
        "fields": "accession,id,sequence,protein_name",
        "size": "10",
    }
    url = UNIPROT + "?" + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "variant-benchmark/1.0"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                payload = json.load(resp)
            results = payload.get("results", [])
            if not results:
                return None
            best = max(results, key=lambda r: r["sequence"]["length"])
            return {
                "gene": gene,
                "accession": best["primaryAccession"],
                "sequence": best["sequence"]["value"],
                "length": best["sequence"]["length"],
            }
        except Exception as exc:                       # noqa: BLE001
            wait = 5 * (attempt + 1)
            print("  %s failed (%s), retry in %ds" % (gene, type(exc).__name__, wait))
            time.sleep(wait)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-drop-rate", type=float, default=0.5,
                    help="warn loudly if a gene loses more than this fraction")
    args = ap.parse_args()

    print("fetching canonical sequences from UniProt\n")
    seqs = {}
    seq_records = []
    for gene in ALL_GENES:
        rec = fetch_sequence(gene)
        if rec is None:
            print("%-10s NOT FOUND" % gene)
            continue
        seqs[gene] = rec["sequence"]
        seq_records.append(rec)
        print("%-10s %-10s %5d aa" % (gene, rec["accession"], rec["length"]))
        time.sleep(0.3)

    write_jsonl(os.path.join(DATA_DIR, "sequences.jsonl"), seq_records)

    pathogenic = [r for r in read_jsonl(os.path.join(DATA_DIR, "clinvar_pathogenic.jsonl"))
                  if r.get("label") == 1]
    benign = read_jsonl(os.path.join(DATA_DIR, "gnomad_benign.jsonl"))

    kept, dropped = [], {}
    for rec in pathogenic + benign:
        gene, pos, wt = rec["gene"], rec["pos"], rec["wt"]
        seq = seqs.get(gene)
        d = dropped.setdefault(gene, {"kept": 0, "no_seq": 0, "out_of_range": 0,
                                      "wt_mismatch": 0})
        if seq is None:
            d["no_seq"] += 1
            continue
        if pos < 1 or pos > len(seq):
            d["out_of_range"] += 1
            continue
        if seq[pos - 1] != wt:
            d["wt_mismatch"] += 1
            continue
        d["kept"] += 1
        kept.append(rec)

    # A variant can appear as pathogenic in ClinVar and also be present in
    # gnomAD. Those are genuinely ambiguous, so drop them from both sides
    # rather than letting the label depend on file order.
    from collections import Counter
    key_counts = Counter((r["gene"], r["pos"], r["wt"], r["mt"]) for r in kept)
    conflicted = {k for k, n in key_counts.items() if n > 1}
    if conflicted:
        kept = [r for r in kept
                if (r["gene"], r["pos"], r["wt"], r["mt"]) not in conflicted]
        print("\ndropped %d variants labelled both ways" % len(conflicted))

    out = os.path.join(DATA_DIR, "variants.jsonl")
    write_jsonl(out, kept)

    counts_path = os.path.join(DATA_DIR, "counts_dataset.jsonl")
    if os.path.exists(counts_path):
        os.remove(counts_path)

    print("\n%-10s %-8s %5s %5s %8s %10s" %
          ("gene", "panel", "path", "benign", "dropped", "drop_rate"))
    warned = False
    for gene in ALL_GENES:
        g = [r for r in kept if r["gene"] == gene]
        n_pos = sum(1 for r in g if r["label"] == 1)
        n_neg = sum(1 for r in g if r["label"] == 0)
        d = dropped.get(gene, {"kept": 0, "no_seq": 0, "out_of_range": 0,
                               "wt_mismatch": 0})
        n_drop = d["no_seq"] + d["out_of_range"] + d["wt_mismatch"]
        total = n_drop + d["kept"]
        rate = (n_drop / total) if total else 0.0
        flag = ""
        if rate > args.max_drop_rate and total > 0:
            flag = "  <-- CHECK ISOFORM"
            warned = True
        append_jsonl(counts_path, {
            "gene": gene, "panel": PANEL[gene], "pathogenic": n_pos,
            "benign": n_neg, "dropped": n_drop, "drop_rate": rate,
            "wt_mismatch": d["wt_mismatch"], "out_of_range": d["out_of_range"],
        })
        print("%-10s %-8s %5d %5d %8d %9.2f%s"
              % (gene, PANEL[gene], n_pos, n_neg, n_drop, rate, flag))

    n_eval = sum(1 for gene in ALL_GENES
                 if sum(1 for r in kept if r["gene"] == gene and r["label"] == 1) >= 5
                 and sum(1 for r in kept if r["gene"] == gene and r["label"] == 0) >= 5)
    print("\nwrote %s  (%d variants)" % (out, len(kept)))
    print("genes with at least 5 of each class: %d of %d" % (n_eval, len(ALL_GENES)))
    if warned:
        print("high drop rate on some genes means UniProt canonical is not the "
              "transcript ClinVar used -- exclude those genes or map isoforms")


if __name__ == "__main__":
    main()
