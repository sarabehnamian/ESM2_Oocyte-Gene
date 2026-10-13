#!/usr/bin/env python3
"""
11_mave_controls.py

Control-gene AUROC against labels that do not come from ClinVar.

The ClinVar labels on the control genes are not independent of the
scorers: ClinVar curation under ACMG/AMP uses in-silico predictions and
conservation (PP3/BP4) as supporting evidence, and control genes have
been curated far more heavily than oocyte genes. Deep mutational scans
measure function directly in the lab. If ESM-2 still separates
functional from non-functional variants on BRCA1, TP53, PTEN and MSH2
about as well as it separates ClinVar pathogenic from gnomAD variants,
the control-panel ceiling is not an artefact of curation.

Data. One published score set per gene from MaveDB (https://mavedb.org).
Which score set and which threshold is a choice that must be written
down, so it lives in data/mave_config.json, one entry per gene:

    urn            the MaveDB score-set URN
    score_column   column holding the score, normally "score"
    rule           how labels are made
                   "below"   score <= lof_max is non-functional (1),
                             score >= func_min is functional (0),
                             anything between is left out
                   "above"   the mirror image, for scores where higher
                             means loss of function
                   "class"   a classification column, with class_map
                             mapping its values to 1 / 0
    lof_max, func_min   thresholds for "below" / "above"
    offset         added to MaveDB protein positions to reach UniProt
                   numbering, normally 0
    source         citation and where the threshold comes from
    verified       true once you have checked urn, columns and
                   thresholds against the paper

The config shipped here was written WITHOUT access to MaveDB, so every
entry has verified=false. Use --search and --inspect to confirm or
replace each entry, then set verified=true. Entries with verified=false
are refused unless --allow-unverified is given.

Every MaveDB variant is checked against the UniProt sequence in
data/sequences.jsonl; a gene whose wild-type match rate is below 90
percent is stopped, because that means the coordinates are off (set
offset) or the target is a different isoform.

Scoring. The MAVE variants are mostly not in results/scores.jsonl, so
this step runs ESM-2 itself, with the same log-likelihood ratio, window
rule and BLOSUM62 table as step 03 (imported from it, not copied).

Output
    data/mave/<urn>.csv                   cached score set
    results/mave_variants.jsonl           labelled, scored MAVE variants
    results/mave_auroc.jsonl              per gene, MAVE vs ClinVar AUROC

Usage
    python 11_mave_controls.py --search BRCA1
    python 11_mave_controls.py --inspect urn:mavedb:00000097-0-1
    python 11_mave_controls.py                       # verified entries only
    python 11_mave_controls.py --method masked-marginals --device cuda
"""

import argparse
import csv
import importlib
import io
import json
import os
import re
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

import numpy as np

from common import (AA3TO1, DATA_DIR, RESULT_DIR, bootstrap_auroc, fmt_ci,
                    read_jsonl, write_jsonl)

API = "https://api.mavedb.org/api/v1"
CONFIG = os.path.join(DATA_DIR, "mave_config.json")
CACHE = os.path.join(DATA_DIR, "mave")
GENES = ["BRCA1", "TP53", "PTEN", "MSH2"]

# Single missense change, with or without parentheses: p.Arg175His, p.(Arg175His)
HGVS = re.compile(r"^p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)?$")
SYN = re.compile(r"^p\.\(?(=|[A-Z][a-z]{2}\d+=)\)?$")
NONSENSE = re.compile(r"^p\.\(?[A-Z][a-z]{2}\d+(Ter|\*)\)?$")


def http(url, data=None):
    headers = {"User-Agent": "variant-benchmark/1.0", "Accept": "application/json"}
    if data is not None:
        data = json.dumps(data).encode("utf-8")
        headers["Content-Type"] = "application/json"
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=120) as resp:
                return resp.read().decode("utf-8")
        except Exception as exc:                       # noqa: BLE001
            if attempt == 2:
                raise SystemExit("MaveDB request failed: %s (%s)" % (url, exc))
            time.sleep(5 * (attempt + 1))


def search(gene):
    """List score sets whose target gene matches, with title and size."""
    payload = json.loads(http(API + "/score-sets/search", {"text": gene}))
    hits = payload if isinstance(payload, list) else payload.get("scoreSets", [])
    print("%-28s %8s  %s" % ("urn", "variants", "title"))
    for s in hits:
        targets = [t.get("name", "") for t in s.get("targetGenes", [])]
        if not any(gene.upper() in t.upper() for t in targets):
            continue
        print("%-28s %8s  %s" % (s.get("urn"), s.get("numVariants", "?"),
                                 (s.get("title") or "")[:90]))


def fetch_scores(urn):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, urn.replace(":", "_") + ".csv")
    if not os.path.exists(path):
        text = http(API + "/score-sets/%s/scores" % urllib.parse.quote(urn))
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        meta = http(API + "/score-sets/%s" % urllib.parse.quote(urn))
        with open(path[:-4] + ".meta.json", "w", encoding="utf-8") as fh:
            fh.write(meta)
    with open(path, encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(io.StringIO("".join(lines))))


def to_float(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if np.isnan(v) else v


def inspect(urn, score_column="score"):
    rows = fetch_scores(urn)
    print("%s: %d rows, columns: %s" % (urn, len(rows), ", ".join(rows[0].keys())))
    groups = defaultdict(list)
    for r in rows:
        s = to_float(r.get(score_column))
        h = (r.get("hgvs_pro") or "").strip()
        if s is None:
            continue
        if SYN.match(h):
            groups["synonymous"].append(s)
        elif NONSENSE.match(h):
            groups["nonsense"].append(s)
        elif HGVS.match(h):
            groups["missense"].append(s)
        else:
            groups["other"].append(s)
    for k, v in groups.items():
        v = np.asarray(v)
        print("  %-10s n=%6d  median %8.3f  5%% %8.3f  95%% %8.3f"
              % (k, len(v), np.median(v), np.percentile(v, 5), np.percentile(v, 95)))
    meta = os.path.join(CACHE, urn.replace(":", "_") + ".meta.json")
    if os.path.exists(meta):
        m = json.load(open(meta, encoding="utf-8"))
        print("  title: %s" % m.get("title"))
        for p in m.get("primaryPublicationIdentifiers", []):
            print("  publication: %s" % (p.get("url") or p.get("identifier")))


def label(row, cfg):
    s = row.get(cfg.get("score_column", "score"))
    if cfg["rule"] == "class":
        return cfg["class_map"].get((row.get(cfg["class_column"]) or "").strip())
    v = to_float(s)
    if v is None:
        return None
    if cfg["rule"] == "below":
        if v <= cfg["lof_max"]:
            return 1
        if v >= cfg["func_min"]:
            return 0
        return None
    if cfg["rule"] == "above":
        if v >= cfg["lof_max"]:
            return 1
        if v <= cfg["func_min"]:
            return 0
        return None
    raise SystemExit("unknown rule %s" % cfg["rule"])


def build_gene(gene, cfg, seq):
    rows = fetch_scores(cfg["urn"])
    offset = int(cfg.get("offset", 0))
    out, seen, n_missense, n_match = [], set(), 0, 0
    for r in rows:
        m = HGVS.match((r.get("hgvs_pro") or "").strip())
        if not m or m.group(1) not in AA3TO1 or m.group(3) not in AA3TO1:
            continue
        wt, pos, mt = AA3TO1[m.group(1)], int(m.group(2)) + offset, AA3TO1[m.group(3)]
        if wt == mt:
            continue
        n_missense += 1
        if not (1 <= pos <= len(seq) and seq[pos - 1] == wt):
            continue
        n_match += 1
        lab = label(r, cfg)
        if lab is None or (pos, wt, mt) in seen:
            continue                   # intermediate, unlabelled, or duplicate codon
        seen.add((pos, wt, mt))
        out.append({"gene": gene, "pos": pos, "wt": wt, "mt": mt, "label": lab,
                    "mave_score": to_float(r.get(cfg.get("score_column", "score"))),
                    "urn": cfg["urn"]})
    rate = n_match / n_missense if n_missense else 0.0
    return out, n_missense, rate


def score_esm(records, seqs, args):
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    step03 = importlib.import_module("03_score_esm")

    device = ("cuda" if torch.cuda.is_available() else "cpu") \
        if args.device == "auto" else args.device
    print("\nscoring %d MAVE variants with %s (%s, %s)"
          % (len(records), args.model, args.method, device))
    tok = AutoTokenizer.from_pretrained(args.model)
    model = AutoModelForMaskedLM.from_pretrained(args.model).to(device).eval()

    @torch.no_grad()
    def logprobs(sub, mask_at=None):
        enc = tok(sub, return_tensors="pt").to(device)
        if mask_at is not None:
            enc["input_ids"][0, mask_at + 1] = tok.mask_token_id
        return torch.log_softmax(model(**enc).logits[0], dim=-1)

    by = defaultdict(list)
    for r in records:
        by[r["gene"]].append(r)
    for gene, group in by.items():
        seq, cache, mcache = seqs[gene], {}, {}
        t0 = time.time()
        for r in group:
            sub, off = step03.window_for(seq, r["pos"])
            local = r["pos"] - 1 - off
            if args.method == "wt-marginals":
                if off not in cache:
                    cache[off] = logprobs(sub)
                row = cache[off][local + 1]
            else:
                key = (off, local)       # one masked pass per position
                if key not in mcache:
                    mcache[key] = logprobs(sub, mask_at=local)[local + 1]
                row = mcache[key]
            llr = float(row[tok.convert_tokens_to_ids(r["mt"])]
                        - row[tok.convert_tokens_to_ids(r["wt"])])
            r["esm_llr"], r["esm_score"] = llr, -llr
            r["blosum_score"] = -step03.BLOSUM62.get((r["wt"], r["mt"]), 0)
        print("  %-6s %5d variants  %.0fs" % (gene, len(group), time.time() - t0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--search", metavar="GENE")
    ap.add_argument("--inspect", metavar="URN")
    ap.add_argument("--allow-unverified", action="store_true",
                    help="run config entries not yet checked against the paper")
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--method", default="wt-marginals",
                    choices=["wt-marginals", "masked-marginals"])
    ap.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    ap.add_argument("--n-boot", type=int, default=2000)
    args = ap.parse_args()

    if args.search:
        search(args.search)
        return
    if args.inspect:
        inspect(args.inspect)
        return

    config = json.load(open(CONFIG, encoding="utf-8"))
    seqs = {r["gene"]: r["sequence"] for r in read_jsonl(os.path.join(DATA_DIR, "sequences.jsonl"))}

    records, stats = [], {}
    for gene in GENES:
        cfg = config.get(gene)
        if not cfg or not cfg.get("urn"):
            print("%-6s no score set configured, skipped" % gene)
            continue
        if not cfg.get("verified") and not args.allow_unverified:
            print("%-6s config not verified, skipped (check it, set verified=true, "
                  "or pass --allow-unverified)" % gene)
            continue
        if cfg["rule"] in ("below", "above") and (cfg.get("lof_max") is None
                                                   or cfg.get("func_min") is None):
            print("%-6s thresholds not set in %s, skipped" % (gene, CONFIG))
            continue
        recs, n_missense, rate = build_gene(gene, cfg, seqs[gene])
        c = Counter(r["label"] for r in recs)
        print("%-6s %s  %d missense, wt match %.1f%%, labelled %d non-functional / %d functional"
              % (gene, cfg["urn"], n_missense, 100 * rate, c[1], c[0]))
        if rate < 0.9:
            print("       wild-type match below 90%%, coordinates do not line up; "
                  "fix offset in %s. Gene skipped." % CONFIG)
            continue
        stats[gene] = {"n_missense": n_missense, "wt_match_rate": rate}
        records += recs

    if not records:
        raise SystemExit("nothing to score")
    score_esm(records, seqs, args)
    write_jsonl(os.path.join(RESULT_DIR, "mave_variants.jsonl"), records)

    clinvar = {r["gene"]: r for r in read_jsonl(os.path.join(RESULT_DIR, "auroc_by_gene.jsonl"))}
    out = []
    print("\n%-6s %11s  %-22s %-22s  %-12s %-12s"
          % ("gene", "nonfun/fun", "ESM-2 (MAVE)", "BLOSUM62 (MAVE)", "ESM ClinVar", "BLO ClinVar"))
    for gene in GENES:
        g = [r for r in records if r["gene"] == gene]
        if not g:
            continue
        lab = [r["label"] for r in g]
        e = bootstrap_auroc([r["esm_score"] for r in g], lab, args.n_boot)
        b = bootstrap_auroc([r["blosum_score"] for r in g], lab, args.n_boot)
        cv = clinvar.get(gene, {})
        cfg = config[gene]
        out.append({"gene": gene, "urn": cfg["urn"], "rule": cfg["rule"],
                    "lof_max": cfg.get("lof_max"), "func_min": cfg.get("func_min"),
                    "source": cfg.get("source"), "verified": bool(cfg.get("verified")),
                    "n_nonfunctional": lab.count(1), "n_functional": lab.count(0),
                    "esm_auroc": e[0], "esm_lo": e[1], "esm_hi": e[2],
                    "blosum_auroc": b[0], "blosum_lo": b[1], "blosum_hi": b[2],
                    "clinvar_esm_auroc": cv.get("esm_auroc"),
                    "clinvar_blosum_auroc": cv.get("blosum_auroc"),
                    "model": args.model, "method": args.method, **stats[gene]})
        print("%-6s %5d/%-5d  %-22s %-22s  %-12.3f %-12.3f"
              % (gene, lab.count(1), lab.count(0), fmt_ci(*e), fmt_ci(*b),
                 cv.get("esm_auroc") or float("nan"), cv.get("blosum_auroc") or float("nan")))
    write_jsonl(os.path.join(RESULT_DIR, "mave_auroc.jsonl"), out)
    print("\nwrote results/mave_variants.jsonl and results/mave_auroc.jsonl")


if __name__ == "__main__":
    main()
