#!/usr/bin/env python3
"""
08_tubb8_paralog_positions.py

Does ESM-2 fail on TUBB8 mainly at the positions where TUBB8 differs
from the other human beta-tubulins?

TUBB8 is primate-specific, but it belongs to the beta-tubulin family,
one of the most conserved protein families there is. Most of its
residues are identical across the nine human beta-tubulin paralogs, and
ESM-2 has seen thousands of tubulins in training. If the model's weak
TUBB8 AUROC comes from positions that are TUBB8-specific (not shared
with the paralogs), that points at the lineage-specific part of the
protein. If it is just as weak at shared positions, the "thin
evolutionary signal" explanation does not fit TUBB8.

Steps
  1. fetch the nine human beta-tubulins from UniProt (reviewed entries)
  2. align each paralog to TUBB8; for every TUBB8 position record the
     residue each paralog carries there (a gap counts as different).
     By default this is a pairwise global alignment of each paralog to
     TUBB8 (Biopython, BLOSUM62, gap open -10, extend -0.5), which gives
     the TUBB8-anchored column directly. With --aligner mafft the nine
     sequences go through MAFFT instead and the TUBB8 row is read off.
  3. label each TUBB8 position "shared" if TUBB8's residue is identical
     in at least --shared-frac of the other eight paralogs, otherwise
     "specific". Default 1.0, i.e. identical in all eight; 0.75 (6 of 8)
     is also reported.
  4. ESM-2 and BLOSUM62 AUROC within each class, class-stratified
     bootstrap intervals, and a permutation test of the AUROC difference
     between classes (class labels shuffled within pathogenic and within
     benign variants, so each class keeps its counts).

TUBB8 coordinates are taken from data/sequences.jsonl, the same sequence
the variants were validated against. The fetched TUBB8 must match it.

Output
    data/tubulin_paralogs.jsonl
    results/tubb8_positions.jsonl          one row per TUBB8 residue
    results/tubb8_position_classes.jsonl   AUROC per class and the tests

Usage
    python 08_tubb8_paralog_positions.py                 # fetch if needed, then analyse
    python 08_tubb8_paralog_positions.py --fetch-only    # just download sequences
    python 08_tubb8_paralog_positions.py --aligner mafft
"""

import argparse
import json
import os
import shutil
import subprocess
import tempfile
import time
import urllib.parse
import urllib.request

import numpy as np

from common import (DATA_DIR, RESULT_DIR, auroc, bootstrap_auroc, fmt_ci,
                    read_jsonl, write_jsonl)

PARALOGS = ["TUBB", "TUBB1", "TUBB2A", "TUBB2B", "TUBB3",
            "TUBB4A", "TUBB4B", "TUBB6", "TUBB8"]
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
PARALOG_FILE = os.path.join(DATA_DIR, "tubulin_paralogs.jsonl")


def fetch_one(gene, retries=3):
    """Reviewed human entry whose primary gene name is exactly this symbol."""
    params = {
        "query": "gene_exact:%s AND organism_id:9606 AND reviewed:true" % gene,
        "format": "json",
        "fields": "accession,gene_primary,sequence",
        "size": "10",
    }
    url = UNIPROT + "?" + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "variant-benchmark/1.0"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                payload = json.load(resp)
            hits = []
            for r in payload.get("results", []):
                names = [g.get("geneName", {}).get("value") for g in r.get("genes", [])]
                if gene in names:
                    hits.append(r)
            if not hits:
                return None
            best = max(hits, key=lambda r: r["sequence"]["length"])
            return {"gene": gene, "accession": best["primaryAccession"],
                    "sequence": best["sequence"]["value"],
                    "length": best["sequence"]["length"]}
        except Exception as exc:                       # noqa: BLE001
            wait = 5 * (attempt + 1)
            print("  %s failed (%s), retry in %ds" % (gene, type(exc).__name__, wait))
            time.sleep(wait)
    return None


def fetch_all():
    print("fetching human beta-tubulins from UniProt")
    out = []
    for gene in PARALOGS:
        rec = fetch_one(gene)
        if rec is None:
            raise SystemExit("could not fetch %s from UniProt. If rest.uniprot.org "
                             "is blocked here, run this step on a machine that "
                             "can reach it." % gene)
        out.append(rec)
        print("  %-7s %-8s %4d aa" % (gene, rec["accession"], rec["length"]))
        time.sleep(0.3)
    write_jsonl(PARALOG_FILE, out)
    print("wrote %s" % PARALOG_FILE)
    return out


def columns_pairwise(ref, others):
    """
    For each other sequence, a list the length of ref giving the residue
    aligned to each ref position, or '-' for a gap.
    """
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    cols = {}
    for name, seq in others.items():
        aln = aligner.align(ref, seq)[0]
        mapped = ["-"] * len(ref)
        # aligned gives matching blocks as ((ref_start, ref_end), (seq_start, seq_end))
        for (a0, a1), (b0, b1) in zip(*aln.aligned):
            for k in range(a1 - a0):
                mapped[a0 + k] = seq[b0 + k]
        cols[name] = mapped
    return cols


def columns_mafft(ref_name, seqs):
    """Run MAFFT on all sequences and read columns off the ref row."""
    exe = shutil.which("mafft")
    if exe is None:
        raise SystemExit("mafft not found on PATH; use the default --aligner pairwise")
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, "in.fa")
        with open(fa, "w") as fh:
            for name, s in seqs.items():
                fh.write(">%s\n%s\n" % (name, s))
        res = subprocess.run([exe, "--auto", "--quiet", fa],
                             capture_output=True, text=True, check=True)
    aligned, name = {}, None
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            aligned[name] = ""
        elif name:
            aligned[name] += line.strip()
    ref_row = aligned[ref_name]
    cols = {n: [] for n in seqs if n != ref_name}
    for i, c in enumerate(ref_row):
        if c == "-":
            continue
        for n in cols:
            cols[n].append(aligned[n][i].upper())
    return cols


def diff_permutation(esm, blo, labels, groups, n_perm, seed):
    """
    Permutation tests for class differences (shared minus specific):
    ESM AUROC, BLOSUM AUROC, and the ESM-minus-BLOSUM advantage.
    Class labels are shuffled within pathogenic and within benign
    variants. Two-sided p with the +1 correction.
    """
    esm, blo = np.asarray(esm, float), np.asarray(blo, float)
    labels, groups = np.asarray(labels, int), np.asarray(groups, bool)

    def stats(g):
        out = []
        for s in (esm, blo):
            a, b = auroc(s[g], labels[g]), auroc(s[~g], labels[~g])
            out.append(None if a is None or b is None else a - b)
        out.append(None if None in out else out[0] - out[1])
        return out

    obs = stats(groups)
    if None in obs:
        return obs, [None, None, None]
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(labels == 1), np.flatnonzero(labels == 0)
    hits, done = np.zeros(3), 0
    for _ in range(n_perm):
        g = groups.copy()
        g[pos] = rng.permutation(groups[pos])
        g[neg] = rng.permutation(groups[neg])
        s = stats(g)
        if None in s:
            continue
        done += 1
        hits += np.abs(s) >= np.abs(obs) - 1e-12
    return obs, list((hits + 1) / (done + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch-only", action="store_true")
    ap.add_argument("--refetch", action="store_true")
    ap.add_argument("--aligner", default="pairwise", choices=["pairwise", "mafft"])
    ap.add_argument("--shared-frac", type=float, nargs="+", default=[1.0, 0.75],
                    help="fraction of the other paralogs that must carry TUBB8's "
                         "residue for a position to count as shared; the first "
                         "value is the primary analysis")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--n-perm", type=int, default=10000)
    args = ap.parse_args()

    if args.refetch or not os.path.exists(PARALOG_FILE):
        paralogs = fetch_all()
    else:
        paralogs = read_jsonl(PARALOG_FILE)
        print("using cached %s" % PARALOG_FILE)
    if args.fetch_only:
        return

    seqs = {r["gene"]: r["sequence"] for r in paralogs}
    ref = {r["gene"]: r for r in read_jsonl(os.path.join(DATA_DIR, "sequences.jsonl"))}["TUBB8"]
    if seqs.get("TUBB8") != ref["sequence"]:
        raise SystemExit("fetched TUBB8 differs from data/sequences.jsonl; "
                         "variant coordinates would not line up")
    tubb8 = ref["sequence"]
    others = {g: s for g, s in seqs.items() if g != "TUBB8"}

    if args.aligner == "pairwise":
        cols = columns_pairwise(tubb8, others)
    else:
        cols = columns_mafft("TUBB8", seqs)

    print("\nidentity to TUBB8 over TUBB8 positions")
    for g in sorted(others):
        ident = sum(1 for i, c in enumerate(cols[g]) if c == tubb8[i]) / len(tubb8)
        print("  %-7s %.1f%%" % (g, 100 * ident))

    # Per-position table.
    pos_rows = []
    for i, aa in enumerate(tubb8):
        residues = {g: cols[g][i] for g in sorted(others)}
        n_same = sum(1 for c in residues.values() if c == aa)
        pos_rows.append({"pos": i + 1, "tubb8": aa, "paralog_residues": residues,
                         "n_identical": n_same, "n_paralogs": len(others),
                         "frac_identical": n_same / len(others)})

    scored = [r for r in read_jsonl(os.path.join(RESULT_DIR, "scores.jsonl"))
              if r["gene"] == "TUBB8"]
    frac_at = {r["pos"]: r["frac_identical"] for r in pos_rows}

    summary = []
    for k, thr in enumerate(args.shared_frac):
        tag = "primary" if k == 0 else "sensitivity"
        for r in pos_rows:
            r["class_%g" % thr] = "shared" if r["frac_identical"] >= thr - 1e-9 else "specific"
        n_shared = sum(1 for r in pos_rows if r["frac_identical"] >= thr - 1e-9)
        print("\nthreshold %.2f (%s): %d shared, %d specific positions"
              % (thr, tag, n_shared, len(tubb8) - n_shared))

        groups = [frac_at[r["pos"]] >= thr - 1e-9 for r in scored]
        labels = [r["label"] for r in scored]
        esm = [r["esm_score"] for r in scored]
        blo = [r["blosum_score"] for r in scored]

        print("%-9s %5s %6s  %-22s %-22s" % ("class", "path", "benign", "ESM-2", "BLOSUM62"))
        for cls, flag in (("shared", True), ("specific", False)):
            idx = [i for i, g in enumerate(groups) if g == flag]
            lab = [labels[i] for i in idx]
            e = bootstrap_auroc([esm[i] for i in idx], lab, args.n_boot)
            b = bootstrap_auroc([blo[i] for i in idx], lab, args.n_boot)
            summary.append({"threshold": thr, "analysis": tag, "class": cls,
                            "n_positions": n_shared if flag else len(tubb8) - n_shared,
                            "n_pathogenic": lab.count(1), "n_benign": lab.count(0),
                            "esm_auroc": e[0], "esm_lo": e[1], "esm_hi": e[2],
                            "blosum_auroc": b[0], "blosum_lo": b[1], "blosum_hi": b[2]})
            print("%-9s %5d %6d  %-22s %-22s"
                  % (cls, lab.count(1), lab.count(0), fmt_ci(*e), fmt_ci(*b)))

        obs, p = diff_permutation(esm, blo, labels, groups, args.n_perm, seed=k)
        summary.append({"threshold": thr, "analysis": tag, "class": "shared_minus_specific",
                        "esm_diff": obs[0], "esm_p": p[0],
                        "blosum_diff": obs[1], "blosum_p": p[1],
                        "adv_diff": obs[2], "adv_p": p[2], "n_perm": args.n_perm})
        if obs[0] is None:
            print("a class lacks pathogenic or benign variants, no difference test")
        else:
            print("shared minus specific: ESM %+.3f (p=%.4f)  BLOSUM %+.3f (p=%.4f)  "
                  "advantage %+.3f (p=%.4f)"
                  % (obs[0], p[0], obs[1], p[1], obs[2], p[2]))

    # Pathogenic and benign variants are unlikely to fall evenly across
    # classes; print the split so the AUROCs can be read in context.
    write_jsonl(os.path.join(RESULT_DIR, "tubb8_positions.jsonl"), pos_rows)
    write_jsonl(os.path.join(RESULT_DIR, "tubb8_position_classes.jsonl"), summary)
    print("\nwrote results/tubb8_positions.jsonl and results/tubb8_position_classes.jsonl")


if __name__ == "__main__":
    main()
