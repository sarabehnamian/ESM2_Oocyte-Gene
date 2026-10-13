#!/usr/bin/env python3
"""
03_score_esm.py

Score every variant with ESM-2, zero shot, plus a cheap BLOSUM62 baseline
so the language model has something to beat.

Scoring. The score is the log-likelihood ratio at the mutated position:

    score = log p(mutant | context) - log p(wild type | context)

More negative means the model finds the substitution less plausible, so
scores are NEGATED before evaluation and a higher final score means more
likely pathogenic.

Two methods:
  wt-marginals    one forward pass per sequence window. Fast. Default.
  masked-marginals  mask each variant position and run a pass per
                  position. More accurate, much slower. Use on a GPU.

Long proteins. ESM-2 handles about 1022 residues. Longer sequences are
cut to a window centred on the variant, and the window offset is tracked
so positions stay correct.

Output
    results/scores.jsonl

Usage
    python 03_score_esm.py
    python 03_score_esm.py --method masked-marginals     # slow, better
    python 03_score_esm.py --model facebook/esm2_t12_35M_UR50D   # quick test
"""

import argparse
import os
import time
from collections import defaultdict

import torch

from common import DATA_DIR, RESULT_DIR, read_jsonl, write_jsonl

MAX_LEN = 1022          # ESM-2 budget once the two special tokens are added

# BLOSUM62, as a substitution-plausibility baseline.
BLOSUM62 = {}
_AAS = "ARNDCQEGHILKMFPSTWYV"
_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
for _i, _line in enumerate(_ROWS.strip().split("\n")):
    for _j, _v in enumerate(_line.split()):
        BLOSUM62[(_AAS[_i], _AAS[_j])] = int(_v)


def window_for(seq, pos):
    """Crop long sequences to a window centred on the variant position.
    Returns (subsequence, offset) where offset is the index of the first
    kept residue, so local_pos = pos - offset."""
    if len(seq) <= MAX_LEN:
        return seq, 0
    half = MAX_LEN // 2
    start = max(0, min(pos - 1 - half, len(seq) - MAX_LEN))
    return seq[start:start + MAX_LEN], start


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--method", default="wt-marginals",
                    choices=["wt-marginals", "masked-marginals"])
    ap.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    args = ap.parse_args()

    from transformers import AutoTokenizer, AutoModelForMaskedLM

    device = ("cuda" if torch.cuda.is_available() else "cpu") \
        if args.device == "auto" else args.device
    print("model  %s" % args.model)
    print("method %s" % args.method)
    print("device %s\n" % device)

    tok = AutoTokenizer.from_pretrained(args.model)
    model = AutoModelForMaskedLM.from_pretrained(args.model).to(device).eval()

    variants = read_jsonl(os.path.join(DATA_DIR, "variants.jsonl"))
    seqs = {r["gene"]: r["sequence"]
            for r in read_jsonl(os.path.join(DATA_DIR, "sequences.jsonl"))}

    by_gene = defaultdict(list)
    for v in variants:
        by_gene[v["gene"]].append(v)

    @torch.no_grad()
    def logprobs_for(subseq, mask_at=None):
        """Log-softmax over the vocabulary at every position of subseq.
        mask_at is a 0-based index into subseq to replace with <mask>."""
        enc = tok(subseq, return_tensors="pt").to(device)
        if mask_at is not None:
            enc["input_ids"][0, mask_at + 1] = tok.mask_token_id   # +1 for <cls>
        out = model(**enc).logits[0]
        return torch.log_softmax(out, dim=-1)

    scored = []
    started = time.time()

    for gi, (gene, group) in enumerate(sorted(by_gene.items()), 1):
        seq = seqs.get(gene)
        if seq is None:
            continue

        cache = {}
        for v in group:
            sub, offset = window_for(seq, v["pos"])
            local = v["pos"] - 1 - offset            # 0-based within sub

            if args.method == "wt-marginals":
                if offset not in cache:
                    cache[offset] = logprobs_for(sub)
                lp = cache[offset]
            else:
                lp = logprobs_for(sub, mask_at=local)

            row = lp[local + 1]                      # +1 for <cls>
            wt_id = tok.convert_tokens_to_ids(v["wt"])
            mt_id = tok.convert_tokens_to_ids(v["mt"])
            llr = float(row[mt_id] - row[wt_id])

            rec = dict(v)
            rec["esm_llr"] = llr
            rec["esm_score"] = -llr                  # higher = more pathogenic
            rec["blosum62"] = BLOSUM62.get((v["wt"], v["mt"]), 0)
            rec["blosum_score"] = -BLOSUM62.get((v["wt"], v["mt"]), 0)
            scored.append(rec)

        print("[%2d/%2d] %-10s %4d variants  (%.0fs elapsed)"
              % (gi, len(by_gene), gene, len(group), time.time() - started))

    out = os.path.join(RESULT_DIR, "scores.jsonl")
    write_jsonl(out, scored)
    print("\nwrote %s  (%d scored variants, %.0fs total)"
          % (out, len(scored), time.time() - started))


if __name__ == "__main__":
    main()
