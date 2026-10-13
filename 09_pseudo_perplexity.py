#!/usr/bin/env python3
"""
09_pseudo_perplexity.py

A direct test of the "thin evolutionary signal" explanation.

If ESM-2 does badly on oocyte genes because it has learned little about
them, it should also find their wild-type sequences hard to predict.
Two per-gene measures of that:

  pseudo-perplexity   mask each position in turn, take the model's
                      probability of the true residue, and report
                      exp(-mean log p). High means the model cannot
                      reconstruct the sequence from its context.
  wt marginal entropy one unmasked forward pass; entropy (nats) of the
                      predicted distribution at each position, over the
                      20 standard amino acids renormalised, averaged
                      along the protein. High means the model is unsure
                      what belongs at each position. The same entropy at
                      the masked position is also recorded.

A model-free reference is reported next to them: the Shannon entropy
(nats) of the protein's own amino-acid composition. A smoke test with a
randomly initialised model already produced a sizeable rank correlation
between pseudo-perplexity and AUROC, which can only come from sequence
composition and untrained-network biases, not from learned evolution.
So run the step twice: once with the trained model and once with
--random-init (same architecture, random weights, fixed seed). A
pseudo-perplexity correlation means something only if it clearly beats
the random-init one.

All are then correlated, across the included genes, with per-gene
ESM-2 AUROC and with the ESM-minus-BLOSUM advantage. Spearman rho with
a two-sided permutation p-value. Under the explanation, higher
perplexity should go with lower AUROC and lower advantage.

TUBB8 is a beta-tubulin, a deeply conserved family. If ESM-2 predicts
TUBB8 easily and still scores its variants poorly, the explanation does
not hold for TUBB8, and the output says so either way.

Long proteins. ESM-2 takes 1022 residues. For proteins longer than that
each position is scored in a window of 1022 residues centred on it,
clipped at the ends, the same rule step 03 uses for variants. The model
therefore never sees the full protein for BRCA1, BRCA2, ATM, CFTR and
other long genes; this is reported per gene in the output
("windowed": true). The unmasked entropy pass on long proteins uses
windows every 256 residues and reads each position from the window in
which it is most central.

Cost. One masked forward pass per residue, plus a few unmasked passes.
The 14 included genes total 16,752 residues, all 28 genes 24,775. With
the 650M model that is roughly 15-40 minutes on a recent GPU and most of
a day on a CPU. Run it twice (trained, then --random-init).

Output
    results/pseudo_perplexity.jsonl               one row per gene
    results/pseudo_perplexity_correlations.jsonl  the correlation tests

Usage
    python 09_pseudo_perplexity.py                       # included genes
    python 09_pseudo_perplexity.py --all-genes
    python 09_pseudo_perplexity.py --batch-size 16 --device cuda
    python 09_pseudo_perplexity.py --random-init         # untrained control, run this too
    python 09_pseudo_perplexity.py --analyse-only        # reuse saved per-gene values
"""

import argparse
import math
import os
import time

import numpy as np

from common import (DATA_DIR, RESULT_DIR, append_jsonl, read_jsonl,
                    spearman_permutation, write_jsonl)

MAX_LEN = 1022
UNMASKED_STRIDE = 256   # spacing of unmasked windows on long proteins
STANDARD = "ACDEFGHIKLMNPQRSTVWY"
OUT_GENES = os.path.join(RESULT_DIR, "pseudo_perplexity.jsonl")


def window_for(seq, pos0):
    """Same rule as step 03: a MAX_LEN window centred on 0-based pos0.
    Returns (subsequence, offset)."""
    if len(seq) <= MAX_LEN:
        return seq, 0
    half = MAX_LEN // 2
    start = max(0, min(pos0 - half, len(seq) - MAX_LEN))
    return seq[start:start + MAX_LEN], start


def score_gene(seq, tok, model, device, batch_size):
    import torch

    aa_ids = torch.tensor(tok.convert_tokens_to_ids(list(STANDARD)), device=device)

    def entropy(logits_rows):
        p = torch.softmax(logits_rows[:, aa_ids], dim=-1)
        return -(p * torch.log(p.clamp_min(1e-12))).sum(-1)

    # Unmasked passes: wild-type marginal entropy at each position. Long
    # proteins use windows every UNMASKED_STRIDE residues, and each position
    # is read from the window where it sits closest to the centre, so this
    # needs a handful of passes rather than one per residue.
    wt_entropy = np.zeros(len(seq))
    if len(seq) <= MAX_LEN:
        offsets = [0]
    else:
        offsets = sorted(set(list(range(0, len(seq) - MAX_LEN + 1, UNMASKED_STRIDE))
                             + [len(seq) - MAX_LEN]))
    cache = {}
    with torch.no_grad():
        for off in offsets:
            sub = seq[off:off + MAX_LEN]
            enc = tok(sub, return_tensors="pt").to(device)
            cache[off] = entropy(model(**enc).logits[0, 1:len(sub) + 1]).cpu().numpy()
    for i in range(len(seq)):
        off = min(offsets, key=lambda o: abs(i - (o + MAX_LEN // 2)))
        wt_entropy[i] = cache[off][i - off]

    # Masked passes, batched. Windows only differ for long proteins, and
    # all windows have the same length, so batches need no padding.
    log_p = np.zeros(len(seq))
    mask_entropy = np.zeros(len(seq))
    with torch.no_grad():
        for b0 in range(0, len(seq), batch_size):
            idx = list(range(b0, min(b0 + batch_size, len(seq))))
            subs, locals_ = [], []
            for i in idx:
                sub, off = window_for(seq, i)
                subs.append(sub)
                locals_.append(i - off)
            enc = tok(subs, return_tensors="pt", padding=True).to(device)
            ids = enc["input_ids"]
            true_ids = ids[torch.arange(len(idx)), torch.tensor(locals_) + 1].clone()
            ids[torch.arange(len(idx)), torch.tensor(locals_) + 1] = tok.mask_token_id
            logits = model(**enc).logits
            rows = logits[torch.arange(len(idx)), torch.tensor(locals_) + 1]
            lp = torch.log_softmax(rows, dim=-1)
            log_p[idx] = lp[torch.arange(len(idx)), true_ids].cpu().numpy()
            mask_entropy[idx] = entropy(rows).cpu().numpy()

    return {
        "pseudo_perplexity": float(math.exp(-log_p.mean())),
        "mean_log_p_true": float(log_p.mean()),
        "wt_marginal_entropy": float(wt_entropy.mean()),
        "masked_entropy": float(mask_entropy.mean()),
    }


def compute(args, genes):
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    device = ("cuda" if torch.cuda.is_available() else "cpu") \
        if args.device == "auto" else args.device
    print("model  %s\ndevice %s\n" % (args.model, device))
    tok = AutoTokenizer.from_pretrained(args.model)
    if args.random_init:
        from transformers import AutoConfig
        torch.manual_seed(args.seed)
        model = AutoModelForMaskedLM.from_config(AutoConfig.from_pretrained(args.model))
        print("random-init control: trained weights NOT loaded, seed %d" % args.seed)
    else:
        model = AutoModelForMaskedLM.from_pretrained(args.model)
    model = model.to(device).eval()

    seqs = {r["gene"]: r for r in read_jsonl(os.path.join(DATA_DIR, "sequences.jsonl"))}
    done = set()
    if os.path.exists(OUT_GENES) and not args.restart:
        done = {r["gene"] for r in read_jsonl(OUT_GENES) if r.get("model") == args.model_tag}
        if done:
            print("resuming, %d genes already done" % len(done))
    elif os.path.exists(OUT_GENES):
        keep = [r for r in read_jsonl(OUT_GENES) if r.get("model") != args.model_tag]
        write_jsonl(OUT_GENES, keep)

    started = time.time()
    for gi, gene in enumerate(genes, 1):
        if gene in done or gene not in seqs:
            continue
        seq = seqs[gene]["sequence"]
        t0 = time.time()
        rec = score_gene(seq, tok, model, device, args.batch_size)
        rec.update({"gene": gene, "length": len(seq), "windowed": len(seq) > MAX_LEN,
                    "window": MAX_LEN, "model": args.model_tag})
        append_jsonl(OUT_GENES, rec)
        print("[%2d/%2d] %-8s %5d aa  PPPL %6.2f  wt entropy %.3f  (%.0fs, %.0fs total)"
              % (gi, len(genes), gene, len(seq), rec["pseudo_perplexity"],
                 rec["wt_marginal_entropy"], time.time() - t0, time.time() - started))


def composition_entropy(seq):
    counts = np.array([seq.count(a) for a in STANDARD], dtype=float)
    p = counts[counts > 0] / counts.sum()
    return float(-(p * np.log(p)).sum())


def analyse(args):
    per_gene = {r["gene"]: r for r in read_jsonl(OUT_GENES) if r.get("model") == args.model_tag}
    seqs = {r["gene"]: r["sequence"]
            for r in read_jsonl(os.path.join(DATA_DIR, "sequences.jsonl"))}
    for gene, rec in per_gene.items():
        rec["composition_entropy"] = composition_entropy(seqs[gene])
    rows = [r for r in read_jsonl(os.path.join(RESULT_DIR, "auroc_by_gene.jsonl"))
            if r["included"] and r["gene"] in per_gene]
    if len(rows) < 4:
        raise SystemExit("need per-gene values for at least 4 included genes")

    print("\n%-8s %-8s %8s %8s %8s %8s %8s" %
          ("gene", "panel", "PPPL", "wt_ent", "comp_ent", "ESM", "adv"))
    for r in sorted(rows, key=lambda r: per_gene[r["gene"]]["pseudo_perplexity"]):
        g = per_gene[r["gene"]]
        print("%-8s %-8s %8.2f %8.3f %8.3f %8.3f %+8.3f"
              % (r["gene"], r["panel"], g["pseudo_perplexity"], g["wt_marginal_entropy"],
                 g["composition_entropy"], r["esm_auroc"],
                 r["esm_auroc"] - r["blosum_auroc"]))

    targets = {"esm_auroc": [r["esm_auroc"] for r in rows],
               "esm_adv": [r["esm_auroc"] - r["blosum_auroc"] for r in rows]}
    out = []
    print("\nSpearman across %d included genes (two-sided permutation p)" % len(rows))
    for measure in ("pseudo_perplexity", "wt_marginal_entropy", "masked_entropy",
                    "composition_entropy"):
        x = [per_gene[r["gene"]][measure] for r in rows]
        for tname, y in targets.items():
            for scope in ("all_included", "control_only"):
                keep = [i for i, r in enumerate(rows)
                        if scope == "all_included" or r["panel"] == "control"]
                if len(keep) < 4:
                    continue
                rho, p = spearman_permutation([x[i] for i in keep], [y[i] for i in keep],
                                              n_perm=args.n_perm)
                out.append({"measure": measure, "target": tname, "scope": scope,
                            "n_genes": len(keep), "spearman_rho": rho,
                            "p_two_sided": p, "n_perm": args.n_perm,
                            "model": args.model_tag})
                print("  %-20s vs %-9s  %-12s n=%2d  rho %+.3f  p %.4f"
                      % (measure, tname, scope, len(keep), rho, p))

    # Where does TUBB8 sit? Under the explanation it should be among the
    # least predictable genes.
    if "TUBB8" in per_gene:
        order = sorted(per_gene, key=lambda g: per_gene[g]["pseudo_perplexity"])
        rank = order.index("TUBB8") + 1
        print("\nTUBB8 pseudo-perplexity rank: %d of %d genes scored (1 = most predictable)"
              % (rank, len(order)))
        out.append({"measure": "tubb8_pppl_rank", "rank": rank, "n_genes": len(order),
                    "tubb8_pseudo_perplexity": per_gene["TUBB8"]["pseudo_perplexity"],
                    "model": args.model_tag})

    path = os.path.join(RESULT_DIR, "pseudo_perplexity_correlations.jsonl")
    others = []
    if os.path.exists(path):
        others = [r for r in read_jsonl(path) if r.get("model") != args.model_tag]
    write_jsonl(path, others + out)
    print("\nwrote %s (rows for %s replaced)" % (path, args.model_tag))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--all-genes", action="store_true",
                    help="score all 28 genes, not just the 14 included ones")
    ap.add_argument("--analyse-only", action="store_true")
    ap.add_argument("--restart", action="store_true",
                    help="discard saved per-gene values instead of resuming")
    ap.add_argument("--n-perm", type=int, default=100000)
    ap.add_argument("--random-init", action="store_true",
                    help="control run: same architecture, random weights")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    # Keep the two runs apart in the output files.
    tag = args.model + (":random-init:seed%d" % args.seed if args.random_init else "")
    args.model_tag = tag

    if not args.analyse_only:
        rows = read_jsonl(os.path.join(RESULT_DIR, "auroc_by_gene.jsonl"))
        genes = [r["gene"] for r in rows if args.all_genes or r["included"]]
        compute(args, genes)
    analyse(args)


if __name__ == "__main__":
    main()
