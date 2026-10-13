#!/usr/bin/env python3
"""
12_figures_extra.py

Figures for the follow-up analyses, in the same style as 05_figures.py
(the style settings and palette are imported from it, not copied).

    fig4_permutation_null.png     exact permutation null distributions for
                                  the ESM-2, BLOSUM62 and advantage gaps,
                                  observed gap marked (step 06)
    fig5_tubb8_position_classes.png  TUBB8 AUROC at paralog-shared versus
                                  TUBB8-specific positions (step 08)
    fig6_auroc_vs_pppl.png        per-gene ESM-2 AUROC and advantage
                                  against pseudo-perplexity (step 09)

A figure whose input file does not exist yet is skipped with a message,
so this can be rerun as each local step is completed.

Usage
    python 12_figures_extra.py
    python 12_figures_extra.py --setting min_af=0.001
"""

import argparse
import importlib
import os

import matplotlib.pyplot as plt
import numpy as np

from common import RESULT_DIR, read_jsonl

fig05 = importlib.import_module("05_figures")      # sets rcParams on import
GRAD, COLOUR, LABEL, save = fig05.GRAD, fig05.COLOUR, fig05.LABEL, fig05.save

QUANTITY_LABEL = {
    "esm": "ESM-2 gap",
    "blosum": "BLOSUM62 gap",
    "esm_adv": "ESM-2 advantage gap",
}


def fig_permutation_null(setting):
    path = os.path.join(RESULT_DIR, "permutation_null.jsonl")
    if not os.path.exists(path):
        print("skipping fig4, run 06_permutation_test.py first")
        return
    rows = {r["quantity"]: r for r in read_jsonl(path) if r["setting"] == setting}
    tests = {r["quantity"]: r for r in read_jsonl(os.path.join(RESULT_DIR, "permutation_test.jsonl"))
             if r["setting"] == setting}
    order = [q for q in ("esm_adv", "esm", "blosum") if q in rows]
    if not order:
        print("skipping fig4, no rows for %s" % setting)
        return

    fig, axes = plt.subplots(1, len(order), figsize=(3.3 * len(order), 3.2),
                             squeeze=False)
    for ax, q in zip(axes[0], order):
        null = np.asarray(rows[q]["null_gaps"])
        obs = rows[q]["observed_gap"]
        t = tests[q]
        bins = np.linspace(null.min(), max(null.max(), obs), 32)
        ax.hist(null, bins=bins, color=GRAD(0.45), alpha=0.75,
                edgecolor="white", linewidth=0.5)
        # Splits at least as extreme as observed, in the darker shade.
        ax.hist(null[null >= obs - 1e-12], bins=bins, color=GRAD(0.18),
                edgecolor="white", linewidth=0.5)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.3)      # headroom for the label
        ax.axvline(obs, color=GRAD(0.18), linewidth=1.6)
        ax.axvline(0.0, color="0.35", linestyle="--", linewidth=1.0)
        ax.text(0.03, 0.97, "observed %.3f\np = %d/%d" % (obs, t["n_at_least"], t["n_splits"]),
                transform=ax.transAxes, ha="left", va="top", fontsize=8, color="0.2",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
        ax.set_title(QUANTITY_LABEL.get(q, q))
        ax.set_xlabel("control minus oocyte, mean AUROC")
    axes[0][0].set_ylabel("label assignments")
    fig.suptitle("Exact permutation null, all %d ways to place %d oocyte labels on %d genes"
                 % (tests[order[0]]["n_splits"], tests[order[0]]["n_oocyte_genes"],
                    tests[order[0]]["n_oocyte_genes"] + tests[order[0]]["n_control_genes"]),
                 fontsize=10)
    fig.tight_layout()
    suffix = "" if setting == "min_af=0.0001" else "_" + setting.replace("=", "").replace(".", "p")
    save(fig, "fig4_permutation_null%s.png" % suffix)


def fig_tubb8_classes():
    path = os.path.join(RESULT_DIR, "tubb8_position_classes.jsonl")
    if not os.path.exists(path):
        print("skipping fig5, run 08_tubb8_paralog_positions.py first")
        return
    rows = [r for r in read_jsonl(path) if r["analysis"] == "primary"
            and r["class"] in ("shared", "specific")]
    test = [r for r in read_jsonl(path) if r["analysis"] == "primary"
            and r["class"] == "shared_minus_specific"][0]

    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    width = 0.36
    for k, (method, shade) in enumerate((("esm", 0.18), ("blosum", 0.72))):
        xs, ys, lo, hi = [], [], [], []
        for i, r in enumerate(rows):
            if r[method + "_auroc"] is None:
                continue
            xs.append(i + (k - 0.5) * width)
            ys.append(r[method + "_auroc"])
            lo.append(r[method + "_auroc"] - (r[method + "_lo"] or r[method + "_auroc"]))
            hi.append((r[method + "_hi"] or r[method + "_auroc"]) - r[method + "_auroc"])
        ax.bar(xs, ys, width=width * 0.92, yerr=[lo, hi], capsize=4,
               color=GRAD(shade), edgecolor="white", linewidth=0.8,
               label="ESM-2" if method == "esm" else "BLOSUM62")
    ax.axhline(0.5, color="0.35", linestyle="--", linewidth=1.0)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(["%s positions\n(%d path / %d benign)"
                        % ("paralog-shared" if r["class"] == "shared" else "TUBB8-specific",
                           r["n_pathogenic"], r["n_benign"]) for r in rows], fontsize=8)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("AUROC within class")
    note = ""
    if test.get("esm_diff") is not None:
        note = "\nESM-2 shared minus specific %+.3f, permutation p = %.3f" % (
            test["esm_diff"], test["esm_p"])
    ax.set_title("TUBB8, by conservation across human beta-tubulins%s" % note, fontsize=9)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    save(fig, "fig5_tubb8_position_classes.png")


def fig_auroc_vs_pppl(model):
    path = os.path.join(RESULT_DIR, "pseudo_perplexity.jsonl")
    if not os.path.exists(path):
        print("skipping fig6, run 09_pseudo_perplexity.py first")
        return
    pppl = {r["gene"]: r for r in read_jsonl(path) if r.get("model") == model}
    cors = {(r["measure"], r["target"], r.get("scope")): r for r in read_jsonl(
        os.path.join(RESULT_DIR, "pseudo_perplexity_correlations.jsonl"))
        if r.get("model") == model and "target" in r}
    rows = [r for r in read_jsonl(os.path.join(RESULT_DIR, "auroc_by_gene.jsonl"))
            if r["included"] and r["gene"] in pppl]
    if not rows:
        print("skipping fig6, no pseudo-perplexity rows for model %s" % model)
        return

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
    for ax, (target, ylabel) in zip(axes, (("esm_auroc", "ESM-2 AUROC"),
                                           ("esm_adv", "ESM-2 minus BLOSUM62 AUROC"))):
        for r in rows:
            x = pppl[r["gene"]]["pseudo_perplexity"]
            y = r["esm_auroc"] if target == "esm_auroc" else r["esm_auroc"] - r["blosum_auroc"]
            ax.plot(x, y, "o", color=COLOUR[r["panel"]], ms=7,
                    markeredgecolor="white", markeredgewidth=0.8)
            ax.annotate(r["gene"], (x, y), xytext=(4, 3), textcoords="offset points",
                        fontsize=7, color="0.25")
        c = cors.get(("pseudo_perplexity", target, "all_included"))
        title = ylabel
        if c:
            title += "\nSpearman rho %+.2f, permutation p = %.3f (n=%d)" % (
                c["spearman_rho"], c["p_two_sided"], c["n_genes"])
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("ESM-2 pseudo-perplexity of the wild-type protein")
        ax.set_ylabel(ylabel)
    handles = [plt.Line2D([], [], color=COLOUR[p], marker="o", linestyle="",
                          label=LABEL[p]) for p in ("oocyte", "control")]
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8)
    save(fig, "fig6_auroc_vs_pppl.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--setting", default="min_af=0.0001",
                    help="which step 06 setting to plot in fig4")
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D",
                    help="which step 09 model run to plot in fig6")
    args = ap.parse_args()
    fig_permutation_null(args.setting)
    fig_tubb8_classes()
    fig_auroc_vs_pppl(args.model)


if __name__ == "__main__":
    main()
