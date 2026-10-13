#!/usr/bin/env python3
"""
05_figures.py

    fig1_auroc_by_gene.png   per-gene AUROC with 95% intervals, oocyte
                             panel against control panel. The main figure.
    fig2_panel_gap.png       panel means and the gap interval.
    fig3_score_distributions.png  ESM-2 score distributions, pathogenic
                             against population, for the genes with the
                             most data in each panel.

Style: Georgia serif, white background, gradient palette, 200 DPI.

Usage
    python 05_figures.py
"""

import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import FIG_DIR, RESULT_DIR, read_jsonl

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Georgia", "DejaVu Serif", "Times New Roman"],
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 200,
    "figure.dpi": 200,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.frameon": False,
})

GRAD = plt.get_cmap("viridis")
COLOUR = {"oocyte": GRAD(0.72), "control": GRAD(0.18)}
LABEL = {"oocyte": "Oocyte / IVF-failure genes", "control": "Conserved control genes"}


def save(fig, name):
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print("wrote %s" % path)


def fig_auroc_by_gene():
    rows = [r for r in read_jsonl(os.path.join(RESULT_DIR, "auroc_by_gene.jsonl"))
            if r["included"] and r["esm_auroc"] is not None]
    if not rows:
        print("skipping fig1, no genes met the minimum class count")
        return

    rows.sort(key=lambda r: (r["panel"] != "oocyte", -r["esm_auroc"]))
    ys = np.arange(len(rows))[::-1]

    fig, ax = plt.subplots(figsize=(7.2, 0.34 * len(rows) + 1.9))
    for y, r in zip(ys, rows):
        lo = r["esm_lo"] if r["esm_lo"] is not None else r["esm_auroc"]
        hi = r["esm_hi"] if r["esm_hi"] is not None else r["esm_auroc"]
        ax.plot([lo, hi], [y, y], color=COLOUR[r["panel"]], linewidth=2.2,
                solid_capstyle="round", alpha=0.85)
        ax.plot(r["esm_auroc"], y, "o", color=COLOUR[r["panel"]], ms=6,
                markeredgecolor="white", markeredgewidth=0.8)

    ax.axvline(0.5, color="0.35", linestyle="--", linewidth=1.0)
    ax.text(0.5, 1.005, "chance", transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=8, color="0.35")

    ax.set_yticks(ys)
    ax.set_yticklabels(["%s  (%d/%d)" % (r["gene"], r["n_pathogenic"], r["n_benign"])
                        for r in rows], fontsize=8)
    ax.set_xlabel("ESM-2 zero-shot AUROC, pathogenic against population variants")
    ax.set_xlim(0.0, 1.0)

    handles = [plt.Line2D([], [], color=COLOUR[p], marker="o", linestyle="-",
                          linewidth=2.2, label=LABEL[p])
               for p in ("oocyte", "control")]
    ax.legend(handles=handles, loc="lower left", fontsize=9)
    ax.set_title("Per-gene discrimination with 95% bootstrap intervals\n"
                 "(counts shown as pathogenic/benign)", fontsize=10)
    save(fig, "fig1_auroc_by_gene.png")


def fig_panel_gap():
    rows = read_jsonl(os.path.join(RESULT_DIR, "panel_summary.jsonl"))
    panels = {(r["method"], r["panel"]): r for r in rows}

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))

    ax = axes[0]
    xs, means, los, his, cols, names = [], [], [], [], [], []
    for i, (method, panel) in enumerate([("esm", "oocyte"), ("esm", "control"),
                                         ("blosum", "oocyte"), ("blosum", "control")]):
        r = panels.get((method, panel))
        if not r or r["auroc"] is None:
            continue
        xs.append(i)
        means.append(r["auroc"])
        los.append(r["auroc"] - (r["lo"] if r["lo"] is not None else r["auroc"]))
        his.append((r["hi"] if r["hi"] is not None else r["auroc"]) - r["auroc"])
        cols.append(COLOUR[panel])
        names.append("%s\n%s" % ("ESM-2" if method == "esm" else "BLOSUM62",
                                 "oocyte" if panel == "oocyte" else "control"))
    ax.bar(xs, means, yerr=[los, his], capsize=4, color=cols, width=0.62,
           edgecolor="white", linewidth=0.8)
    ax.axhline(0.5, color="0.35", linestyle="--", linewidth=1.0)
    ax.set_xticks(xs)
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("mean AUROC across genes")
    ax.set_title("Panel means")

    ax = axes[1]
    ys, pts, errlo, errhi, ticks = [], [], [], [], []
    for i, method in enumerate(["esm", "blosum"]):
        r = panels.get((method, "control_minus_oocyte"))
        if not r or r["auroc"] is None:
            continue
        ys.append(i)
        pts.append(r["auroc"])
        errlo.append(r["auroc"] - (r["lo"] if r["lo"] is not None else r["auroc"]))
        errhi.append((r["hi"] if r["hi"] is not None else r["auroc"]) - r["auroc"])
        ticks.append("ESM-2" if method == "esm" else "BLOSUM62")
    ax.errorbar(pts, ys, xerr=[errlo, errhi], fmt="o", color=GRAD(0.45),
                capsize=4, ms=7, markeredgecolor="white", markeredgewidth=0.8)
    ax.axvline(0.0, color="0.35", linestyle="--", linewidth=1.0)
    ax.set_yticks(ys)
    ax.set_yticklabels(ticks)
    ax.set_ylim(-0.6, len(ys) - 0.4)
    ax.set_xlabel("control minus oocyte AUROC")
    ax.set_title("Gap, 95% interval\n(crossing zero means no difference)",
                 fontsize=10)

    fig.tight_layout()
    save(fig, "fig2_panel_gap.png")


def fig_score_distributions(n_per_panel=2):
    scored = read_jsonl(os.path.join(RESULT_DIR, "scores.jsonl"))
    rows = {r["gene"]: r for r in read_jsonl(
        os.path.join(RESULT_DIR, "auroc_by_gene.jsonl"))}

    by_gene = defaultdict(list)
    for r in scored:
        by_gene[r["gene"]].append(r)

    picks = []
    for panel in ("oocyte", "control"):
        cands = [(g, v) for g, v in by_gene.items()
                 if v[0]["panel"] == panel and rows.get(g, {}).get("included")]
        cands.sort(key=lambda kv: -min(
            sum(1 for r in kv[1] if r["label"] == 1),
            sum(1 for r in kv[1] if r["label"] == 0)))
        picks.extend(cands[:n_per_panel])

    if not picks:
        print("skipping fig3, no included genes")
        return

    fig, axes = plt.subplots(1, len(picks), figsize=(3.1 * len(picks), 3.2),
                             squeeze=False)
    for ax, (gene, recs) in zip(axes[0], picks):
        path = [r["esm_score"] for r in recs if r["label"] == 1]
        beni = [r["esm_score"] for r in recs if r["label"] == 0]
        bins = np.linspace(min(path + beni), max(path + beni), 18)
        ax.hist(beni, bins=bins, alpha=0.62, color=GRAD(0.18),
                label="population (n=%d)" % len(beni))
        ax.hist(path, bins=bins, alpha=0.62, color=GRAD(0.78),
                label="pathogenic (n=%d)" % len(path))
        a = rows[gene]["esm_auroc"]
        ax.set_title("%s  (%s)\nAUROC %.2f" % (gene, recs[0]["panel"], a))
        ax.set_xlabel("ESM-2 score")
        ax.legend(fontsize=7)
    axes[0][0].set_ylabel("variants")
    fig.tight_layout()
    save(fig, "fig3_score_distributions.png")


if __name__ == "__main__":
    fig_auroc_by_gene()
    fig_panel_gap()
    fig_score_distributions()
    print("\nfigures in %s" % FIG_DIR)
