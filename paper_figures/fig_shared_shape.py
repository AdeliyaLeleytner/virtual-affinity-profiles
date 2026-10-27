"""Figure: the panels agree about molecules and disagree about targets.

a  Each molecule's mean score on one broad panel against its mean on the
   other, over the molecules both panels score.
b  How much of the joint 102-column panel's variance the leading directions
   hold, before and after each molecule's own mean is removed.
c  For the nine genes both panels contain, how well the two panels agree
   about that one target, before and after the same subtraction.

Everything is recomputed from the frozen score arrays in data/datasets/.
No value is transcribed.
"""
import csv
import os

import numpy as np
import matplotlib.pyplot as plt

from figstyle import TEXTWIDTH_IN, PALETTE, use_paper_style, figsize, save

HERE = os.path.dirname(os.path.abspath(__file__))
DATASETS = os.path.join(HERE, os.pardir, "data", "datasets")
INK, BLUE, ORANGE, GREY, PALE = (PALETTE["ink"], PALETTE["blue"],
                                 PALETTE["orange"], PALETTE["grey"],
                                 PALETTE["pale"])


def load():
    """The two broad panels on the molecules both of them score completely."""
    a = np.load(os.path.join(DATASETS, "dockstring.npz"), allow_pickle=True)
    b = np.load(os.path.join(DATASETS, "own44.npz"), allow_pickle=True)
    ia = {str(x): i for i, x in enumerate(a["molecule_ids"])}
    ib = {str(x): i for i, x in enumerate(b["molecule_ids"])}
    shared = sorted(set(ia) & set(ib))
    A = a["scores"][[ia[s] for s in shared]]
    B = b["scores"][[ib[s] for s in shared]]
    keep = ~(np.isnan(A).any(1) | np.isnan(B).any(1))
    return (A[keep], B[keep],
            [str(x) for x in a["column_genes"]],
            [str(x) for x in b["column_ids"]])


def fractions(Z):
    """Share of variance per direction, columns standardised first."""
    mu, sd = Z.mean(0), Z.std(0)
    sd[sd == 0] = 1.0
    s = np.linalg.svd((Z - mu) / sd / np.sqrt(len(Z)), compute_uv=False)
    f = s ** 2 / (s ** 2).sum()
    return f


A, B, GENES_A, IDS_B = load()
N = len(A)
RAW = np.hstack([A, B])
CEN = np.hstack([A - A.mean(1, keepdims=True), B - B.mean(1, keepdims=True)])
f_raw, f_cen = fractions(RAW), fractions(CEN)
pr = lambda f: 1.0 / (f ** 2).sum()
r_means = float(np.corrcoef(A.mean(1), B.mean(1))[0, 1])

# the nine genes both panels contain, one column each side
meta = list(csv.DictReader(open(os.path.join(DATASETS, "broad_columns.csv"))))
shared_genes = sorted({m["nominal_gene"] for m in meta
                       if m["shared_gene_between_panels"] == "True"})
col_a = {m["nominal_gene"]: int(m["source_column_index"]) for m in meta
         if m["source"] == "DOCKSTRING" and m["shared_gene_between_panels"] == "True"}
col_b = {m["nominal_gene"]: int(m["source_column_index"]) for m in meta
         if m["source"] == "Docking44" and m["shared_gene_between_panels"] == "True"}

Ac, Bc = A - A.mean(1, keepdims=True), B - B.mean(1, keepdims=True)
agree = {g: (float(np.corrcoef(A[:, col_a[g]], B[:, col_b[g]])[0, 1]),
             float(np.corrcoef(Ac[:, col_a[g]], Bc[:, col_b[g]])[0, 1]))
         for g in shared_genes}
order = sorted(shared_genes, key=lambda g: -agree[g][1])

use_paper_style()
fig = plt.figure(figsize=figsize(TEXTWIDTH_IN, ratio=0.92))
fig.set_layout_engine("none")

# ── a: the panels agree about molecules ─────────────────────────────────
ax = fig.add_axes([0.105, 0.635, 0.300, 0.275])
ax.plot(A.mean(1), B.mean(1), "o", ms=2.0, color=BLUE, alpha=0.45, mew=0)
lo = min(A.mean(1).min(), B.mean(1).min()) - 0.15
hi = max(A.mean(1).max(), B.mean(1).max()) + 0.15
ax.plot([lo, hi], [lo, hi], color=GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_aspect("equal")
ax.set_xlabel("mean score, DOCKSTRING")
ax.set_ylabel("mean score, Docking44")
ax.text(0.04, 0.955, f"r = {r_means:.3f}\nn = {N}", transform=ax.transAxes,
        fontsize=7.5, color=INK, va="top", linespacing=1.5)
ax.set_title("a  Agreement about molecules", loc="left", pad=7,
             fontsize=9, fontweight="bold")

# ── b: and that agreement is one direction ──────────────────────────────
ax = fig.add_axes([0.600, 0.635, 0.360, 0.275])
k = np.arange(1, len(f_raw) + 1)
for f, col, lab, dy in ((f_raw, BLUE, "as scored", 5.5),
                        (f_cen, ORANGE, "molecule mean removed", -6.0)):
    ax.plot(k, 100 * np.cumsum(f), color=col, lw=1.3)
    j = int(len(f) * 0.60)
    ax.text(k[j], 100 * np.cumsum(f)[j] + dy, lab, color=col, fontsize=7.2,
            ha="center", va="bottom" if dy > 0 else "top")
ax.plot([1], [100 * f_raw[0]], "o", color=BLUE, ms=3.4, mew=0, zorder=3)
ax.annotate(f"one direction holds {100 * f_raw[0]:.0f}%", (1, 100 * f_raw[0]),
            textcoords="offset points", xytext=(9, -11), fontsize=7.2,
            color=BLUE, va="center")
ax.set_xlim(0, len(f_raw))
ax.set_ylim(0, 103)
ax.set_xlabel("directions kept")
ax.set_ylabel("variance held (%)")
ax.text(0.97, 0.06, f"effective dimensionality {pr(f_raw):.2f}, "
        f"then {pr(f_cen):.2f}", transform=ax.transAxes, fontsize=7.2,
        color=INK, ha="right")
ax.set_title(f"b  {RAW.shape[1]} columns, one direction", loc="left",
             pad=7, fontsize=9, fontweight="bold")

# ── c: and they disagree about targets ──────────────────────────────────
ax = fig.add_axes([0.115, 0.075, 0.845, 0.370])
y = np.arange(len(order))
for i, g in enumerate(order):
    raw, cen = agree[g]
    ax.plot([cen, raw], [i, i], color=PALE, lw=3.6, solid_capstyle="round",
            zorder=1)
    ax.plot([raw], [i], "o", color=BLUE, ms=4.0, mew=0, zorder=3)
    ax.plot([cen], [i], "s", color=ORANGE, mfc="white", mec=ORANGE, mew=1.1,
            ms=3.8, zorder=3)
    ax.text(cen - 0.022, i, f"{cen:.2f}", fontsize=7, color=ORANGE,
            ha="right", va="center")
    ax.text(raw + 0.022, i, f"{raw:.2f}", fontsize=7, color=BLUE,
            ha="left", va="center")
ax.set_yticks(y)
ax.set_yticklabels(order)
ax.set_ylim(len(order) - 0.4, -1.5)
ax.set_xlim(0, 1.0)
ax.set_xlabel("agreement between the two panels about that one target "
              "(Pearson r)")
for side in ("left", "right", "top"):
    ax.spines[side].set_visible(False)
ax.text(0.015, -1.08, "molecule mean removed", color=ORANGE, fontsize=7.4,
        va="center", ha="left")
ax.text(0.985, -1.08, "as scored", color=BLUE, fontsize=7.4,
        va="center", ha="right")
ax.set_title("c  Agreement about one target at a time", loc="left",
             pad=20, fontsize=9, fontweight="bold")

save(fig, "fig_shared_shape")
print(f"n={N}  r_means={r_means:.4f}  PC1_raw={100 * f_raw[0]:.1f}%  "
      f"PR {pr(f_raw):.3f} -> {pr(f_cen):.3f}  PC1_cen={100 * f_cen[0]:.1f}%")
for g in order:
    print(f"  {g:9s} raw {agree[g][0]:+.3f}   centred {agree[g][1]:+.3f}")
