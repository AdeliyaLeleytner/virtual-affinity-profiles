"""Figure: what cheap chemistry already knows about a panel.

a  Each molecule's mean score across the joint 102-column panel against a
   prediction of it from six ordinary descriptors, cross-validated by
   chemical group. Colour is the ring count, the single descriptor that
   carries most of the signal.
b  The same model fitted to what is left once each molecule's own mean is
   removed, once per column. The part a panel shares is cheap; the part that
   distinguishes one target from another is not.
c  How each descriptor on its own tracks the shared part.

Everything is recomputed from data/datasets/ and data/molecules.csv.
"""
import csv
import os

import numpy as np
import matplotlib.pyplot as plt

from figstyle import TEXTWIDTH_IN, PALETTE, use_paper_style, figsize, save

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, "data")
DATASETS = os.path.join(DATA, "datasets")
INK, BLUE, ORANGE, GREY, PALE = (PALETTE["ink"], PALETTE["blue"],
                                 PALETTE["orange"], PALETTE["grey"],
                                 PALETTE["pale"])
FOLDS = 5

a = np.load(os.path.join(DATASETS, "dockstring.npz"), allow_pickle=True)
b = np.load(os.path.join(DATASETS, "own44.npz"), allow_pickle=True)
ia = {str(x): i for i, x in enumerate(a["molecule_ids"])}
ib = {str(x): i for i, x in enumerate(b["molecule_ids"])}
shared = sorted(set(ia) & set(ib))
A = a["scores"][[ia[s] for s in shared]]
B = b["scores"][[ib[s] for s in shared]]
Dsc = a["source_descriptor_values"][[ia[s] for s in shared]]
DNAME = [str(x) for x in a["source_descriptor_names"]]
keep = ~(np.isnan(A).any(1) | np.isnan(B).any(1) | np.isnan(Dsc).any(1))
ids = [s for s, k in zip(shared, keep) if k]
A, B, Dsc = A[keep], B[keep], Dsc[keep]

group = {r["standard_molecule_id"]: r["chemical_group"]
         for r in csv.DictReader(open(os.path.join(DATA, "molecules.csv")))}
g = np.array([group.get(i, "?") for i in ids])
X = np.column_stack([Dsc, np.ones(len(Dsc))])


def cross_validated(y, seed=0):
    """Prediction and R^2 with whole chemical groups held out."""
    keys = np.array(sorted(set(g)))
    np.random.default_rng(seed).shuffle(keys)
    fold = {k: i % FOLDS for i, k in enumerate(keys)}
    pred = np.empty_like(y)
    for f in range(FOLDS):
        te = np.array([fold[x] == f for x in g])
        beta, *_ = np.linalg.lstsq(X[~te], y[~te], rcond=None)
        pred[te] = X[te] @ beta
    return pred, 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()


MEAN = np.hstack([A, B]).mean(1)
pred_mean, r2_mean = cross_validated(MEAN)
RES = np.hstack([A - A.mean(1, keepdims=True), B - B.mean(1, keepdims=True)])
r2_res = np.array([cross_validated(RES[:, j])[1] for j in range(RES.shape[1])])
rings = Dsc[:, DNAME.index("RingCount")]
single = np.array([np.corrcoef(MEAN, Dsc[:, j])[0, 1] for j in range(len(DNAME))])

use_paper_style()
fig = plt.figure(figsize=figsize(TEXTWIDTH_IN, ratio=0.84))
fig.set_layout_engine("none")

# ── a: the shared part is cheap ─────────────────────────────────────────
ax = fig.add_axes([0.085, 0.400, 0.345, 0.460])
sc = ax.scatter(pred_mean, MEAN, c=rings, cmap="viridis", s=7, linewidths=0,
                alpha=0.9)
lo = min(pred_mean.min(), MEAN.min()) - 0.15
hi = max(pred_mean.max(), MEAN.max()) + 0.15
ax.plot([lo, hi], [lo, hi], color=GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_aspect("equal")
ax.set_xlabel("predicted from six descriptors")
ax.set_ylabel("mean score across the panel")
ax.text(0.04, 0.955, f"$R^2$ = {r2_mean:.2f}\nn = {len(MEAN)}",
        transform=ax.transAxes, fontsize=7.5, color=INK, va="top",
        linespacing=1.5)
ax.set_title("a  The part the panel shares", loc="left", pad=7, fontsize=9,
             fontweight="bold")

cax = ax.inset_axes([0.46, 0.115, 0.49, 0.036])
cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
cb.outline.set_visible(False)
cb.ax.xaxis.set_label_position("top")
cb.set_label("rings in the molecule", labelpad=3, fontsize=6.8)
cb.ax.tick_params(length=2, pad=1.2, labelsize=6.5)

# ── b: the specific part is not ─────────────────────────────────────────
ax = fig.add_axes([0.605, 0.400, 0.360, 0.460])
rng = np.random.default_rng(1)
ax.scatter(rng.normal(1, 0.055, len(r2_res)), r2_res, s=7, linewidths=0,
           color=ORANGE, alpha=0.75)
ax.plot([0.80, 1.20], [np.median(r2_res)] * 2, color=ORANGE, lw=1.4)
ax.plot([-0.20, 0.20], [r2_mean] * 2, color=BLUE, lw=1.4)
ax.plot([0], [r2_mean], "o", color=BLUE, ms=4.5, mew=0)
ax.axhline(0, color=GREY, lw=0.7, ls=(0, (3, 3)), zorder=0)
ax.text(0.26, r2_mean, f"{r2_mean:.2f}", color=BLUE, fontsize=7.5,
        va="center", ha="left")
ax.text(1.26, np.median(r2_res), f"median {np.median(r2_res):.2f}",
        color=ORANGE, fontsize=7.5, va="center", ha="left")
ax.set_xlim(-0.45, 1.95)
ax.set_ylim(-0.12, 0.92)
ax.set_xticks([0, 1])
ax.set_xticklabels(["what the panel\nshares", "what distinguishes\none column"])
ax.set_ylabel("variance the six descriptors explain ($R^2$)")
for side in ("right", "top"):
    ax.spines[side].set_visible(False)
ax.set_title("b  and the part it does not", loc="left", pad=7, fontsize=9,
             fontweight="bold")

# ── c: which descriptor carries it ──────────────────────────────────────
ax = fig.add_axes([0.230, 0.075, 0.560, 0.190])
o = np.argsort(np.abs(single))
y = np.arange(len(o))
ax.barh(y, single[o], height=0.58, color=[BLUE if v > 0 else ORANGE
                                          for v in single[o]], zorder=2)
for i, j in enumerate(o):
    v = single[j]
    ax.text(v + (0.03 if v > 0 else -0.03), i,
            f"{v:+.2f}".replace("-", "\u2212"), fontsize=7,
            color=INK, va="center", ha="left" if v > 0 else "right")
ax.axvline(0, color=INK, lw=0.7)
ax.set_yticks(y)
ax.set_yticklabels([DNAME[j] for j in o])
ax.set_xlim(-0.62, 0.95)
ax.set_xticks([-0.5, 0.0, 0.5])
ax.set_xlabel("correlation with the part the panel shares")
for side in ("left", "right", "top"):
    ax.spines[side].set_visible(False)
ax.set_title("c  Ring count carries most of it", loc="left", pad=7,
             fontsize=9, fontweight="bold")

save(fig, "fig_chemistry_knows")
print(f"shared R2 {r2_mean:.3f} | residual median {np.median(r2_res):.3f} "
      f"range {r2_res.min():.3f}..{r2_res.max():.3f} | "
      f"above 0.3: {(r2_res > 0.3).sum()}/{len(r2_res)}")
