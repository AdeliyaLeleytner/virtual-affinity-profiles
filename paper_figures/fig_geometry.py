"""Figure 2. How concentrated a score matrix is, and how much of it transfers.

Panel a is recomputed from the frozen score arrays by geometry.py. Panel b is
the measured recovery on a protein withheld from calibration, read from
data/decisions_paired_uncertainty.csv. Rows share the order of panel a, so the
agreement of the two orderings is read straight off the right-hand column.
"""
import numpy as np
import matplotlib.pyplot as plt

from decisions import HELD_OUT_PROTEIN, trade_offs
from geometry import NAME, spd_table, joint_broad
from figstyle import TEXTWIDTH_IN, PALETTE, use_paper_style, figsize, save

tab, n_mol = spd_table()
gain = trade_offs(HELD_OUT_PROTEIN, "native_score")
J = joint_broad()
ks = [k for k in tab if k in gain]
order = sorted(ks, key=lambda k: -tab[k][0])     # most concentrated first
y = np.arange(len(order))

rho = np.corrcoef(
    np.argsort(np.argsort([tab[k][0] for k in ks])),
    np.argsort(np.argsort([gain[k]["d"] for k in ks])))[0, 1]

use_paper_style()
fig = plt.figure(figsize=figsize(TEXTWIDTH_IN, ratio=0.80))
fig.set_layout_engine("none")
axL = fig.add_axes([0.285, 0.600, 0.690, 0.330])
axR = fig.add_axes([0.285, 0.115, 0.690, 0.330])
COOL, HOT, GREY, INK = (PALETTE["blue"], PALETTE["orange"],
                        PALETTE["grey"], PALETTE["ink"])

# -- a: how concentrated each score matrix is -----------------------------
axL.barh(y, [tab[k][0] for k in order], height=0.55, color=COOL, zorder=2)
for i, k in enumerate(order):
    axL.text(tab[k][0] + 0.09, i, f"{tab[k][0]:.2f}", va="center", ha="left",
             color=INK)
axL.axvline(1.0, color=INK, lw=0.8, zorder=3)
axL.set_xlim(0, 3.9)
axL.set_xticks([0, 1, 2, 3])
axL.set_yticks(y)
axL.set_yticklabels([NAME[k] for k in order])
axL.set_ylim(len(order) - 0.4, -0.8)
axL.set_title("a  Variance sits in one direction", loc="left", pad=14)
axL.set_xlabel("effective dimensionality")
axL.text(0.0, 1.015, f"{n_mol} molecules × 10 proteins; "
                     f"joint 102-column panel {J['pr_raw']:.2f}",
         transform=axL.transAxes, ha="left", va="bottom", color=GREY)

# -- b: what that protein's own score recovers ----------------------------
axR.axvline(0, color=INK, lw=0.8, zorder=1)
for i, k in enumerate(order):
    g = gain[k]
    axR.plot([g["lo"], g["hi"]], [i, i], color=HOT, lw=1.1,
             solid_capstyle="butt", zorder=2)
    axR.plot([g["d"]], [i], "o", color=HOT, ms=3.4, mec="white", mew=0.6,
             zorder=3)
    axR.text(g["hi"] + 1.2, i, f"{g['d']:+.1f}", va="center", ha="left",
             color=INK)
axR.set_xlim(-2, 47)
axR.set_xticks([0, 10, 20, 30, 40])
axR.set_yticks(y)
axR.set_yticklabels([NAME[k] for k in order])
axR.set_ylim(len(order) - 0.4, -0.8)
axR.set_xlabel("trade-offs recovered, points vs the chemistry baseline")
axR.set_title("b  and how concentrated it is tracks what transfers",
              loc="left", pad=14)
axR.text(0.0, 1.015, f"rows ordered by a; Spearman ρ = {rho:.2f}",
         transform=axR.transAxes, ha="left", va="bottom", color=GREY)

save(fig, "fig_geometry")
