"""Graphical abstract: the reversal, in one picture.

Bars are the median across the eight scores of the fraction of measured
trade-offs each input recovers, read from data/decisions_paired_uncertainty.csv
(exact four-cell subset, equal weight per actual assay pair). The white rule
inside each bar spans the eight scores. Nothing is synthesised.
"""
import statistics as st

import numpy as np
import matplotlib.pyplot as plt

from decisions import HELD_OUT_PROTEIN, NEW_MOLECULE, trade_offs
from figstyle import TEXTWIDTH_IN, PALETTE, use_paper_style, figsize, save


def inputs(regime):
    nat = trade_offs(regime, "native_score")
    prof = trade_offs(regime, "full_profile")
    return [("chemistry", [d["base"] for d in nat.values()]),
            ("the whole\npanel", [d["est"] for d in prof.values()]),
            ("that protein's\nown score", [d["est"] for d in nat.values()])]


PANELS = [("New protein,\nnever calibrated", inputs(HELD_OUT_PROTEIN)),
          ("New molecule,\nprotein calibrated", inputs(NEW_MOLECULE))]

use_paper_style(base_pt=9)
fig = plt.figure(figsize=figsize(TEXTWIDTH_IN, ratio=0.58))
fig.set_layout_engine("none")
axes = [fig.add_axes([0.065, 0.330, 0.400, 0.455]),
        fig.add_axes([0.565, 0.330, 0.400, 0.455])]
COOL, HOT, GREY, INK = (PALETTE["blue"], PALETTE["orange"],
                        PALETTE["grey"], PALETTE["ink"])

for ax, (title, rows) in zip(axes, PANELS):
    x = np.arange(len(rows))
    med = [st.median(v) for _, v in rows]
    ax.bar(x, med, width=0.62, color=[GREY, HOT, COOL], zorder=2)
    for i, (_, v) in enumerate(rows):
        ax.plot([i, i], [min(v), max(v)], color="white", lw=0.9, zorder=3)
        ax.text(i, med[i] + 1.4, f"{med[i]:.0f}%", ha="center", va="bottom",
                color=INK, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([r[0] for r in rows])
    ax.set_ylim(0, 42)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.set_title(title, loc="center", pad=8)

fig.text(0.5, 0.965, "Which input recovers a measured trade-off?",
         ha="center", va="top", fontsize=11, fontweight="bold", color=INK)
fig.text(0.5, 0.090,
         "On a new protein the score itself carries the trade-off;\n"
         "pooling it into a panel loses most of it.",
         ha="center", va="bottom", fontsize=9.5, color=INK)
fig.text(0.5, 0.025, "Median of eight scores; the white rule spans them.",
         ha="center", va="bottom", color=GREY)

save(fig, "fig_graphical_abstract")
