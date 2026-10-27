"""Figure 1. The reversal between the two transfer settings.

Every plotted value is read from data/decisions_paired_uncertainty.csv: the
exact four-cell measurement subset, with equal weight per actual assay pair,
and 2.5th/97.5th percentiles from 400 paired chemical-group resampling draws.
Nothing is synthesised.
"""
import numpy as np
import matplotlib.pyplot as plt

from decisions import NAME, HELD_OUT_PROTEIN, NEW_MOLECULE, trade_offs
from figstyle import TEXTWIDTH_IN, PALETTE, use_paper_style, figsize, share_limits, save

nat_h = trade_offs(HELD_OUT_PROTEIN, "native_score")
prof_h = trade_offs(HELD_OUT_PROTEIN, "full_profile")
nat_n = trade_offs(NEW_MOLECULE, "native_score")
prof_n = trade_offs(NEW_MOLECULE, "full_profile")

# order by the held-out-protein gain: the axis the argument runs along
order = sorted(NAME, key=lambda c: nat_h[c]["d"])

use_paper_style()
fig, (axA, axB) = plt.subplots(
    2, 1, figsize=figsize(TEXTWIDTH_IN, ratio=0.88), sharex=True, sharey=True)
COOL, HOT, GREY, INK = (PALETTE["blue"], PALETTE["orange"],
                        PALETTE["grey"], PALETTE["ink"])


def panel(ax, nat, prof, title):
    ax.axvline(0, color=INK, lw=0.8, zorder=1)
    for i, ch in enumerate(order):
        for d, col, off, mk, fill in ((prof[ch], HOT, +0.17, "s", "none"),
                                      (nat[ch], COOL, -0.17, "o", None)):
            ax.plot([d["lo"], d["hi"]], [i + off] * 2, color=col, lw=1.0,
                    solid_capstyle="butt", zorder=2)
            ax.plot([d["d"]], [i + off], mk, color=col, ms=3.2,
                    mfc=fill or col, mec=col, mew=0.9, zorder=3)
    base = [d["base"] for d in nat.values()]
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([NAME[c] for c in order])
    ax.set_ylim(-1.4, len(order) - 0.4)
    ax.invert_yaxis()
    ax.set_title(title, loc="left", pad=11)
    ax.text(0.0, 1.015,
            f"chemistry baseline recovers "
            f"{min(base):.1f}–{max(base):.1f}% of trade-offs",
            transform=ax.transAxes, ha="left", va="bottom", color=GREY)


panel(axA, nat_h, prof_h, "a  Protein withheld from calibration")
panel(axB, nat_n, prof_n, "b  New molecule, protein calibrated")

for ax in (axA, axB):
    ax.set_xlim(-26, 39)
axA.text(-24.0, -1.18, "the protein's own score", color=COOL,
         ha="left", va="center")
axA.text(-24.0, -0.62, "the whole profile", color=HOT,
         ha="left", va="center")

share_limits([axA, axB])
fig.supxlabel("correct trade-offs recovered, points vs the chemistry baseline")
save(fig, "fig_reversal")
