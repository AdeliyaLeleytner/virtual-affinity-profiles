"""Figure 3. What buying score columns costs.

Curves are read from data/acquisition_curves.csv, for new molecules at
calibrated proteins, on the cohort that matches the support of the main
decision analysis. The y axis is the change from each panel's own k = 0 value,
so both panels share one scale and the claim is read as the shape of the curve.
Each line is one acquisition order; the lines describe composition variation,
not uncertainty.
"""
import matplotlib.pyplot as plt

from decisions import growth
from figstyle import TEXTWIDTH_IN, PALETTE, use_paper_style, figsize, share_limits, save

PANELS = [("dockstring", "a  DOCKSTRING, 58 cols"),
          ("own44", "b  Docking44, 44 cols")]

use_paper_style()
fig = plt.figure(figsize=figsize(TEXTWIDTH_IN, ratio=0.44))
fig.set_layout_engine("none")
axes = [fig.add_axes([0.115, 0.255, 0.395, 0.56]),
        fig.add_axes([0.590, 0.255, 0.395, 0.56])]
GREY, HOT, INK = PALETTE["grey"], PALETTE["orange"], PALETTE["ink"]

for ax, (ds, title) in zip(axes, PANELS):
    cur = growth(ds, "full_profile")
    base = growth(ds, "global_chemistry")
    anchor = list(base.values())[0][0]          # k = 0: chemistry alone
    for order, series in sorted(cur.items()):
        ks = sorted(series)
        fam = order.startswith("family")
        ax.plot(ks, [series[k] - anchor for k in ks],
                color=HOT if fam else GREY,
                ls=(0, (3.5, 2.0)) if fam else "-",
                lw=0.9, alpha=0.9, zorder=2)
    ax.axhline(0, color=INK, lw=0.8, zorder=3)
    ax.set_xlim(0, max(max(v) for v in cur.values()))
    ax.set_title(title, loc="left", pad=14)
    ax.set_xlabel("score columns purchased")
    ax.text(0.0, 1.015, f"chemistry alone: {anchor:.1f}%",
            transform=ax.transAxes, ha="left", va="bottom", color=GREY)

share_limits(axes)
axes[0].set_ylabel("change from chemistry alone\n(percentage points)")
axes[1].set_yticklabels([])
axes[1].text(0.97, 0.06, "solid grey: random order\ndashed orange: stratified",
             transform=axes[1].transAxes, ha="right", va="bottom", color=INK)

save(fig, "fig_acquisition")
