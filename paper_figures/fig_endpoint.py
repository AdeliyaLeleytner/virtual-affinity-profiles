"""Figure: one number is not enough.

a  The anatomy of a single measured decision. Two near-identical molecules,
   the four measured cells behind the contrast, and where each input places
   the pair in the plane of the two changes it has to get right.
b  The same two criteria across all eight scores on proteins withheld from
   calibration. Agreeing with the sign of D separates the inputs barely;
   recovering both directions separates them more than tenfold.

Panel a reads the released worked case in figures/data/figure6_choices/.
Panel b is recomputed from data/decisions_paired_uncertainty.csv (Additional
file 5) under the exact four-cell subset with actual-assay-pair weighting.
Nothing here is transcribed.
"""
import csv
import json
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyArrowPatch

from figstyle import TEXTWIDTH_IN, PALETTE, use_paper_style, figsize, save

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)
CASES = os.path.join(ROOT, "figures", "data", "figure6_choices")
ASSETS = os.path.join(ROOT, "figures", "assets", "figure6_choices")
DECISIONS = os.path.join(HERE, "data", "decisions_paired_uncertainty.csv")

CASE = "B"            # proteins withheld from calibration, Nesso-1
INK, BLUE, ORANGE, GREY, PALE = (PALETTE["ink"], PALETTE["blue"],
                                 PALETTE["orange"], PALETTE["grey"],
                                 PALETTE["pale"])

# the three inputs the paper names, in the colours Figures 1 and 2 give them
INPUTS = [("global_chemistry", "chemistry alone", GREY, "D"),
          ("full_profile", "the whole panel", ORANGE, "s"),
          ("native_score", "that protein's own column", BLUE, "o")]
SHORT = {"global_chemistry": "chemistry", "full_profile": "whole panel",
         "native_score": "own column"}

NAME = {"raw_ad4": "AutoDock4", "raw_gnina_own": "GNINA (own search)",
        "vina": "Vina", "gnina": "GNINA (Vina pose)",
        "boltz_affinity": "Boltz-2", "boltzina_vina_far_affinity": "Boltzina (far)",
        "boltzina_vina_zero_affinity": "Boltzina (zero)",
        "nesso_affinity": "Nesso-1"}
ORDER = ["raw_ad4", "raw_gnina_own", "vina", "gnina", "boltz_affinity",
         "boltzina_vina_zero_affinity", "boltzina_vina_far_affinity",
         "nesso_affinity"]


def rows(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh))


# ── panel a: the worked case ────────────────────────────────────────────
case = [r for r in rows(os.path.join(CASES, "selected_cases.csv"))
        if r["case"] == CASE][0]
cells = [r for r in rows(os.path.join(CASES, "actual_four_cells.csv"))
         if r["case"] == CASE]
bound = [r for r in rows(os.path.join(CASES, "observed_delta_bounds.csv"))
         if r["case"] == CASE][0]
preds = {r["model"]: (float(r["delta_t"]), float(r["delta_u"]))
         for r in rows(os.path.join(CASES, "all_source_variant_decisions.csv"))
         if r["case"] == CASE}

GENE_T = case["gene_a"]
GENE_U = case["gene_b"]
ic50 = {(c["molecule_end"], c["target"]): (c["ic50_uM"], c["source_relation"])
        for c in cells}

dt_lo, dt_hi = float(bound["delta_t_lower"]), float(bound["delta_t_upper"])
du_hi = float(bound["delta_u_upper"])
du_open = bound["delta_u_upper_open"] == "True"


# ── panel b: the same two criteria over the panel ───────────────────────
DEC = rows(DECISIONS)
BASE = dict(measurement_stratum="exact_four_cell",
            weighting="actual_assay_pair_equal",
            domain="all_source_operational",
            resampling="chemical_group_fixed_panel",
            dataset="spd", regime="new_proteins_known_molecules",
            comparison_support="native_matched_common")


def value(channel, model, metric, field="estimate"):
    for r in DEC:
        if (all(r[k] == v for k, v in BASE.items())
                and r["channel"] == channel and r["model"] == model
                and r["metric"] == metric):
            return 100 * float(r[field])
    raise KeyError((channel, model, metric))


def series(metric):
    """chemistry / whole panel / own column, one triple per score."""
    out = {}
    for ch in ORDER:
        out[ch] = (value(ch, "native_score", metric, "baseline_estimate"),
                   value(ch, "full_profile", metric),
                   value(ch, "native_score", metric))
    return out


SIGN = series("D_direction_concordance")
BOTH = series("pareto_correct_fraction")
span = lambda d: (min(min(v) for v in d.values()), max(max(v) for v in d.values()))

use_paper_style()
fig = plt.figure(figsize=figsize(TEXTWIDTH_IN, ratio=1.16))
fig.set_layout_engine("none")

# ── a, left: the molecules and what was measured ────────────────────────
for j, end in enumerate(("a", "b")):
    ax = fig.add_axes([0.035 + j * 0.215, 0.790, 0.185, 0.150])
    ax.imshow(mpimg.imread(os.path.join(ASSETS, f"{CASE}_{end}.png")))
    ax.axis("off")
    ax.set_title(case[f"name_{end}"], fontsize=7.5, color=INK, pad=1)

arrow = fig.add_axes([0.221, 0.790, 0.030, 0.150])
arrow.axis("off")
arrow.annotate("", xy=(0.95, 0.5), xytext=(0.05, 0.5),
               arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.9,
                               mutation_scale=7))

tab = fig.add_axes([0.035, 0.575, 0.430, 0.200])
tab.set_xlim(0, 1)
tab.set_ylim(0, 1)
tab.axis("off")
tab.text(0.0, 0.93, f"Morgan Tanimoto {float(case['chemical_tanimoto']):.2f}, "
                    "one chemical group", fontsize=7.5, color=GREY, va="center")
tab.text(0.0, 0.70, "measured IC$_{50}$ (µM)", fontsize=7.5, color=INK,
         va="center")
COLX = (0.44, 0.76)
for x, g in zip(COLX, (GENE_T, GENE_U)):
    tab.text(x, 0.48, g, fontsize=7.5, color=INK, ha="center", va="center")
tab.plot([0.0, 1.0], [0.37, 0.37], color=GREY, lw=0.6, clip_on=False)
for i, end in enumerate(("a", "b")):
    y = 0.24 - i * 0.19
    tab.text(0.0, y, case[f"name_{end}"], fontsize=7.5, color=INK, va="center")
    for x, g in zip(COLX, (GENE_T, GENE_U)):
        v, rel = ic50[(end, "t" if g == GENE_T else "u")]
        txt = ("> " if rel == ">" else "") + f"{float(v):g}"
        tab.text(x, y, txt, fontsize=7.5, ha="center", va="center",
                 color=ORANGE if rel == ">" else INK,
                 fontweight="bold" if rel == ">" else "normal")
tab.text(0.0, -0.20, "bold: censored, the exact value is unknown",
         fontsize=6.8, color=ORANGE, ha="left", va="center", clip_on=False)

# ── a, right: the plane of the two changes ──────────────────────────────
ax = fig.add_axes([0.580, 0.610, 0.385, 0.330])
XL, XR, YB, YT = -0.28, 0.38, -0.33, 0.11
ax.add_patch(plt.Rectangle((0, YB), XR, -YB, facecolor=PALE, edgecolor="none",
                           zorder=0))
ax.text(0.014, YB + 0.013, "correct", fontsize=6.8,
        color=GREY, ha="left", va="bottom", zorder=1)
ax.plot([XL, XR], [XL, XR], ls=(0, (4, 3)), color=GREY, lw=0.8, zorder=1)
ax.text(-0.252, -0.238, "D = 0", fontsize=6.8, color=GREY, rotation=33)
ax.axhline(0, color=INK, lw=0.7, zorder=2)
ax.axvline(0, color=INK, lw=0.7, zorder=2)

# what was observed: one exact change, one open bound
ax.add_patch(FancyArrowPatch((dt_lo, du_hi), (dt_lo, YB + 0.012),
                             arrowstyle="-|>", mutation_scale=7,
                             color=INK, lw=1.3, zorder=5))
ax.plot([dt_lo], [du_hi], "o", mfc="white", mec=INK, mew=1.0, ms=4.4, zorder=6)
ax.text(dt_lo - 0.014, du_hi + 0.006, "measured", fontsize=7, color=INK,
        ha="right", va="bottom")
ax.text(dt_lo - 0.014, du_hi - 0.028, "anywhere\nbelow", fontsize=6.8,
        color=GREY, ha="right", va="top", linespacing=1.4)

for key, label, col, mk in INPUTS:
    x, y = preds[key]
    ax.plot([x], [y], mk, color=col, ms=4.2, mec=col,
            mfc=col if mk != "s" else "none", mew=1.0, zorder=6)
    off = {"global_chemistry": (7, 3), "full_profile": (7, -3),
           "native_score": (7, 3)}[key]
    ax.annotate(SHORT[key], (x, y), textcoords="offset points", xytext=off,
                fontsize=7, color=col, ha="left", va="center")

ax.set_xlim(XL, XR)
ax.set_ylim(YB, YT)
ax.set_xlabel(f"change at {GENE_T}  (Δ pIC$_{{50}}$)")
ax.set_ylabel(f"change at {GENE_U}  (Δ pIC$_{{50}}$)")
ax.set_xticks([-0.2, 0.0, 0.2])
ax.set_yticks([-0.3, -0.2, -0.1, 0.0, 0.1])

fig.text(0.035, 0.968, "a  One decision, " + GENE_T + " withheld from "
         "calibration", fontsize=9, fontweight="bold", color=INK, va="center")
fig.text(0.035, 0.520,
         "All three inputs put the pair on the correct side of D = 0. "
         "Only the queried column puts\nboth changes on the correct side of "
         "zero, which is what the endpoint asks.",
         fontsize=7.6, color=INK, va="top", linespacing=1.5)

# ── b: the same two criteria across the panel ───────────────────────────
fig.text(0.035, 0.432, "b  The two criteria across all eight scores, "
         "proteins withheld", fontsize=9, fontweight="bold", color=INK,
         va="center")

y = np.arange(len(ORDER))
axes = []
for k, (data, title) in enumerate(((SIGN, "sign of D recovered"),
                                   (BOTH, "both directions recovered"))):
    ax = fig.add_axes([0.300 + k * 0.345, 0.115, 0.300, 0.250])
    axes.append(ax)
    lo, hi = span(data)
    for i, ch in enumerate(ORDER):
        trip = data[ch]
        ax.plot([min(trip), max(trip)], [i, i], color=PALE, lw=3.4, zorder=1,
                solid_capstyle="round")
        for (key, _, col, mk), v in zip(INPUTS, trip):
            ax.plot([v], [i], mk, color=col, ms=3.6, mec=col,
                    mfc=col if mk != "s" else "white", mew=0.9, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([NAME[c] for c in ORDER] if k == 0 else [])
    ax.set_ylim(len(ORDER) - 0.4, -0.9)
    ax.set_title(f"{title}\nspan {hi - lo:.0f} points, {hi / lo:.1f}×",
                 loc="left", pad=6, fontsize=8)
    ax.set_xlabel("% of identified trade-offs" if k else "% of pairs")
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)

axes[0].set_xlim(48, 68)
axes[0].set_xticks([50, 55, 60, 65])
axes[1].set_xlim(0, 38)
axes[1].set_xticks([0, 10, 20, 30])

for (key, label, col, mk) in INPUTS:
    pass
handles = [plt.Line2D([], [], color=c, marker=m, ls="none", ms=3.6,
                      mfc=c if m != "s" else "white", mec=c, label=lab)
           for _, lab, c, m in INPUTS]
fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.135, 0.005),
           ncol=3, frameon=False, fontsize=7.2, handletextpad=0.4,
           columnspacing=2.0)

save(fig, "fig_endpoint")
print(f"sign of D   span {span(SIGN)[0]:.1f}-{span(SIGN)[1]:.1f}")
print(f"both dirs   span {span(BOTH)[0]:.1f}-{span(BOTH)[1]:.1f}")
