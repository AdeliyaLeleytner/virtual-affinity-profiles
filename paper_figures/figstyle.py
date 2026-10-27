"""Minimal, self-contained figure style for the manuscript figures.

Everything the figure scripts need, with no dependency beyond matplotlib and
numpy, so the figures can be regenerated from a clean checkout.

Figures are built at their final printed width and placed with
\\includegraphics{...} without a width= or scale= option. Scaling a figure in
LaTeX also scales its text, which is how labels end up below the legible floor.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

# Journal of Cheminformatics (Springer Nature sn-jnl class): \textwidth = 372 pt
TEXTWIDTH_IN = 372.0 / 72.27

PALETTE = {
    "ink": "#1A1A1A",      # axes, rules, body labels
    "blue": "#2C6FAD",     # series 1: the score as computed
    "orange": "#E07A3E",   # series 2: a fitted model over the profile
    "grey": "#9AA0A6",     # scaffolding and secondary curves
    "pale": "#EEF2F6",
    "white": "#FFFFFF",
}

_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "mathtext.fontset": "stixsans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "axes.edgecolor": PALETTE["ink"],
    "axes.labelcolor": PALETTE["ink"],
    "axes.titlelocation": "center",
    "axes.titlepad": 4.0,
    "axes.labelpad": 2.5,
    "axes.axisbelow": True,
    "axes.facecolor": "white",
    "axes.xmargin": 0.02,
    "axes.ymargin": 0.05,
    "axes.unicode_minus": True,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.color": PALETTE["ink"],
    "ytick.color": PALETTE["ink"],
    "text.color": PALETTE["ink"],
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.bbox": None,          # never crop: the width must stay exact
    "savefig.pad_inches": 0.0,
    "figure.constrained_layout.use": True,   # keep labels inside the fixed width
    "figure.constrained_layout.h_pad": 0.02,
    "figure.constrained_layout.w_pad": 0.02,
    "pdf.fonttype": 42,            # embed real fonts, not outlines
    "ps.fonttype": 42,
}


def use_paper_style(base_pt=8.0):
    """Apply the style. base_pt is the in-figure body size in points."""
    mpl.rcParams.update(_RC)
    mpl.rcParams.update({
        "font.size": base_pt,
        "axes.titlesize": base_pt,
        "axes.labelsize": base_pt,
        "xtick.labelsize": base_pt - 1.0,
        "ytick.labelsize": base_pt - 1.0,
        "legend.fontsize": base_pt - 1.0,
        "figure.titlesize": base_pt + 1.0,
        "figure.labelsize": base_pt,
    })
    mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
        color=[PALETTE["blue"], PALETTE["orange"], PALETTE["grey"]])


def figsize(width_in=TEXTWIDTH_IN, ratio=0.6):
    """(width, height) in inches. Default width is the journal text width."""
    return (width_in, width_in * ratio)


def share_limits(axes):
    """Force a common x and y range over axes that are meant to be compared."""
    xs = [a.get_xlim() for a in axes]
    ys = [a.get_ylim() for a in axes]
    xlim = (min(x[0] for x in xs), max(x[1] for x in xs))
    ylim = (min(y[0] for y in ys), max(y[1] for y in ys))
    for a in axes:
        a.set_xlim(xlim)
        a.set_ylim(ylim)


def save(fig, stem, outdir="."):
    """Write <stem>.pdf (vector, for the manuscript) and <stem>.png (preview)."""
    import os
    os.makedirs(outdir, exist_ok=True)
    pdf = os.path.join(outdir, stem + ".pdf")
    png = os.path.join(outdir, stem + ".png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=400)
    w, h = fig.get_size_inches()
    print(f"{pdf}  {w:.3f} x {h:.3f} in")
    plt.close(fig)
