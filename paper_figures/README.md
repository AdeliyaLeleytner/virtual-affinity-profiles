# Figures for *Affinity panels dilute the target-specific signal that transfers to a held-out protein*

Regenerates the three manuscript figures and the graphical abstract from the
tables published with the article. Nothing here contains a transcribed number:
every value is read from a file in `data/` or recomputed from the frozen score
arrays in `../data/datasets/`.

## Run

```
cd paper_figures
python3 geometry.py            # prints the geometry table and checks it against the paper
python3 fig_reversal.py        # Figure 1
python3 fig_geometry.py        # Figure 2
python3 fig_acquisition.py     # Figure 3
python3 fig_graphical_abstract.py
```

Requires only `numpy` and `matplotlib`. Each script writes `<name>.pdf`
(vector, the version placed in the manuscript) and `<name>.png` (preview).

`geometry.py` run on its own prints its results beside the published values, so
a reader can see immediately whether the checkout reproduces them:

```
  recomputed: 1.43-2.93, PC1 55.1-83.3%
  manuscript: 1.43-2.93, PC1 55.1-83.3%
```

## Files

| File | What it is |
|---|---|
| `figstyle.py` | the figure style: palette, sizing, saving. No external style dependency |
| `geometry.py` | effective dimensionality from `../data/datasets/*.npz` |
| `decisions.py` | loaders for the measured decision outcomes |
| `fig_reversal.py`, `fig_geometry.py`, `fig_acquisition.py` | Figures 1-3 of the article |
| `fig_graphical_abstract.py` | a 920 x 300 px graphical abstract to the journal's specification |
| `fig_chemistry_knows.py` | how much of a panel's leading direction six ordinary descriptors carry; prints the two values quoted in Section 2.4 |
| `fig_endpoint.py`, `fig_shared_shape.py` | two further figures built from the same released tables, not placed in the article |
| `data/decisions_paired_uncertainty.csv` | Additional file 5 of the article |
| `data/acquisition_curves.csv` | Additional file 4 of the article |

`fig_chemistry_knows.py` prints its two headline numbers when run:

```
shared R2 0.769 | residual median 0.136 range -0.057..0.514 | above 0.3: 15/102
```

The first is the out-of-fold coefficient of determination for a molecule's mean
score across the joint 102-column panel, regressed on the six source
descriptors over five folds that hold out whole chemical groups. The second is
the median of the same quantity fitted once per column to what is left after
that molecule's own panel mean is removed. `fig_endpoint.py` additionally reads
the released worked case in `../figures/data/figure6_choices/` and the molecular
depictions in `../figures/assets/figure6_choices/`.

## Two things worth knowing before you change anything

**Figures are built at their final printed width.** The journal's text block is
372 pt, which is `figstyle.TEXTWIDTH_IN`. They are placed with
`\includegraphics{...}` and no `width=` or `scale=`. Scaling a figure in LaTeX
scales its text too, which is how 8 pt labels end up below the legible floor.

**The acquisition tables hold two cohorts.** `decisions.growth` selects
`all_source_operational`, which matches the support of the main decision
analysis. The smaller `metadata_matched_operational` cohort is in the same file
and gives different numbers. The two are never averaged, and anything you
report should say which one it used.
