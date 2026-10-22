# Geometry from scientific inputs

Run from any working directory:

```sh
python scripts/replay_geometry.py
```

This recomputes geometry CSVs from `data/datasets/{dockstring,own44,spd}.npz`,
then compares all 14 corresponding figure-data tables with the frozen numerical
references in `reference/geometry/`. The references are verification targets,
not the calculation inputs. Generated tables and the replay receipt are written
to `outputs/geometry/data/`. Figure scripts can consume those regenerated
CSV files. The source hashes, original OpenResearch identifiers and support-ID
hash are in `provenance/geometry_sources.json`.

## Joint DOCKSTRING/Docking44 geometry

The prepared arrays retain all 890 DOCKSTRING source rows (885 identities),
668 Docking44 identities and every score column. The supervised data loader is
not used: biological-anchor eligibility must not filter geometric support.

Join by standard molecular identity, retaining DOCKSTRING source-row order and
every DOCKSTRING identity variant. Docking44 has one row per identity. This gives
585 joined rows, 582 identities and 102 columns (58 + 44). Both original source-row
index sequences and the complete joined score array, including NaNs, were compared
with the original joint input and matched exactly. The comparison is recorded in
provenance. The common 102-column finite slice contains 384 rows/384 identities.
No variant scores are averaged and no missing score is imputed. All joined rows
remain in the support table; the common slice is explicitly labelled.

For a supported row with identity i, assign mass `1 / count(i)` and normalize masses
to sum to one. Thus each identity has equal total weight. This applies separately
on complete or pair-available support. The common slice happens to have no repeated
identities; the weighting rule still remains explicit.

The original four views are retained:

- `raw`: higher-is-stronger native scores; no clipping.
- `row_centered`: subtract each molecule's arithmetic score mean **within each
  source panel**, before column standardization.
- `rank_columns`: within-column average ranks on common support, then standardize.
- `rank_centered`: rank within columns, subtract within-source row means, then
  standardize.

Each transformed column is centered and divided by its identity-weighted population
standard deviation. SVD is applied to `sqrt(weight) * standardized_matrix`.
Eigenvalue fractions are squared singular values divided by their sum;
participation ratio is `1 / sum(fraction**2)` and PC1 is the leading fraction.
Numerical rank uses the original machine-precision SVD tolerance, not a scientific
choice of dimensionality. The primary joint view gives every column equal weight.
The saved equal-block sensitivity multiplies source columns by
`sqrt(102 / (2 * source_column_count))`.

Same-gene Pearson correlations use the weighted standardized columns on the same
384 identities; Spearman correlations use average ranks and the same identity
weights. Nine genes are shared, with raw and row-centered results retained separately.
Pair-available raw sensitivities use each gene pair's own finite support and are
not substituted for common support. The cross-source mean correlation compares
the two **raw source-specific row means**, with the same common-slice weights.
Two-way variance components retain the original source-specific weighted
molecule/target/interaction decomposition.

The reproduced joint PR values are 1.6256759428876772 (raw) and 20.351886033012843
(row-centered); raw source means correlate at 0.9728854229131486. These are
software replay targets, not performance gates or biological mechanism claims.

## SPD geometry

The prepared SPD NPZ includes all 12 archived heads, eight explicitly designated
primary affinity/energy heads, 930 identities and 10 targets. The original
`common_complete_rows_all12` mask is checked against finiteness across **all 12
heads and every target**. The ordered 909 included IDs are retained separately in
`geometry_inputs/spd_common_all12_ids.csv` and hash-bound to that NPZ. The script
also checks exact ordered identity equality. The fact that the saved eight-head
mask has the same count here does not redefine how this support was constructed.

Report the eight primary heads on this frozen common-all-12 support. For raw and
molecule-centered views, column-center the scores; compute sample covariance
(`n - 1` denominator), convert it to correlation, and diagonalize. Only negative
floating-point eigenvalue roundoff is truncated at zero under the original
numerical PSD check. The original PC1-sign convention, PR, PC1 variance fraction,
entropy rank and row-mean association definitions are preserved. Source score
orientations are those in `provenance/METHODS_TABLE.md`; no outcomes or recalibration
are involved. Inherited energy-clipping sensitivities are not substituted for
unclipped primary values.

## Verification scope

Replay checks schemas, exact categorical identities/order and every numeric entry
of the 14 retained reference tables, with absolute/relative numerical tolerance
`2e-10`. This is an arithmetic correctness tolerance, not an effect-size or
scientific acceptance threshold. No model is fit and no experimental outcome value is read. Binding-record
identity membership is used only to regenerate the support-count panel; no
concentration, bound or label value is loaded. No docking/inference is launched. Family-label permutations and earlier
experimental associations are outside this figure-data replay. The original joint
geometry is recorded as OpenResearch project `42ba66d2-a533-4aa8-9864-ecb2b5986391`,
experiment `49f96c0a-7229-4474-990d-d28eeaa51654`, completed run
`8d42f171-3639-4139-850a-67309cec8dcd`, commit `667e161`.

## Render directly from fresh geometry tables

Run the numerical replay in the pinned analysis environment, then run the wrapper
in the separate pinned figure environment:

```sh
python scripts/replay_geometry.py
# Activate the Python 3.12 figure environment.
python scripts/render_geometry_from_replay.py
```

The wrapper copies only the10 required fresh CSVs and the two unchanged plotting
/style sources into `outputs/geometry_figures/`, then renders there. It never
falls back to cached CSVs. All three geometry PNG pixel arrays were checked
against the final170mm cache and are identical in the validated environment; see
`provenance/fresh_geometry_figure_check.json`. Normal cached rendering remains
available under `figures/RENDERING.md`.
