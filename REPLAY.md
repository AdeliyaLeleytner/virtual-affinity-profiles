# Frozen prediction replay and independent refitting

Two reproduction paths are provided. Prediction replay uses archived scientific
arrays and reruns descriptive calculations without fitting. Refitting uses the
same prepared inputs, exact fitting source and portable configurations to create
new fits; it requires CPU time and does not assert identical platform arithmetic.
Original GPU docking/co-folding inference is not rerun in either path. Model,
checkpoint, orientation and runtime provenance/known gaps are recorded separately.

## Exact stage packages

| Stage | Package/source | Entry point |
|---|---|---|
| Core, all-binding, both-cold fits | `src/profile_information/` | `scripts/run_fit.py` |
| Canonical growth fits | `growth_src/profile_information/` | `scripts/run_fit.py growth` |
| Secondary decoder fits | `decoder_src/profile_information/` | `scripts/run_fit.py decoder` |
| Fixed prediction replay | `replay_src/profile_information/`, `postprocessing/` | `scripts/replay_results.py postprocessing` |
| Measurement strata | same fixed replay package | `scripts/replay_results.py measurement-strata` |
| Growth aggregation/replay | `growth_src/`, `growth_analysis/`, `growth_analysis_support/` | `scripts/replay_results.py growth` |

The fixed replay helper preserves model names containing `__` when deriving
compact contexts. Its isolated package retains the authoritative parser repair;
original fit packages are byte-identical to their source worktrees. Growth uses
the exact numerical aggregation support recorded for its own analysis commit.
`all` runs stages in fresh interpreters to avoid cross-stage module caching.
Source hashes/commits are in `provenance/replay_source_packages.json`.

## Assets, manifests and metadata derivatives

`tools/stage_frozen_results.py` checks terminal run receipts and every original
completion artifact hash before copying. The complete core audit pack retains
all original fit artifacts, including coefficients and record-level predictions;
this lets the exact decoder module verify its frozen inheritance. Other packs
omit optional coefficients/record tables while retaining held-out predictions,
splits, compact metrics, diagnostics and source completion records.

Original absolute paths are not public dependencies. Portable effective configs
and historical signature references are labelled as metadata derivatives; the
historical signature ID does not claim to be a hash of a reduced public JSON
body. Original input/module byte manifests remain available for exact checks.

Decoder `baseline_inheritance.json` contains a machine-specific frozen-directory
pointer. Its public derivative changes only that pointer to a repository-relative
asset path. The byte-identical original completion is retained as
`source_complete.original.json`; public `complete.json` records both original
and derivative metadata hashes. Original predictions, splits and compact metrics
are checked against the original completion, not merely the derivative manifest.
The original path-bearing metadata file is not redistributed. This transformation
does not change any model, source row, index, prediction, tuning or metric.
`tools/verify_frozen_assets.py` checks this chain explicitly for all 60 splits.

`tools/package_assets.py create` builds deterministic archives with relative paths,
zeroed machine owner/time headers and a checksum catalog. `restore ARCHIVE`
requires the catalog hash and restores only regular files/directories under the
repository asset location; unexpected paths are rejected. Restoring all five
archives into a fresh relocated repository, verifying every artifact/module/input,
replaying geometry and running the exact decoder inheritance check succeeded;
see `provenance/relocation_restore_qa.json`.

After restoration, `tools/build_replay_manifests.py` verifies original scientific
hashes and pins the exact portable file membership/bytes. Public replay checks
those hashes before and after analysis. Each generated output receipt also
records its local sources/code. Outputs are local generated artifacts; do not
publish runtime receipts containing machine paths without an explicit metadata
sanitization record.

## Calculations and outputs

| Stage | Output directory |
|---|---|
| Fixed core/all-binding/both-cold/decoder replay | `outputs/profile_postprocessing_run02/` |
| Exact four-cell / censor-involving identified contrasts | `outputs/profile_measurement_strata_run02/` |
| Growth curves, acquisition costs, paired zero/query/full endpoints | `outputs/panel_growth_analysis/` |
| Geometry tables and figures | `outputs/geometry/data/`, `outputs/figures/` |

Within each actual assay pair, rectangle contributions are pooled before equal
weighting across defined assay pairs. Raw rectangle weighting is a separately
labelled sensitivity. Source variants are integrated after each source-row pair's
decision; squared moments remain separate from squared mean predictions. Every
arm-minus-baseline comparison uses paired finite support. Native scores retain
directional/association metrics and do not receive uncalibrated pIC50 RMSE.

Conditional descriptive uncertainties resample observed chemical-group endpoints,
observed gene endpoints, or both, separately. A dyad of distinct groups receives
the product of endpoint counts; a same-group dyad receives one count. No rectangle
is treated as an independent sampling unit. Association comparisons use the
intersection of defined assay contexts for the arm and baseline in each draw.
Gene resampling is panel-composition sensitivity; these distributions do not
include refitting uncertainty, new assay calibration or pretrained exposure.
The fixed 400 resamples/seed and numerical replay tolerances reproduce the
calculation and do not impose scientific acceptance rules.

Growth retains every integer k, all five folds and six predeclared acquisition
orders. k=0 is an exact score-free chemistry anchor. Task-specific matching-score
and training-selected-score comparators are not labelled as a common fixed k=1
panel; their training search and per-query/profile acquisition costs remain
explicit. Full-panel equal-input comparisons are numerical invariance checks,
not biological composition variability. No best k is selected by this release.
