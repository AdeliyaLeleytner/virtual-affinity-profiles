# Final figure rendering handoff

The eight numerical figure PDFs use exactly 170 mm width through the single `orx_figstyle.py` initializer. Their text, axes, intervals and matrix cells remain vector; Figure 6 intentionally embeds four cached molecular pictures. The optional graphical abstract follows its separate 920×300-pixel, white-background, ≤150 KB specification. Archived drafts are not normal render targets.

Use the separate validated **Python 3.12** figure environment in `requirements-figures-validated.txt`. It is independent of the **Python 3.9** fitting/replay environment. The normal commands below require Matplotlib, NumPy, pandas and Pillow, and use cached molecular/protein pictures. RDKit and SciPy are optional source-picture regeneration dependencies only. Recorded pixel comparisons establish relocation identity within the validated environment, not across Matplotlib/library/OS versions.

## Commands

Run from the directory containing the figure scripts (`figures/` in a public package). Outputs are always written beside the scripts. Paths below are relative; no upstream analysis workspace is accessed by these rendering modes.

```sh
python build_geometry_figures.py
python build_primary_result_figures.py --figure broad --from-cache
python build_primary_result_figures.py --figure spd --from-cache
python build_growth_figure.py --analysis-dir data/figure5_sources --orx-receipt data/figure5_sources/23e36b7b-fc22-450a-8f6e-01b20fa8e64e.log
python build_growth_figure.py --analysis-dir data/figure5_sources --orx-receipt data/figure5_sources/23e36b7b-fc22-450a-8f6e-01b20fa8e64e.log --model global_chemistry
python build_choice_examples.py --mode render
python build_graphical_abstract_final.py
```

`write_captions.py` optionally regenerates the historical geometry-caption draft. Every main result/choice/growth/GA renderer writes its own caption and source receipt. The growth commands deliberately reuse the same cached source tree for both arms; `growth_level_sources` is an equivalent authoring cache, not a second data requirement.

## Required files beside scripts

- `orx_figstyle.py` and the five `build_*.py` scripts used above. The project width/font modifications are in this one style initializer.
- Geometry inputs: `data/spd_readouts.csv`, `panel_composition.csv`, `cohort_support.csv`, `shared_means.csv`, `shared_mean_correlations.csv`, `full_spectra.csv`, `spectral_summary.csv`, `shared_gene_main_figure.csv`, `column_metadata.csv`, and `correlation_row_centered_complete.csv`. Optional geometry captions also use `geometry_support.csv`.
- Figure 3/4: their existing `*_render_receipt.json` plus `data/<figure-name>_values.csv`, `_within_assay_absolute.csv`, `_native_decision_absolute.csv`, and `_all_weightings_domains.csv`. Figure 4's caption also reads `data/spd_readouts.csv`. The `data/figure3_sources` and `data/figure4_sources` snapshots preserve full source-table provenance; the normal cached renderer does not read an absolute upstream path.
- Growth: `data/figure5_sources/growth_analysis_receipt.json`, `all_growth_preferences.csv`, the retained registered run log, and for each `dockstring`/`own44` subdirectory: `all_measurements_paired_uncertainty.csv`, `full_budget_invariance_checks.csv`, `query_assay_pair_costs.csv`, `query_profile_union_costs.csv`.
- Figure 6: `data/figure6_choices/selected_cases.csv`, `actual_four_cells.csv`, `all_source_variant_decisions.csv`, `selection_receipt.json`, `nesso_unit_identity.json`; all files in `assets/figure6_choices`. Reference SMILES, exact experimental records, native algebra and every selected source-variant weight remain explicit. `--mode prepare` is the historical source-census selection stage, not required for public cached rendering.
- Final GA: `data/molecule_illustrations.csv`, `illustration_profiles.csv`, `protein_illustrations.csv`; cached butenafine/indomethacin molecule PNG/SVG files, 6CM4/6PDJ/3B66 backbone PNG/SVG files and retained PDB files; `assets/structure_retrieval.json`. Retain coordinate entry JSON files and the original `graphical_abstract_render_receipt.json` for attribution/chain provenance.
- `ASSET_PROVENANCE.md`, environment requirements, captions and verification JSONs for the complete provenance trail. No third-party paper figure is needed.

## Fresh saved-prediction replay → plotted inputs

These commands run from the **public repository root**, after restoring the frozen prediction assets and rebuilding replay manifests as described in `REPLAY.md`. The public replay configs read restored predictions under `artifacts/frozen_predictions/` and write new analysis tables under `outputs/`. First use the separate core/replay environment; capture the completed growth-analysis stdout log:

```sh
python scripts/replay_results.py postprocessing
python scripts/replay_results.py growth > growth-analysis.log
```

Then use the validated Python 3.12 figure environment:

```sh
python figures/build_primary_result_figures.py --figure broad --results-root outputs/profile_postprocessing_run02 --prediction-root artifacts/frozen_predictions/core_run02
python figures/build_primary_result_figures.py --figure spd --results-root outputs/profile_postprocessing_run02 --prediction-root artifacts/frozen_predictions/core_run02
python figures/build_growth_figure.py --analysis-dir outputs/panel_growth_analysis --run-log growth-analysis.log
python figures/build_growth_figure.py --analysis-dir outputs/panel_growth_analysis --run-log growth-analysis.log --model global_chemistry
python figures/build_choice_examples.py --mode prepare --input-dir data --prediction-root artifacts/frozen_predictions/core_run02 --postprocessing-root outputs/profile_postprocessing_run02 --data-loader replay_src/profile_information/data.py
python figures/build_choice_examples.py --mode render
```

Figures 3/4 refresh all plotted tables directly from the newly completed replay; `--prediction-root` supplies each source input description without trusting an old `dataset_state.source` path. Figure 5 reads the newly aggregated full growth tables; `--run-log` accepts the completed analysis log (`--orx-receipt` remains an archival alias). Figure 6 regenerates deterministic selected cases, actual source cells and every variant-specific decision from the complete replay census, prepared inputs and saved predictions; all four source paths are mandatory for prepare/both mode. It has no private-parent fallback. If new predictions select different molecules, explicitly regenerate their optional RDKit pictures before rendering; fixed archived predictions reproduce the supplied selected numeric records exactly.

For newly fitted prediction directories rather than restored archives, pass the actual core dataset root in `--prediction-root` and update the replay config/manifests accordingly. These plotting/selection commands themselves perform no predictive fits, resampling or new performance estimation. Cached rendering remains the default short path. The fresh-path check uses read-only aliases of the completed outputs at these public relative locations and writes only an owned temporary directory.

## Optional source-picture regeneration

```sh
python build_choice_examples.py --mode render --regenerate-molecules
python regenerate_illustration_assets.py
```

These require `requirements-figure-assets-optional.txt` and retained exact SMILES/PDB inputs. They regenerate source pictures without predictive fitting or structure prediction. Normal rendering uses the already verified cached picture files. Source software/data attribution and PDB publication citations are in `ASSET_PROVENANCE.md`.

The manuscript-local `verify_*` scripts also check original absolute provenance sources where available. Those source-audit paths are metadata, not rendering dependencies; a public integration should use the relative cached-input checks when the original workspace is absent.
