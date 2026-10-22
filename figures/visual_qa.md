# Final figure verification and visual review

All six main figures, both supplementary figures and the final graphical abstract are complete. This review concerns source fidelity, rendering and format; it does not assert author approval or add a scientific success threshold.

## Sources and scientific content

- Figures 3/4: all 76 cached source tables/cohort descriptions exactly match the complete refreshed postprocessing run c133ae0c-9f6a-4da6-a321-7ed9f79e585a. Every unrounded estimate, baseline, interval endpoint and support value in all 24 broad-source and 96 SPD plotted rows was independently matched to the final selectors. Original cache provenance and its reconciliation to final run02 remain explicit.
- Figure 6: the entire decoded candidate censuses for both broad sources and Nesso are byte-identical between the original and refreshed analyses. Selection, all eight actual IC50 source cells, strict relations/units/assay IDs, original predictions, source representations and variant weights are preserved. Both chosen cases have one source-variant pair; none was replaced by a mean over changing categories.
- The Nesso-only native differences share the case-B pIC50 axes because the documented oriented scalar equals nominal predicted pIC50 minus 6, and the additive constant cancels exactly. This is not empirical recalibration or a unit assignment to other readouts. The earlier separate-axis working figure and its source data are archived with hashes.
- Figure 5 explicitly labels all six full-profile curves. The supplementary level-plus-chemistry figure retains all six corresponding curves from the same saved analysis. All integer panel sizes, score-free anchors, separate query comparisons and acquisition costs remain visible; no best k or budget was chosen. Composition variation is not drawn as a sampling confidence interval. Paired uncertainty is conditional chemical-cluster resampling from the registered completed analysis.
- Score geometry, experimentally eligible identities and evaluated decision support remain distinct. The 909-identity common-channel geometry cohort is not conflated with the 930 score inputs, 928 experimentally eligible identities or readout-specific native-matched decision cohorts.

## Format and actual visual inspection

The recorded Journal of Cheminformatics requirements specify 85/170 mm PDF widths and a maximum 225 mm height. The shared initializer now uses exactly 170 mm for all eight numerical figures, replacing the generic 171.45 mm preset. Fonts, data, numerical values and vector plot primitives were preserved. Every numerical PDF is below 10 MB and within the height limit.

All eight rerendered numerical PDFs were rendered by Poppler at their final 170 mm width and inspected directly. Gene/readout labels, panel letters, legends, zero references, interval endpoints, matrix source/class strips, molecule drawings and actual IC50 tables are readable; no text collision or canvas clipping was observed. All eight Matplotlib style audits are clean. Numerical SVGs have no raster plot layers; Figure 6's four molecular picture layers are intentional, with SVG originals retained. PDF text embeds TrueType fonts without Type 3 fonts.

The optional graphical abstract retains its separate 920×300-pixel policy: white background, 34,419 bytes. Its question, comparisons, measured preferences, both generalization regimes, real molecular/protein pictures and concluding sentences are visible. Protein images are coordinate-derived illustrations, not attestations of experimental assay constructs. Draft GA files remain unchanged.

## Environment and relocated rendering

The figure environment is separately validated with Python 3.12.11, Matplotlib 3.10.6, NumPy 1.26.4, pandas 2.3.2 and Pillow 12.3.0. Core fitting/replay remains in its separate Python 3.9 environment. Default rendering uses cached molecular/protein pictures and does not import RDKit or SciPy; optional source-picture regeneration has separate pinned requirements and exact reference SMILES/PDB inputs. Source IDs, structure-publication attribution and license notes are in ASSET_PROVENANCE.md.

The compact authoring bundle passed 40 format/relocation checks. All seven normal commands regenerated all nine final artworks in an unrelated temporary location with RDKit/SciPy imports deliberately blocked; PNGs and generated captions were exact in the validated environment. This is not a cross-Matplotlib, cross-image-library or cross-platform identity claim.

The sanitized public release bundle was independently read and copied only to an owned temporary location. All seven normal commands passed; all nine PNGs matched both final authoring and public artwork, and public captions matched exactly. All 75 declared dependency hashes were valid before/after; scientific CSVs and source pictures remained exact. All 227 public-bundle predicates passed. No release file was written by this check.

All further source/export suites passed: 17 geometry/draft checks, 16 Figure 3 checks, 68 Figure 4 checks and 117 growth/example/GA checks. These are technical predicates, not independent scientific replications.

## Durable records

- `refreshed_final_figure_sources_verification.json`: complete final-source table and plotted-row equality.
- `data/figure6_choices/refreshed_source_crosscheck.json`: equality of all three complete candidate censuses.
- `final_figure_verification.json`: exact PDF dimensions, file hashes and compact relocation proof.
- `public_bundle_relocation_verification.json`: read-only sanitized public-bundle command, data and image proof.
- `figure_bundle_manifest.json`: 75 minimal relative rendering dependencies; public integration maintains derivative hashes for sanitized metadata.
- `figure_environment.json`, `FIGURE_PORTABILITY.md`, and the two figure requirements files: exact environment, relative commands, dependencies and optional asset regeneration.
- Per-figure render/verification receipts and `ASSET_PROVENANCE.md`: input/output hashes and actual data/image provenance.

## Final PDF hashes

| File | SHA-256 |
| --- | --- |
| `figure_geometry_readouts_support.pdf` | `f5e5d6220e7b078bd4da590e37a7acf6e534bacd029646559ba9e1aed9bfba04` |
| `figure_geometry_shared_shape.pdf` | `c5692cd16116a811241e24216d6779b412becc87c73bbc6f053e75fa19bcdca2` |
| `figure3_broad_paired_value.pdf` | `26407de21d912fcad4d49e01cc13506deb4f9c40c11b7011a995486ec42cb742` |
| `figure4_spd_paired_value.pdf` | `e4346126127083ffbefb1735ccc02283f523166bdfe84e627eb25c5167760710` |
| `figure5_panel_acquisition.pdf` | `71c7a40a807beb5a3a3a142cd4886d4df3d84c5272503f8610ae94a25a591d9f` |
| `figure6_measured_choice_examples.pdf` | `d63f3e3fac80d4398153d3222a9766c274556ef1e044fc4a03bc4501a891ef8a` |
| `figure_supplement_residual_correlation.pdf` | `93a8cea72e605506d1d2e2513562286ebdbd89c746e87e4bdad4b8234ae1f1c3` |
| `figure_supplement_growth_level_chemistry.pdf` | `bcdcc55a2b6f8388055ef268ad248f662a0a1bed019f99bf124dddf5a8105138` |
| `graphical_abstract_final.pdf` | `87e39e3b01557ff85c7168d59d70d4f1a81d1426eba6d366a5e6de9637838464` |

Archived draft references concern authoring evidence outside the normal public render bundle.
