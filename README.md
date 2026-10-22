# When do virtual affinity profiles improve biological decisions?

The studied score panels share a dominant axis associated with the overall
molecular score level. Additional residual dimensions do not consistently improve
experimentally measured preference changes beyond a strong chemistry/global
baseline. Useful signals depend on the query-target column, decoder and transfer
regime, with positive cases retained alongside the limitations.

These are the reproducibility materials prepared for a **Journal of
Cheminformatics** manuscript by **Adeliya Leleytner and Maxim Fedorov**. There is
no published article DOI. The manuscript and submission documents remain separate
from this code/data/results/figures release.

Use chemistry without any scores as the acquisition-cost anchor. Then compare
the chemistry-plus-global-level baseline, matching query-target score and full
residual profile at their respective acquisition costs. Evaluate
them on the intended measured decision and transfer regime. The results distinguish
new molecules/known proteins, new proteins/known molecules and supplementary
both-cold evaluation; they do not certify pretrained-model non-exposure or equate
operational IC50 contrasts with thermodynamic Kd selectivity.

## Inspect the results and figures

Compact completed results are included under [`results/`](results/), including
[paired primary estimates](results/profile_postprocessing_run02/),
[exact/censor-involving measurement strata](results/profile_measurement_strata_run02/),
[canonical acquisition curves](results/panel_growth_analysis/) and
[actual decision denominators](results/decision_support/).
[`figures/`](figures/) contains the six main figures, graphical abstract,
supplementary figures, renderer code and cached numerical/image inputs.

```sh
python tools/verify_compact_results.py
```

This verifies the bundled compact subset. Raw rectangle and resampling-draw
censuses are regenerated from the archived predictions, rather than duplicated
in the repository. Source/private-path metadata derivatives are explicit in
[`provenance/completed_output_assets.json`](provenance/completed_output_assets.json).
The counts in the support tables are denominators, not independent replicates.

## Install the analysis environment

The verified analysis runtime is Python 3.9.23 with exact versions in
`requirements.txt`. The declared range is Python 3.9–3.12; numerical reproduction
on untested Python/BLAS combinations is not assumed.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install --no-deps -e .
python scripts/verify_bundle.py
python tools/run_checks.py
python scripts/replay_geometry.py
```

The tests include small synthetic mathematical fits, not scientific-study
refits. Geometry is recomputed from scientific arrays and checked against 14
reference tables. No GPU, checkpoint download or original inference environment
is needed for these checks. Exact verification scope is recorded in
[`RELEASE_STATUS.md`](RELEASE_STATUS.md) and [`provenance/`](provenance/).

## Download archived predictions and replay analyses

The five source archives are versioned assets of
[v1.0.0](https://github.com/AdeliyaLeleytner/virtual-affinity-profiles/releases/tag/v1.0.0).
Their filenames, sizes, download addresses and SHA-256 checksums are pinned in
[`provenance/archive_assets.json`](provenance/archive_assets.json).

```sh
python tools/download_assets.py
for pack in release_assets/*.tar.gz; do python tools/package_assets.py restore "$pack"; done
python tools/verify_frozen_assets.py
python tools/build_replay_manifests.py
python scripts/replay_results.py all
```

The downloader verifies existing files and skips downloading them again. Restoring
archives uses relative paths and rejects unexpected paths. The complete core
archive retains the original fit audit needed for decoder inheritance; other
archives contain compact prediction replay assets. Expanded archives and full
per-rectangle outputs require substantially more storage than the compact repo.

A small representative public replay can use only the fingerprint archive:

```sh
python tools/download_assets.py --run fingerprint_run01
python tools/package_assets.py restore release_assets/fingerprint_run01.tar.gz
python tools/build_replay_manifests.py --stage postprocessing --config configs/postprocessing_smoke.json
python scripts/replay_results.py postprocessing --config configs/postprocessing_smoke.json
python scripts/replay_results.py measurement-strata --config configs/measurement_strata_smoke.json
```

This uses all 10 saved splits of the selected Docking44 all-binding dataset, both
regimes and 400 paired draws. It is a bounded usability example, not the full study.

Separate replay stages are `postprocessing`, `measurement-strata` and `growth`.
Outputs are written under `outputs/`; scientific predictive models are not fit by
these commands. [`REPLAY.md`](REPLAY.md) describes the frozen versus patched
packages, portable completion manifests, weighting and conditional resampling.
Original fit packages retain their exact hashes. The repaired descriptive parser
lives separately in `replay_src/`. Gene resampling describes panel-composition
sensitivity; no nominal population confidence coverage is asserted.

The authoritative completed analyses and the bounded public-package checks have
different scopes. A previous duplicate public full replay encountered storage
exhaustion; its exact limitation remains in provenance. Complete authoritative
run02 outputs are bundled as compact summaries. A new full public replay is not
claimed by the bounded installation/source-hash/geometry checks.

## Render figures

Figure rendering has its own pinned Python 3.12 environment; see
[`figures/RENDERING.md`](figures/RENDERING.md) and `requirements-figures.txt`.
Default rendering uses cached molecule/structure images and needs no RDKit,
structure prediction, image generation, network retrieval or model fitting.
[`GEOMETRY_REPLAY.md`](GEOMETRY_REPLAY.md) explains the geometry calculation and
its distinct observation cohorts. For the raw-score-to-fresh-table-to-figure
chain, run `scripts/replay_geometry.py` in the analysis environment, then
`scripts/render_geometry_from_replay.py` in the figure environment.

## Explicit scientific refits

These commands are optional and can require substantial CPU time:

```sh
python scripts/run_fit.py core --cost-only
python scripts/run_fit.py core
python scripts/run_fit.py fingerprint
python scripts/run_fit.py both-cold
python scripts/run_fit.py growth
python scripts/run_fit.py decoder
```

`--cost-only` checks inputs and the contract without fitting. Core/all-binding/
both-cold use frozen `src/`; canonical growth uses `growth_src/`; decoder uses
`decoder_src/` and the full core audit assets. Refits write fresh signatures and
do not resume the original private run directories. No best panel size or
scientific success threshold is selected by the release.

## Citation, attribution and licences

Repository: [AdeliyaLeleytner/virtual-affinity-profiles](https://github.com/AdeliyaLeleytner/virtual-affinity-profiles).
Version 1.0.0 citation metadata and both authors' ORCIDs are in
[`CITATION.cff`](CITATION.cff). No publication DOI is invented.

New analysis software is [MIT-licensed](LICENSE). New author-created numerical
outputs, documentation and figures use [CC BY 4.0](LICENSE-DATA.md). Source data,
the alphaXiv/OpenResearch plotting helper, and any separate manuscript templates
retain their own terms and attribution; see
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) and
[`LICENSES/`](LICENSES/).
