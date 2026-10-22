# Source-derived figure assets

All scientific plots are original renders from the accompanying frozen CSV tables. No published paper figure, RCSB gallery image, stock image or generated image was incorporated. Molecular pictures are original RDKit depictions; protein pictures are original coordinate renders. Software licenses do not replace the upstream dataset licenses documented by the project.

## Molecules and reference structures

The graphical abstract depicts butenafine and indomethacin from the exact prepared reference SMILES in `data/molecule_illustrations.csv`. Their actual illustrative Docking44 profile values are in `data/illustration_profiles.csv`. Figure 6 depicts fluphenazine, trifluoperazine, amoxapine and loxapine from `data/figure6_choices/selected_cases.csv`, preserving standardized molecular IDs, names and reference SMILES. The full source representations and variant weights are separately retained in `all_source_variant_decisions.csv`; images do not replace those source records.

`assets/figure6_choices/molecular_depiction_receipt.json` records the RDKit version, reference SMILES and image hashes. Normal rendering reads the cached PNGs and preserves the SVG originals. Optional `build_choice_examples.py --mode render --regenerate-molecules` requires the separate RDKit environment. RDKit software uses a BSD license; its authors recommend citing “RDKit: Open-source cheminformatics” and the version used. See the [official RDKit overview and citation guidance](https://www.rdkit.org/docs/Overview.html). The graphics were created with RDKit 2025.03.6; this is a software attribution, not an attempt to relicense chemical-source data.

## Protein coordinate illustrations

The final graphical abstract uses source-deposited coordinate examples, not new structure predictions or a certification of experimental assay constructs. `assets/structure_retrieval.json` preserves requested/observed RCSB URLs, UTC retrieval and SHA-256 hashes. `graphical_abstract_render_receipt.json` records the first-model longest-chain choices and C-alpha atom counts. Source mapping is in `data/protein_illustrations.csv`. Traces are generated separately across coordinate gaps, oriented deterministically and smoothed for illustration; engineered/fusion residues may be present.

| Source ID and illustrated label | Coordinate record | Original structure publication |
| --- | --- | --- |
| 6CM4, DRD2 | [PDB 6CM4](https://www.rcsb.org/structure/6CM4) | Wang et al. (2018), *Structure of the D2 dopamine receptor bound to the atypical antipsychotic drug risperidone*. [DOI](https://doi.org/10.1038/nature25758) |
| 6PDJ, LCK | [PDB 6PDJ](https://www.rcsb.org/structure/6PDJ) | O'Malley et al. (2019), *Discovery of Pyridazinone and Pyrazolo[1,5-a]pyridine Inhibitors of C-Terminal Src Kinase*. [DOI](https://doi.org/10.1021/acsmedchemlett.9b00354). The deposited PDB TITLE/COMPND/SOURCE records explicitly identify the illustrated structure as LCK. |
| 3B66, AR | [PDB 3B66](https://www.rcsb.org/structure/3B66) | Bohl et al. (2008), *Effect of B-ring substitution pattern on binding mode of propionamide selective androgen receptor modulators*. [DOI](https://doi.org/10.1016/j.bmcl.2008.09.002) |

PDB archive data are available under CC0 1.0; wwPDB encourages attribution to the original structural-data authors. This was checked against the [official wwPDB usage policy](https://www.wwpdb.org/about/usage-policies) on 22 October 2026. The table provides that attribution. These coordinate-derived pictures are not copied publication artwork. Optional `regenerate_illustration_assets.py` uses the retained PDB files and reference SMILES with the separately pinned RDKit/SciPy environment.

## Figure style and environment

The vendored `orx_figstyle.py` retains the OpenResearch MIT attribution supplied with the public release. Project modifications select bundled DejaVu Sans and the Journal of Cheminformatics 85/170 mm widths. Preserve its original copyright/license notice. Springer class/BST copyrights are independent and are not altered by this figure workflow.

Normal figures were validated with Python 3.12.11, Matplotlib 3.10.6, NumPy 1.26.4, pandas 2.3.2 and Pillow 12.3.0, separately from the Python 3.9 scientific analysis environment. Optional asset regeneration additionally used RDKit 2025.03.6 and SciPy 1.16.2. Portable pixel comparisons test relocation within this validated environment; they do not assert pixel identity across Matplotlib, font-library, image-library or operating-system versions.
