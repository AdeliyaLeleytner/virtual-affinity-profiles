# Source notices

The prepared arrays are selected/reorganized derivatives, not original publisher
archives. Missingness, quantitative censoring, source identities and source variants
are preserved. Primary metadata and source hashes are in `provenance/`.

- **SPD:** Jeffrey Sutherland, Dimitar Yonchev, Alexander Fekete and Laszlo Urban,
  *A novel preclinical secondary pharmacology resource illuminates target-adverse
  drug reaction associations of marketed drugs*. Original supplementary dataset:
  [Zenodo 8103950](https://doi.org/10.5281/zenodo.8103950), CC BY 4.0; associated
  [Nature Communications article](https://doi.org/10.1038/s41467-023-40064-9).
  `binding_records.csv.gz`, `human_all_records.csv.gz` and their assay/identity
  metadata are canonicalized derivatives. The SPD GitHub code licence is not the
  dataset licence.
- **DOCKSTRING:** retain attribution to the DOCKSTRING dataset authors listed in
  `provenance/source_licenses.json`; dataset
  [Figshare 16511577 v1](https://doi.org/10.6084/m9.figshare.16511577.v1), Apache 2.0.
  The licence is retained in `LICENSES/Apache-2.0.txt`. Associated paper:
  [10.1021/acs.jcim.1c01334](https://doi.org/10.1021/acs.jcim.1c01334).
  `datasets/dockstring.*` select and reorient source scores; duplicated source
  molecular variants retain identity weights. Negation changes orientation only.
- **Docking44:** Nikitin et al., [Pharmaceutics 2025, 17,
  1573](https://doi.org/10.3390/pharmaceutics17121573), CC BY 4.0 as recorded in the
  original project data manifest. `datasets/own44.*` derives from the frozen
  2026-07-30 project export `df_final_v4.csv.gz`, hash recorded in provenance.
  This export is not claimed byte-identical to a publisher supplement.
- **UniProt:** canonical sequences and annotations from the UniProt Consortium,
  retrieved 2026-10-22, [CC BY 4.0](https://www.uniprot.org/help/license).
  Canonical gene sequences do not certify the actual assay or docked constructs.
- **Cached model scores:** author-generated outputs from Boltz-2/Boltzina, Nesso-1,
  Vina-GPU, AutoDock-GPU and GNINA. Cite their sources and preserve the channel
  versions, transforms, lineage and known gaps in `provenance/METHODS_TABLE.md`
  and `provenance/method_provenance.json`. Software/checkpoint licences are not
  replaced by the new analysis licence. No upstream software or weights are
  redistributed here.
- **RDKit:** molecular fingerprints/descriptors were generated using RDKit; the
  saved features allow replay without installing RDKit. Structures retain their
  source terms. Source-preparation versions are in `provenance/prepared_inputs.json`.

No third-party paper PDFs, article prose, published figures, model weights or
molecular pose archives are bundled. New analyses and figures are covered by
`LICENSE-DATA.md`, not by a blanket override of source terms.

## OpenResearch plotting helper and manuscript templates

`figures/orx_figstyle.py` derives from the installed alphaXiv/OpenResearch skill
helper, MIT License, Copyright (c) 2026 alphaXiv. Its full licence is retained in
`LICENSES/OpenResearch-MIT.txt`; the upstream licence is available at
https://raw.githubusercontent.com/alphaXiv/OpenResearch/main/LICENSE.
The author adaptation changes the publication width for the Journal of
Cheminformatics figure layout (170 mm) and records those changes in figure
provenance. This attribution also applies to any matching local helper copy.

Springer Nature LaTeX class and bibliography-style files, when distributed with
the separate manuscript/submission materials, retain their upstream copyright
and licence notices. The analysis MIT licence does not replace those notices.
