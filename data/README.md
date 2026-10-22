# Prepared inputs

All 20 original input files are retained byte-for-byte and verified by
`python scripts/verify_bundle.py`. This README is an added explanation; it is not
part of the frozen original input manifest. The public configurations write fresh
signatures for replay/refits; original run signatures remain separate provenance.

- `molecules.csv`: 1,399 identities, including the original 930 SPD score identities;
  469 additional score-only identities do not acquire supervised labels.
- `binding_records.csv.gz`: 33,050 actual-assay human binding records, 4,247 exact
  and 28,803 strict censored, 929 identities, 105 assays and 74 genes. The numerical
  coordinate is `6 - log10(source IC50 in micromolar)`, not Kd/pKd. Qualifiers,
  original values, strict bounds and source IDs are retained.
- `human_all_records.csv.gz`: 54,268 original human records retained for provenance;
  functional measurements are not pooled into the binding analysis.
- `assays.csv`, `proteins.csv`: actual assay metadata and 132 canonical human gene
  sequences; repeated assay contexts remain separate. A canonical sequence does
  not qualify an assay or docked construct.
- `molecule_features.npz`: saved Morgan radius-3/chirality fingerprints, six
  descriptors and Tanimoto matrix. `protein_features.npz`: amino-acid/dipeptide
  spectra and kernel. RDKit is not required to load these saved representations.
- `datasets/dockstring.npz`: 890 source rows, 885 identities, 58 columns; source
  identity variants remain distinct and have normalized weights.
- `datasets/own44.npz`: 668 source identities, 44 columns, missing entries preserved.
- `datasets/spd.npz`: 12 archived score heads, 930 identities and 10 targets;
  eight primary affinity/energy channels are explicitly named in the contract.
- `datasets/*_rows.csv`, `*_columns.csv`, `broad_columns.csv`: source identity and
  receptor/gene/context mappings. `splits/`: inherited molecular folds and native
  target-pair ordering. Unknown measurements are not negative labels.

NPZ loading uses `allow_pickle=False`. CSV scientific floats should be read with
`float_precision="round_trip"`. Input derivatives retain their original source
terms; see `../THIRD_PARTY_NOTICES.md`. These saved inputs support downstream
replay/refitting, not an assertion of end-to-end original GPU inference replay.
