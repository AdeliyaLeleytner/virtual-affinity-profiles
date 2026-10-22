# Reproducibility release v1.0.0

The authoritative scientific source analyses and fixed prediction replay are
complete. The repository contains code, prepared data, compact results and
figures for **When do virtual affinity profiles improve biological decisions?**
The manuscript and submission documents are separate local materials.

## Included completed material

- Core, all-binding, supplementary both-cold, canonical growth and decoder source
  predictions: five checksum-pinned downloadable archives.
- Fixed replay and exact/censor-involving strata: completed ORX run
  `c133ae0c-9f6a-4da6-a321-7ed9f79e585a`, commit `5807705`, exit 0;24 datasets
  (22 primary and 2 decoder sources) and 72 measurement-stratum source tables.
- Compact authoritative fixed/growth results under `results/`, with source/public
  hashes and explicit metadata derivatives. Raw large censuses are regenerated.
- Six main figures, two supplements and the final graphical abstract, cached
  scientific inputs, original renderers, image-source attribution and pinned
  figure environment. All eight numerical figures are170 mm wide.
- Exact frozen `src/`, `growth_src/`, `decoder_src/` packages and separately
  patched `replay_src/`; original fitting hashes remain unchanged.

## Verification and its limits

All five source archives were previously restored and checked in a fresh relocated
repository, including scientific/module hashes, geometry and the decoder's exact
inheritance check. That copy was removed; the verified archives remain unchanged.
See `provenance/relocation_restore_qa.json`.

The final bounded public check used a new relocated code/data copy with installed
pinned analysis dependencies. Source hashes,66 existing technical/synthetic tests
and deterministic geometry passed. It then restored only the small fingerprint
archive and used the actual public CLI/config to replay every saved split of one
representative Docking44 all-binding dataset, both transfer regimes, 400 draws,
and its exact/censor strata. All32 compact point/quantile/support tables matched
the authoritative run02 values within the stated numerical rounding tolerance.
Peak additional temporary storage was178,964,301 bytes (below 500MiB); the temporary
copy was removed. No scientific-study fit or inference was performed.
See `provenance/final_bounded_public_qa.json`.

The final sanitized public figure bundle passed 227 relocation predicates in its
separate validated Python 3.12 environment with RDKit/SciPy imports blocked. All
normal renderer commands reproduced nine PNGs and captions exactly; required
input hashes were unchanged. This establishes relocation within that environment,
not pixel identity across library/OS versions.

A previous duplicate public **full** replay encountered storage exhaustion. Its
44 completed growth-table comparisons and partial outputs are historical checks,
not a complete public full-replay claim. The authoritative analyses are complete;
the bounded public checks do not replace or claim another full-study replay.
The exact limitation remains in
`provenance/disk_pressure_verification_receipt.json`.

The source-archive checks, bounded checks and figure checks did not test network
asset transfer. Versioned addresses and checksums are in the archive catalog.
No article DOI, original GPU-inference reconstruction or universal absence of
useful profile information is claimed.
