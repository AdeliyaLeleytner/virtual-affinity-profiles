# Fixed-score panel methods provenance

Original cached fixed-score panel provenance; no new inference, no environment installation, no checkpoint download. Original runtime assertions are distinguished from recovered run records. Unknowns remain unknown.

| Channel | Original version / source | Raw field and analysis orientation | Reproducibility limit |
|---|---|---|---|
| boltz_affinity | boltz 2.2.1 | affinity_pred_value; multiply by -1, no clipping; higher is stronger | Pinned preflight assertions and runtime-proof hashes in export; full per-host runtime proof bodies unavailable in inspected local export |
| boltzina_vina_{zero,far}_affinity | boltz 2.2.1; Boltzina git 8c949647ccdd0fc5948f39004f1766f1e0bd9d73 | affinity_pred_value; multiply by -1, no clipping | Pinned preflight assertions; requirements snapshot records installed source dependency |
| nesso_affinity | nesso 1.0.0; git 6c72f66720d9d3447fd73c515cda963e39128b1f | affinity_pred_value; multiply by -1, no clipping; equals predicted pIC50 minus 6 algebraically, not a recalibration to the present assays | See retained metadata and hashes |
| raw_gnina_own | 1.3.3, gnina.cuda12.8.static | CNNaffinity; identity; higher predicts stronger affinity | rescore of Vina poses and independent GNINA search are separate workflows; GNINA own rank0 is highest CNNscore, not lowest energy |
| gnina | Per-job version link not recovered; project installed GNINA receipt is v1.3.3 | CNNaffinity; identity; higher predicts stronger affinity | CNNaffinity on original Vina pose0; pose-coordinate audit max fixed-frame RMSD=0; independent GNINA search is a separate workflow |
| vina | 2.1 | Vina energy; negative energy; higher is stronger | See retained metadata and hashes |
| raw_ad4 | v1.6 CUDA12 official release | reported AD4 estimated free energy of binding; negative energy; higher is stronger | See retained metadata and hashes |
| DOCKSTRING | Figshare dataset v1 doi:10.6084/m9.figshare.16511577.v1 | published source docking energy; negative energy; source variants preserved | original executable version not recovered in this audit |
| Docking44 | frozen project export2026-07-30; Nikitin etal Pharmaceutics2025 doi:10.3390/pharmaceutics17121573 | published source docking energy; negative energy; source variants preserved | original executable version not recovered in this audit |

Nesso original documentation defines its affinity scalar as predicted log10(IC50/µM). Negation is an orientation transform; it is not empirical calibration to the present experimental assays. The analysis keeps binary binder probabilities separate.

The JSON companion records checkpoint hashes, HF revisions, exact raw matrix hashes, relative source-code pointers, and copied evidence-file hashes. Original inference was not rerun.
