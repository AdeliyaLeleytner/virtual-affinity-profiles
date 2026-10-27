"""Panel geometry, recomputed from the frozen score arrays.

Reads data/datasets/*.npz from this repository and reproduces the effective
dimensionality and first-component fractions reported in the manuscript. Run it
on its own to print the table and the check against the published values.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATASETS = os.path.join(HERE, os.pardir, "data", "datasets")

NAME = {
    "vina": "Vina",
    "raw_ad4": "AutoDock4",
    "gnina": "GNINA (Vina pose)",
    "raw_gnina_own": "GNINA (own search)",
    "boltz_affinity": "Boltz-2",
    "boltzina_vina_far_affinity": "Boltzina (far)",
    "boltzina_vina_zero_affinity": "Boltzina (zero)",
    "nesso_affinity": "Nesso-1",
}


def spectrum(Z, w):
    """Effective dimensionality and first-component fraction of a weighted matrix.

    Columns are standardised under the row weights, then the singular values of
    W^(1/2) Z give eigenvalues lambda_j. With f_j = lambda_j / sum(lambda), the
    effective number of directions is 1 / sum(f_j^2), also called the
    participation ratio.
    """
    w = w / w.sum()
    mu = (w[:, None] * Z).sum(0)
    sd = np.sqrt((w[:, None] * (Z - mu) ** 2).sum(0))
    sd[sd == 0] = 1.0
    Zs = (Z - mu) / sd
    s = np.linalg.svd(np.sqrt(w)[:, None] * Zs, compute_uv=False)
    f = s ** 2 / (s ** 2).sum()
    return 1.0 / (f ** 2).sum(), f[0]


def spd_table():
    """{channel: (effective dimensionality, first-component fraction)} and n."""
    d = np.load(os.path.join(DATASETS, "spd.npz"), allow_pickle=True)
    ch = [str(c) for c in d["channel_names"]]
    rows = d["common_complete_rows_affinity8"].astype(bool)
    w = d["variant_weights"][rows]
    out = {}
    for name in NAME:
        if name not in ch:
            continue
        Z = d["scores"][ch.index(name)][rows]
        if np.isnan(Z).any():
            continue
        out[name] = spectrum(Z, w)
    return out, int(rows.sum())


def joint_broad():
    """Geometry of the two broad panels, jointly, on their shared molecules."""
    a = np.load(os.path.join(DATASETS, "dockstring.npz"), allow_pickle=True)
    b = np.load(os.path.join(DATASETS, "own44.npz"), allow_pickle=True)
    ida = {str(x): i for i, x in enumerate(a["molecule_ids"])}
    idb = {str(x): i for i, x in enumerate(b["molecule_ids"])}
    shared = sorted(set(ida) & set(idb))
    A = a["scores"][[ida[s] for s in shared]]
    B = b["scores"][[idb[s] for s in shared]]
    keep = ~(np.isnan(A).any(1) | np.isnan(B).any(1))
    A, B = A[keep], B[keep]
    w = np.ones(len(A))
    pr_raw, f_raw = spectrum(np.hstack([A, B]), w)
    centred = np.hstack([A - A.mean(1, keepdims=True), B - B.mean(1, keepdims=True)])
    pr_cen, f_cen = spectrum(centred, w)
    return dict(n=int(keep.sum()), cols=A.shape[1] + B.shape[1],
                pr_raw=pr_raw, f_raw=f_raw, pr_cen=pr_cen, f_cen=f_cen,
                r_means=float(np.corrcoef(A.mean(1), B.mean(1))[0, 1]))


if __name__ == "__main__":
    tab, n = spd_table()
    print(f"Ten-target panel, {n} molecules with complete scores\n")
    print(f"  {'score':<20}{'eff. dim.':>10}{'PC1 (%)':>10}")
    for k, (pr, f1) in sorted(tab.items(), key=lambda x: x[1][0]):
        print(f"  {NAME[k]:<20}{pr:>10.2f}{100 * f1:>10.1f}")
    prs = [v[0] for v in tab.values()]
    f1s = [v[1] for v in tab.values()]
    print(f"\n  recomputed: {min(prs):.2f}-{max(prs):.2f}, "
          f"PC1 {100 * min(f1s):.1f}-{100 * max(f1s):.1f}%")
    print(f"  manuscript: 1.43-2.93, PC1 55.1-83.3%")

    j = joint_broad()
    print(f"\nJoint broad panel: {j['n']} molecules x {j['cols']} columns")
    print(f"  raw       {j['pr_raw']:.2f}  PC1 {100 * j['f_raw']:.1f}%"
          f"   (manuscript 1.63 / 78.3%)")
    print(f"  centred   {j['pr_cen']:.2f}  PC1 {100 * j['f_cen']:.1f}%"
          f"   (manuscript 20.35 / 15.9%)")
    print(f"  r(molecular means) {j['r_means']:.3f}   (manuscript 0.973)")
