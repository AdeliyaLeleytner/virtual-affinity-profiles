"""Recompute published descriptive geometry from saved scientific score arrays.

No assay outcomes, fitting, docking, GPU inference, clipping or imputation.
Original mathematical transforms are documented in GEOMETRY_REPLAY.md.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity_weights(ids):
    _, inverse, counts = np.unique(ids, return_inverse=True, return_counts=True)
    weights = 1.0 / counts[inverse]
    return weights / weights.sum()


def corr(a, b, weights):
    a = a - weights @ a
    b = b - weights @ b
    denominator = np.sqrt((weights @ (a * a)) * (weights @ (b * b)))
    return float(weights @ (a * b) / denominator) if denominator > 0 else np.nan


def standardize(matrix, weights):
    centered = matrix - weights @ matrix
    sd = np.sqrt(weights @ (centered * centered))
    if np.any(sd <= 0):
        raise ValueError('Constant column needs explicit handling')
    return centered / sd


def transform(matrix, view):
    if view.startswith('rank'):
        matrix = rankdata(matrix, axis=0)
    if view.endswith('centered'):
        matrix = matrix - matrix.mean(axis=1, keepdims=True)
    return matrix


def svd_info(z, weights):
    x = np.sqrt(weights[:, None]) * z
    _, singular, _ = np.linalg.svd(x, full_matrices=False)
    tolerance = max(x.shape) * np.finfo(x.dtype).eps * singular[0]
    fractions = singular ** 2 / np.sum(singular ** 2)
    return singular, fractions, int(np.sum(singular > tolerance)), float(tolerance)


def join_broad(inputs):
    # Read score arrays directly. The supervised loader can restrict anchor eligibility.
    with np.load(inputs / 'datasets/dockstring.npz', allow_pickle=False) as a, np.load(inputs / 'datasets/own44.npz', allow_pickle=False) as b:
        if len(set(b['molecule_ids'])) != len(b['molecule_ids']):
            raise ValueError('Docking44 identity map is no longer one-to-one')
        lookup = {str(identity): i for i, identity in enumerate(b['molecule_ids'])}
        left = np.flatnonzero(np.isin(a['molecule_ids'], b['molecule_ids']))
        right = np.array([lookup[str(identity)] for identity in a['molecule_ids'][left]])
        return dict(matrix=np.c_[a['scores'][left], b['scores'][right]],
                    ids=a['molecule_ids'][left].copy(),
                    source_rows=a['source_row_indices'][left].copy(),
                    own44_source_rows=b['source_row_indices'][right].copy(),
                    column_ids=np.r_[a['column_ids'], b['column_ids']],
                    genes=np.r_[a['column_genes'], b['column_genes']],
                    blocks=np.array(['DOCKSTRING'] * a['scores'].shape[1] + ['Docking44'] * b['scores'].shape[1]))


def broad_geometry(inputs, output):
    joined = join_broad(inputs)
    raw, ids = joined['matrix'], joined['ids']
    complete = np.isfinite(raw).all(axis=1)
    matrix, complete_ids = raw[complete], ids[complete]
    weights = identity_weights(complete_ids)
    sources = ['DOCKSTRING', 'Docking44']
    indices = {name: np.flatnonzero(joined['blocks'] == name) for name in sources}
    pd.DataFrame(dict(standard_molecule_id=ids, dockstring_source_row=joined['source_rows'],
                      complete_102=complete, n_observed_scores=np.isfinite(raw).sum(axis=1),
                      identity_weight_all=identity_weights(ids))).to_csv(output / 'geometry_support.csv', index=False)
    pd.DataFrame(dict(standard_molecule_id=ids, dockstring_source_row=joined['source_rows'],
                      own44_source_row=joined['own44_source_rows'])).to_csv(output / 'joint_source_identity_map.csv', index=False)
    summaries, spectra, shared, variance = [], [], [], []
    for source, ix in indices.items():
        block = matrix[:, ix]
        grand = float(np.mean(weights @ block))
        molecule = block.mean(axis=1) - grand
        target = weights @ block - grand
        interaction = block - grand - molecule[:, None] - target[None, :]
        total = float(np.sum(weights[:, None] * (block - grand) ** 2))
        components = [len(ix) * (weights @ (molecule * molecule)), np.sum(target * target), np.sum(weights[:, None] * interaction ** 2)]
        np.testing.assert_allclose(sum(components), total, rtol=1e-12, atol=1e-12)
        variance.append(dict(source=source, n_molecules=len(complete_ids), n_targets=len(ix),
                             molecule_fraction=components[0] / total, target_fraction=components[1] / total,
                             interaction_fraction=components[2] / total, total_sum_squares=total))
    shared_genes = sorted(set(joined['genes'][indices[sources[0]]]) & set(joined['genes'][indices[sources[1]]]))
    for view in ['raw', 'row_centered', 'rank_columns', 'rank_centered']:
        blocks = {name: standardize(transform(matrix[:, ix], view), weights) for name, ix in indices.items()}
        joint = np.concatenate([blocks[name] for name in sources], axis=1)
        equal = np.concatenate([blocks[name] * np.sqrt(len(joined['column_ids']) / (2 * len(indices[name]))) for name in sources], axis=1)
        for name, z in list(blocks.items()) + [('Joint_equal_columns', joint), ('Joint_equal_blocks', equal)]:
            singular, fractions, rank, tolerance = svd_info(z, weights)
            summaries.append(dict(view=view, panel=name, n_rows=len(z), n_identities=len(set(complete_ids)),
                                  n_columns=z.shape[1], pr=float(1 / np.sum(fractions ** 2)), pc1=float(fractions[0]),
                                  numerical_rank=rank, numerical_rank_tolerance=tolerance))
            for k, (value, fraction) in enumerate(zip(singular, fractions), 1):
                spectra.append(dict(view=view, panel=name, component=k, eigenvalue=float(value * value),
                                    fraction=float(fraction), cumulative=float(fractions[:k].sum())))
        correlations = joint.T @ (weights[:, None] * joint)
        pd.DataFrame(correlations, index=joined['column_ids'], columns=joined['column_ids']).to_csv(output / ('correlation_' + view + '_complete.csv'))
        for gene in shared_genes:
            a = np.flatnonzero((joined['genes'] == gene) & (joined['blocks'] == sources[0]))[0]
            b = np.flatnonzero((joined['genes'] == gene) & (joined['blocks'] == sources[1]))[0]
            shared.append(dict(view=view, support='all102_complete', gene=gene,
                               column_a=joined['column_ids'][a], column_b=joined['column_ids'][b],
                               n_rows=len(matrix), n_identities=len(set(complete_ids)), pearson=float(correlations[a, b]),
                               spearman=corr(rankdata(joint[:, a]), rankdata(joint[:, b]), weights)))
    for gene in shared_genes:
        a = np.flatnonzero((joined['genes'] == gene) & (joined['blocks'] == sources[0]))[0]
        b = np.flatnonzero((joined['genes'] == gene) & (joined['blocks'] == sources[1]))[0]
        available = np.isfinite(raw[:, a]) & np.isfinite(raw[:, b])
        pair_weights = identity_weights(ids[available])
        shared.append(dict(view='raw', support='pair_available', gene=gene,
                           column_a=joined['column_ids'][a], column_b=joined['column_ids'][b],
                           n_rows=int(available.sum()), n_identities=len(set(ids[available])),
                           pearson=corr(raw[available, a], raw[available, b], pair_weights),
                           spearman=corr(rankdata(raw[available, a]), rankdata(raw[available, b]), pair_weights)))
    shared_frame = pd.DataFrame(shared)
    shared_frame.to_csv(output / 'shared_gene_comparisons.csv', index=False)
    shared_frame[(shared_frame.support == 'all102_complete') & shared_frame['view'].isin(['raw', 'row_centered'])].to_csv(output / 'shared_gene_main_figure.csv', index=False)
    pd.DataFrame(summaries).to_csv(output / 'spectral_summary.csv', index=False)
    pd.DataFrame(spectra).to_csv(output / 'full_spectra.csv', index=False)
    pd.DataFrame(variance).to_csv(output / 'variance_components.csv', index=False)
    means = [matrix[:, indices[name]].mean(axis=1) for name in sources]
    pd.DataFrame(dict(standard_molecule_id=complete_ids, **dict(zip(sources, means)))).to_csv(output / 'shared_means.csv', index=False)
    pearson = corr(means[0], means[1], weights)
    pd.DataFrame([[1., pearson], [pearson, 1.]], index=sources, columns=sources).to_csv(output / 'shared_mean_correlations.csv')
    agreement = dict(n_rows=len(matrix), n_identities=len(set(complete_ids)), pearson=pearson,
                     spearman=corr(rankdata(means[0]), rankdata(means[1]), weights))
    (output / 'common_mean_agreement.json').write_text(json.dumps(agreement, indent=2) + '\n')
    return dict(joined_rows=len(raw), joined_identities=len(set(ids)), columns=raw.shape[1],
                complete_rows=len(matrix), complete_identities=len(set(complete_ids)), shared_genes=shared_genes,
                mean_agreement=agreement, joint_pr={row['view']: row['pr'] for row in summaries if row['panel'] == 'Joint_equal_columns'})


def composition_and_support(inputs, output):
    metadata = pd.read_csv(inputs / 'datasets/broad_columns.csv', float_precision='round_trip')
    metadata.to_csv(output / 'column_metadata.csv', index=False)
    metadata.groupby(['broad_class', 'source']).size().rename('n_columns').reset_index().to_csv(output / 'panel_composition.csv', index=False)
    molecules = pd.read_csv(inputs / 'molecules.csv', float_precision='round_trip')
    # Only record identity membership is used; no IC50 value, bound or outcome is loaded.
    label_ids = set(pd.read_csv(inputs / 'binding_records.csv.gz', usecols=['molecule_index']).molecule_index)
    rows = []
    for name, flag in [('DOCKSTRING', 'in_dockstring'), ('Docking44', 'in_own44'), ('Shared sources', 'in_joint')]:
        use = molecules[flag]
        rows.append(dict(cohort=name, score_identities=int(use.sum()), binding_anchors=int((use & molecules.molecule_index.isin(label_ids)).sum())))
    pd.DataFrame(rows).to_csv(output / 'cohort_support.csv', index=False)


def spd_spectrum(matrix):
    # The original SPD implementation uses sample covariance then converts it to correlation.
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    covariance = centered.T @ centered / (len(centered) - 1)
    sd = np.sqrt(np.diag(covariance))
    if np.any(sd == 0):
        raise ValueError('Constant target: correlation undefined')
    eigenvalues, vectors = np.linalg.eigh(covariance / np.outer(sd, sd))
    eigenvalues, vectors = eigenvalues[::-1], vectors[:, ::-1]
    if np.min(eigenvalues) < -1e-9 * max(1., abs(eigenvalues[0])):
        raise ValueError('Non-PSD numerical correlation matrix')
    eigenvalues = np.maximum(eigenvalues, 0)  # Floating-point negative roundoff only.
    fractions = eigenvalues / eigenvalues.sum()
    axis = vectors[:, 0]
    if axis.sum() < 0:
        axis = -axis
    pc = (centered / centered.std(axis=0, ddof=1)) @ axis
    rowmean = matrix.mean(axis=1)
    return dict(pr=float(1 / (fractions @ fractions)), pc1=float(fractions[0]),
                pc2_cumulative=float(fractions[:2].sum()),
                entropy_rank=float(np.exp(-np.sum(fractions[fractions > 0] * np.log(fractions[fractions > 0])))),
                pc1_uniform_cosine=float(axis.sum() / np.sqrt(len(axis))),
                pc1_positive_loadings=int((axis > 0).sum()),
                pc1_rowmean_correlation=float(np.corrcoef(pc, rowmean)[0, 1]) if rowmean.std() > 1e-12 else None)


def spd_geometry(inputs, output):
    with np.load(inputs / 'datasets/spd.npz', allow_pickle=False) as z:
        scores = z['scores']
        mask = z['common_complete_rows_all12']
        np.testing.assert_array_equal(mask, np.isfinite(scores).all(axis=(0, 2)))
        identifiers = pd.read_csv(ROOT / 'geometry_inputs/spd_common_all12_ids.csv').molecule_id.to_numpy(str)
        np.testing.assert_array_equal(z['molecule_ids'][mask], identifiers)
        records = []
        for k in np.flatnonzero(z['primary_affinity_mask']):
            matrix = scores[k, mask]
            for mode in ['raw', 'molecule_centered']:
                values = matrix if mode == 'raw' else matrix - matrix.mean(axis=1, keepdims=True)
                records.append(dict(channel=str(z['channel_names'][k]), channel_type=str(z['channel_types'][k]),
                                    variant='unclipped', support='common_all12', mode=mode, spectrum='correlation',
                                    n_molecules=len(matrix), n_targets=matrix.shape[1], **spd_spectrum(values)))
    pd.DataFrame(records).to_csv(output / 'spd_readouts.csv', index=False)
    return dict(source_heads=scores.shape[0], reported_primary_heads=len(records) // 2,
                common_all12_identities=len(identifiers), explicit_identity_list_sha256=sha(ROOT / 'geometry_inputs/spd_common_all12_ids.csv'))


def verify_references(output, reference):
    # Numerical equality tolerance only; no scientific success threshold or hypothesis test.
    results = []
    expected_names = set(json.loads((ROOT / 'provenance/geometry_sources.json').read_text())['reference_outputs'])
    actual_names = {p.name for p in reference.glob('*.csv')}
    if actual_names != expected_names:
        raise AssertionError('Missing or unexpected geometry reference files: ' + str(actual_names ^ expected_names))
    for expected in sorted(reference.glob('*.csv')):
        actual = output / expected.name
        a = pd.read_csv(actual, float_precision='round_trip')
        b = pd.read_csv(expected, float_precision='round_trip')
        if list(a.columns) != list(b.columns) or a.shape != b.shape:
            raise AssertionError('Reference schema/support mismatch: ' + expected.name)
        maximum = 0.
        for column in a:
            if pd.api.types.is_numeric_dtype(a[column]) and pd.api.types.is_numeric_dtype(b[column]):
                av, bv = a[column].to_numpy(float), b[column].to_numpy(float)
                np.testing.assert_allclose(av, bv, rtol=2e-10, atol=2e-10, equal_nan=True, err_msg=expected.name + ':' + column)
                finite = np.isfinite(av) & np.isfinite(bv)
                if finite.any():
                    maximum = max(maximum, float(np.max(np.abs(av[finite] - bv[finite]))))
            else:
                np.testing.assert_array_equal(a[column].fillna('').astype(str), b[column].fillna('').astype(str), err_msg=expected.name + ':' + column)
        results.append(dict(file=expected.name, rows=len(a), maximum_absolute_numeric_difference=maximum,
                            reference_sha256=sha(expected), regenerated_sha256=sha(actual)))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=ROOT / 'data')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'outputs/geometry/data')
    parser.add_argument('--reference-dir', type=Path, default=ROOT / 'reference/geometry')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=1):
        broad = broad_geometry(args.input_dir, args.output_dir)
        spd = spd_geometry(args.input_dir, args.output_dir)
        composition_and_support(args.input_dir, args.output_dir)
        checks = verify_references(args.output_dir, args.reference_dir)
    receipt = dict(status='replayed and compared', python=platform.python_version(), numpy=np.__version__,
                   broad=broad, spd=spd, reference_checks=checks, numerical_comparison=dict(atol=2e-10, rtol=2e-10),
                   code_sha256=sha(Path(__file__)),
                   input_sha256={str(p.relative_to(args.input_dir)): sha(p) for p in [args.input_dir / 'datasets/dockstring.npz', args.input_dir / 'datasets/own44.npz', args.input_dir / 'datasets/spd.npz', args.input_dir / 'datasets/broad_columns.csv', args.input_dir / 'molecules.csv', args.input_dir / 'binding_records.csv.gz']},
                   source_provenance='provenance/geometry_sources.json', experimental_outcome_values_read=False, record_identity_membership_used_for_support=True,
                   new_model_fit=False, new_inference=False, scientific_pass_fail_threshold=None,
                   geometry_scope='Same saved score support/transforms as source geometry figures; no family-label permutations, experimental associations or biological mechanism claims are added.')
    (args.output_dir / 'geometry_replay_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], broad=broad, spd=spd, reference_files_checked=len(checks)), indent=2))


if __name__ == '__main__':
    main()
