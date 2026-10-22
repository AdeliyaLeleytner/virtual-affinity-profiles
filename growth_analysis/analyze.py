"""Read complete saved growth outputs and replay endpoint predictions without fits.
Scientific execution is intended for fixed_study.run / a registered ORX node.
python -m growth_analysis.analyze --config growth_analysis/config.portable.json
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys
import time

import numpy as np
import pandas as pd
from .manifest import verify_manifest

ROOT = Path(__file__).resolve().parents[1]
PRIMARY_METRICS = ['D_direction_concordance', 'pareto_correct_fraction']
SUPPORT = ['n_assay_pairs', 'n_gene_pairs', 'n_rectangles', 'n_exact', 'n_identified_D', 'n_pareto']


def same_bits(a, b):
    """NaNs count as identical only when their actual binary representations match."""
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=str) + '\n')


class Sources:
    def __init__(self):
        self.paths = {}

    def add(self, path):
        path = Path(path).resolve()
        if path not in self.paths:
            self.paths[path] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    def csv(self, path):
        return pd.read_csv(self.add(path), float_precision='round_trip')

    def json(self, path):
        return json.loads(self.add(path).read_text())

    def npz(self, path):
        return np.load(self.add(path), allow_pickle=False)


def validate_curve_summary(frame, n_columns, orders):
    """Technical completeness and common-support checks, with no effect cutoff."""
    grouping = ['model', 'domain', 'weighting', 'comparison_support']
    for key, group in frame.groupby(grouping, dropna=False):
        if set(group.order) != set(orders):
            raise ValueError(f'Incomplete order set for {key}')
        for order, curve in group.groupby('order'):
            if curve.k.duplicated().any() or set(curve.k) != set(range(n_columns + 1)):
                raise ValueError(f'Incomplete or duplicate k values for {key}/{order}')
        for col in SUPPORT:
            if group[col].nunique(dropna=False) != 1:
                raise ValueError(f'Support changes across growth cells: {key}/{col}')
        zero = group[group.k.eq(0)]
        for metric in PRIMARY_METRICS:
            if zero[metric].nunique(dropna=False) != 1:
                raise ValueError(f'k0 metric changes across orders: {key}/{metric}')
    return {'all_k_complete': True, 'support_constant': True, 'n_summary_rows': len(frame)}


def replay_arrays(saved, shape):
    rr, aa = saved['test_rows'], saved['test_assays']
    pred = {}
    for name in saved.files:
        if name in ['test_rows', 'test_assays']:
            continue
        result = np.full(shape, np.nan)
        result[np.ix_(rr, aa)] = saved[name]
        pred[name] = result
    return rr, aa, pred


def selected_costs(data, split, diagnostics, score_acquisition, support):
    """Exact per-query unions, distinguishing deployment from label-based search."""
    cols = data['column_genes']
    declared = score_acquisition.get('chemistry_plus_selected_score', {}).get('training_candidate_columns')
    if declared is not None and int(declared) != len(cols):
        raise ValueError('Selected-score candidate acquisition metadata disagrees with the saved panel')
    for row in diagnostics['chemistry_plus_selected_score']['column_selection']:
        if 'candidate_columns' in row and int(row['candidate_columns']) != len(cols):
            raise ValueError('Per-assay selected-score diagnostics disagree with the full candidate panel')
    matched = {j: int(np.flatnonzero(cols == gene)[0]) if np.any(cols == gene) else -1
               for j, gene in enumerate(data['assays'].gene)}
    selected = {int(r['assay_index']): int(r['column_index']) for r in
                diagnostics['chemistry_plus_selected_score']['column_selection']}
    choices = {'chemistry_plus_target_score': matched, 'chemistry_plus_selected_score': selected}
    rr, aa = np.array(split['train_rows']), np.array(split['train_assays'])
    labelled = [int(j) for j in aa if (data['weights'][rr, j] > 0).any()]
    evaluated = list(map(int, split['test_assays']))
    defined_pairs = support[support.n_saved_rectangles.gt(0)]
    contrast_evaluated = sorted(set(defined_pairs.assay_index_a.astype(int)) | set(defined_pairs.assay_index_b.astype(int)))
    assay_rows, pair_rows, unions = [], [], []
    for model, mapping in choices.items():
        union = sorted({mapping[j] for j in evaluated if mapping.get(j, -1) >= 0})
        train_union = sorted({mapping[j] for j in labelled if mapping.get(j, -1) >= 0})
        candidate_count = (len(cols) if model == 'chemistry_plus_selected_score' else len(train_union))
        metadata = dict(split=split['id'], model=model,
                        training_score_candidate_columns=candidate_count,
                        deployment_whole_profile_union_count=len(union),
                        deployment_whole_profile_union_columns=json.dumps(union),
                        deployment_evaluated_contrast_union_count=len({mapping[j] for j in contrast_evaluated if mapping.get(j, -1) >= 0}),
                        training_selected_union_count=len(train_union),
                        training_search_not_free=(model == 'chemistry_plus_selected_score'),
                        full_panel_mean_used=False)
        unions.append(metadata)
        for j in evaluated:
            col = mapping.get(j, -1)
            assay_rows.append(dict(**metadata, assay_index=j, assay_id=data['assays'].iloc[j].assay_id,
                                   gene=data['assays'].iloc[j].gene, query_column=col,
                                   query_column_gene=cols[col] if col >= 0 else '',
                                   deployment_query_column_count=int(col >= 0)))
        for row in support.itertuples():
            queried = sorted({mapping.get(int(j), -1) for j in [row.assay_index_a, row.assay_index_b]} - {-1})
            pair_rows.append(dict(**metadata, assay_index_a=row.assay_index_a,
                                  assay_index_b=row.assay_index_b, assay_id_a=row.assay_id_a,
                                  assay_id_b=row.assay_id_b, gene_a=row.gene_a, gene_b=row.gene_b,
                                  deployment_query_column_count=len(queried),
                                  deployment_query_columns=json.dumps(queried),
                                  n_saved_rectangles=row.n_saved_rectangles))
    return assay_rows, pair_rows, unions


def run(config, config_path=None):
    started = time.monotonic()
    options = config.get('growth_analysis', config)
    source = Sources()
    resolve = lambda p: (ROOT / p).resolve()
    output = resolve(options['output_dir'])
    output.mkdir(parents=True, exist_ok=True)
    pinned_manifest = resolve(options['source_manifest'])
    before_verification = verify_manifest(pinned_manifest)
    print(json.dumps(dict(source_verification='before', **before_verification)), flush=True)
    code_root = resolve(options['frozen_code_dir'])
    frozen_data = load_module('growth_saved_data', source.add(code_root / 'data.py'))
    frozen_evaluation = load_module('growth_saved_evaluation', source.add(code_root / 'evaluation.py'))
    sys.path.insert(0, str(resolve(options['postprocessing_code_dir'])))
    from fast_resampling import summary
    for filename in ['fast_resampling.py', 'postprocess.py']:
        source.add(resolve(options['postprocessing_code_dir']) / filename)
    from metric_support import audit_frame
    all_curves, all_assays, receipts = [], [], []
    folds = options.get('expected_folds', 5)
    orders = [f'{policy}_{seed}' for policy in ['random', 'family_stratified']
              for seed in options['seeds']]
    for spec in options['datasets']:
        folder = resolve(spec['growth_dir'])
        target = output / spec['dataset']
        target.mkdir(exist_ok=True)
        data = frozen_data.load(resolve(options['input_dir']), spec['dataset'],
                                scope=spec.get('scope', 'native_genes'), policy=spec.get('policy', 'all_records'))
        n_columns = data['scores'].shape[1]
        # Hash the exact prepared files read by data.load.
        for filename in ['molecules.csv', 'assays.csv', 'binding_records.csv.gz', 'proteins.csv',
                         'molecule_features.npz', 'protein_features.npz', f'datasets/{spec["dataset"]}.npz']:
            source.add(resolve(options['input_dir']) / filename)
        curve = source.csv(folder / 'summary_preferences.csv')
        checks = validate_curve_summary(curve, n_columns, orders)
        assay_curve = source.csv(folder / 'summary_assays.csv')
        curve['source_label'] = spec['label']
        all_curves.append(curve)
        assay_curve['source_label'] = spec['label']
        all_assays.append(assay_curve)
        comparator_summary = source.csv(folder / 'query_score_comparator_summary.csv')
        comparator_summary.to_csv(target / 'query_comparator_saved_summary.csv', index=False)
        costs, assay_costs, pair_costs, union_costs, details, all_support = [], [], [], [], [], []
        endpoint_checks = []
        class_table = source.csv(resolve(options['broad_columns_csv']))
        source_name = {'dockstring': 'DOCKSTRING', 'own44': 'Docking44'}[spec['dataset']]
        class_map = class_table[class_table.source.eq(source_name)].set_index('source_column_index').broad_class.to_dict()
        for fold in range(folds):
            split_dir = folder / f'mol_fold_{fold}'
            # Verify every completed cell and inspect its declared acquisition columns.
            canonical = None
            endpoint = None
            for order in orders:
                ordering = source.json(split_dir / order / 'order.json')
                columns = ordering['columns']
                if sorted(columns) != list(range(n_columns)):
                    raise ValueError(f'Invalid acquired-column permutation: {split_dir}/{order}')
                for k in range(n_columns + 1):
                    cell = split_dir / order / f'k{k}'
                    complete = source.json(cell / 'complete.json')
                    declared = complete.get('metadata', complete)
                    if declared.get('status') != 'complete' or declared.get('purchased_columns') != columns[:k]:
                        raise ValueError(f'Completion metadata disagrees with acquired columns: {cell}')
                    if not (cell / 'heldout_predictions.npz').is_file():
                        raise ValueError(f'Missing saved predictions: {cell}')
                    # complete.json can wrap the run metadata; column costs use order.json,
                    # not an inferred budget from a model feature count.
                    for model in ['global_chemistry', 'full_profile']:
                        costs.append(dict(dataset=spec['dataset'], split=f'mol_fold_{fold}',
                                          order=order, k=k, model=model,
                                          acquired_column_count=k, acquired_columns=json.dumps(columns[:k]),
                                          acquired_broad_class_count=len({class_map[c] for c in columns[:k]}),
                                          acquired_broad_classes=json.dumps(sorted({class_map[c] for c in columns[:k]})),
                                          score_free=(k == 0), uses_unacquired_full_panel_mean=False,
                                          completion_file=str(cell / 'complete.json')))
                    if k not in [0, n_columns]:
                        continue
                    saved = source.npz(cell / 'heldout_predictions.npz')
                    if k == 0:
                        if not same_bits(saved['global_chemistry'], saved['full_profile']):
                            raise ValueError(f'k0 arms differ bitwise: {cell}')
                        if canonical is None:
                            canonical = {name: saved[name].copy() for name in saved.files}
                        elif not all(same_bits(saved[name], canonical[name]) for name in saved.files):
                            raise ValueError(f'k0 order copies differ bitwise: {cell}')
                    else:
                        if endpoint is None:
                            endpoint = {name: saved[name].copy() for name in saved.files}
                        for model in ['global_chemistry', 'full_profile']:
                            finite = np.isfinite(saved[model]) & np.isfinite(endpoint[model])
                            endpoint_checks.append(dict(split=f'mol_fold_{fold}', order=order, model=model,
                                input_column_set_identical=True, predictions_bitwise_identical=same_bits(saved[model], endpoint[model]),
                                finite_support_identical=np.array_equal(np.isfinite(saved[model]), np.isfinite(endpoint[model])),
                                max_absolute_prediction_difference=float(np.max(np.abs(saved[model][finite] - endpoint[model][finite]))) if finite.any() else None,
                                interpretation='same-input-set technical numerical check; never composition variability'))
            comp_dir = split_dir / 'query_score_comparators'
            saved_comp = source.npz(comp_dir / 'heldout_predictions.npz')
            for saved in [saved_comp, endpoint]:
                for name in ['test_rows', 'test_assays']:
                    if not same_bits(saved[name], canonical[name]):
                        raise ValueError(f'Comparator indices do not match k0: {comp_dir}/{name}')
            split = source.json(comp_dir / 'split.json')
            rr, aa = canonical['test_rows'], canonical['test_assays']
            predictions = {'global_chemistry': canonical['global_chemistry']}
            predictions.update({name: saved_comp[name] for name in saved_comp.files if name not in ['test_rows', 'test_assays']})
            predictions.update({'full_panel_profile': endpoint['full_profile'],
                                'full_panel_level_chemistry': endpoint['global_chemistry']})
            arrays = {}
            for name, values in predictions.items():
                a = np.full(data['lower'].shape, np.nan)
                a[np.ix_(rr, aa)] = values
                arrays[name] = a
            metadata = dict(dataset=spec['dataset'], channel=data['channel'], scope=data['scope'],
                            policy=data['policy'], regime='new_molecules_known_proteins',
                            split=split['id'], variant_set='growth_zero_score_paired_replay')
            frames = frozen_evaluation.evaluate(data, arrays, rr, aa, metadata,
                                                compact=False, include_records=False)
            frames[2].to_csv(target / f'{split["id"]}_paired_rectangles.csv.gz', index=False, compression='gzip')
            frames[1].to_csv(target / f'{split["id"]}_paired_assay_metrics.csv.gz', index=False, compression='gzip')
            details.append(frames[2])
            all_support.append(frames[3])
            diag = source.json(comp_dir / 'fit_diagnostics.json')
            acquisition = source.json(comp_dir / 'score_acquisition.json')
            a, p, u = selected_costs(data, split, diag, acquisition, frames[3])
            assay_costs.extend(a);pair_costs.extend(p);union_costs.extend(u)
        detail = pd.concat(details, ignore_index=True)
        audit_frame(detail).to_csv(target / 'paired_metric_denominators.csv', index=False)
        for stratum, subset in [('all_measurements', detail),
                                ('all_four_exact', detail[detail.observed_D_exact])]:
            if subset.empty:
                continue
            points, draws, intervals = summary(subset, resamples=options['resamples'], seed=options['seed'])
            if stratum == 'all_measurements':
                # Same estimator on the saved compact sufficient statistics
                # and the common-index replay must agree to numerical precision.
                mappings = [('global_chemistry', curve[curve.k.eq(0) & curve.order.eq(orders[0]) & curve.model.eq('global_chemistry')]),
                            ('full_panel_profile', curve[curve.k.eq(n_columns) & curve.order.eq(orders[0]) & curve.model.eq('full_profile')]),
                            ('full_panel_level_chemistry', curve[curve.k.eq(n_columns) & curve.order.eq(orders[0]) & curve.model.eq('global_chemistry')])]
                mappings.extend((model, comparator_summary[comparator_summary.model.eq(model)]) for model in
                                ['chemistry_plus_target_score', 'chemistry_plus_selected_score'])
                for model, expected in mappings:
                    for _, cached in expected.iterrows():
                        for metric in PRIMARY_METRICS:
                            got = points[points.model.eq(model) & points.domain.eq(cached.domain) &
                                         points.weighting.eq(cached.weighting) & points.metric.eq(metric)]
                            if len(got) != 1 or not np.isclose(got.estimate.iloc[0], cached[metric], atol=1e-12, rtol=0, equal_nan=True):
                                raise ValueError(f'Saved-compact versus replay aggregation mismatch: {model}/{metric}/{cached.domain}/{cached.weighting}')
            for frame, filename in [(points, 'paired_estimates.csv'), (draws, 'paired_draws.csv.gz'),
                                    (intervals, 'paired_uncertainty.csv')]:
                frame['measurement_stratum'] = stratum
                frame.to_csv(target / (stratum + '_' + filename), index=False,
                             compression='gzip' if filename.endswith('.gz') else None)
        for frame, filename in [(pd.DataFrame(costs), 'growth_acquisition_costs.csv'),
                                (pd.DataFrame(assay_costs), 'query_assay_costs.csv'),
                                (pd.DataFrame(pair_costs), 'query_assay_pair_costs.csv'),
                                (pd.DataFrame(union_costs), 'query_profile_union_costs.csv'),
                                (pd.DataFrame(endpoint_checks), 'full_budget_invariance_checks.csv'),
                                (pd.concat(all_support, ignore_index=True), 'paired_assay_pair_support.csv')]:
            frame.to_csv(target / filename, index=False)
        receipts.append(dict(dataset=spec['dataset'], n_columns=n_columns, expected_folds=folds,
                             expected_orders=orders, **checks, k0_bitwise_identical=True,
                             paired_replayed_rectangles=len(detail), query_comparator_not_shared_k1=True))
    pd.concat(all_curves, ignore_index=True).to_csv(output / 'all_growth_preferences.csv', index=False)
    pd.concat(all_assays, ignore_index=True).to_csv(output / 'all_growth_assays.csv', index=False)
    after_verification = verify_manifest(pinned_manifest)
    if before_verification != after_verification:
        raise ValueError('Frozen source manifest changed during analysis')
    print(json.dumps(dict(source_verification='after', **after_verification)), flush=True)
    write_json(output / 'growth_analysis_receipt.json', dict(
        scope='Frozen-output aggregation and common-index prediction replay only; no new predictive fitting.',
        frozen_sources_before=before_verification, frozen_sources_after=after_verification,
        elapsed_seconds=time.monotonic() - started, datasets=receipts,
        config=config, config_path=str(config_path) if config_path else None,
        primary_domain='all_source_operational', primary_weighting='actual_assay_pair_equal',
        curve_variability='Six compositions, three seeds per policy; not statistical confidence intervals.',
        exact_all_k_scope='D_mse/D_rmse use exact rectangles in frozen compact summaries; exact direction/choice are regenerated only for zero-score/query/full-panel endpoint comparisons.',
        endpoint_order='First predeclared random order; full set of columns, no outcome-based choice.',
        source_files=[dict(path=str(path), sha256=digest) for path, digest in sorted(source.paths.items())]))
    print(json.dumps({'growth_analysis_complete': True, 'datasets': receipts, 'output_dir': str(output)}), flush=True)


def main(config_path=None):
    if config_path is None:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--config', type=Path, required=True)
        config_path = parser.parse_args().config
    config_path = Path(config_path).resolve()
    run(json.loads(config_path.read_text()), config_path)


if __name__ == '__main__':
    main()
