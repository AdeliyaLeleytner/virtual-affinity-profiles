"""Prediction-only measurement strata; primary all-identified outputs stay intact."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys

import numpy as np
import pandas as pd
import scipy

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from fast_resampling import summary
STRATA = ('exact_four_cell', 'censor_involving_identified_D')
METRICS = (
    'D_direction_concordance', 'D_predicted_tie_fraction',
    'pareto_correct_fraction', 'pareto_reversed_fraction',
    'pareto_tradeoff_fraction', 'pareto_tie_fraction',
)
LABELS = ('source_run', 'dataset', 'channel', 'scope', 'policy', 'regime', 'comparison_support')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def strata(frame):
    """Measurement classes, with exact no-change retained as undefined direction."""
    exact = frame.observed_D_exact
    if exact.dtype != bool:
        raise ValueError('observed_D_exact must be boolean')
    direction = frame.observed_D_direction
    if not direction.isin([-1, 0, 1]).all():
        raise ValueError('D direction must be source-identified -1/0/1')
    return {
        STRATA[0]: frame.loc[exact].copy(),
        STRATA[1]: frame.loc[~exact & direction.ne(0)].copy(),
    }


def support_counts(frame):
    rows = []
    for name, sub in strata(frame).items():
        for domain, g in [('all_source_operational', sub),
                          ('metadata_matched_operational', sub[sub.quantitative_context_matched])]:
            def pairs(mask):
                return len(g.loc[mask, ['assay_id_a', 'assay_id_b']].drop_duplicates())
            identified = g.observed_D_direction.ne(0)
            choice = g.pareto_identified
            rows.append(dict(
                measurement_stratum=name, domain=domain,
                observed_rectangle_rows=len(g), observed_exact_rows=int(g.observed_D_exact.sum()),
                observed_identified_D_rows=int(identified.sum()),
                observed_exact_nochange_D_rows=int((g.observed_D_exact & ~identified).sum()),
                observed_joint_choice_rows=int(choice.sum()),
                observed_actual_assay_pairs=len(g[['assay_id_a', 'assay_id_b']].drop_duplicates()),
                observed_direction_defined_assay_pairs=pairs(identified),
                observed_choice_defined_assay_pairs=pairs(choice),
                observed_molecules=len(set(g.molecule_a) | set(g.molecule_b)),
                observed_chemical_groups=len(set(g.chemical_group_a.astype(str)) | set(g.chemical_group_b.astype(str))),
                observed_genes=len(set(g.gene_a) | set(g.gene_b)),
                undefined_reason='' if len(g) else 'no eligible rectangles in this measurement stratum/domain',
            ))
    return pd.DataFrame(rows)


def analyse(frame, resamples=400, seed=20261022):
    estimates, draws, uncertainty, statuses = [], [], [], []
    for name, sub in strata(frame).items():
        statuses.append(dict(measurement_stratum=name, eligible_rectangle_rows=len(sub),
                             status='defined_observation_pool' if len(sub) else 'empty_observation_pool'))
        if sub.empty:
            continue
        results = summary(sub, resamples=resamples, seed=seed)
        for target, result in zip([estimates, draws, uncertainty], results):
            if result.empty:
                continue
            selected = result.loc[result.metric.isin(METRICS)].copy()
            selected.insert(0, 'measurement_stratum', name)
            target.append(selected)
    def combine(frames):
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=['measurement_stratum', 'domain', 'weighting', 'model', 'metric'])
    return combine(estimates), combine(draws), combine(uncertainty), support_counts(frame), pd.DataFrame(statuses)


def source_labels(source, input_root):
    relative = source.relative_to(input_root)
    if len(relative.parts) < 5:
        raise ValueError('Expected source-run/dataset/regime/support/continuous_rectangle_data.csv.gz: ' + str(relative))
    run = '/'.join(relative.parts[:-4])
    dataset_key, regime, comparison_support, _ = relative.parts[-4:]
    parts = dataset_key.split('__')
    if len(parts) != 4:
        raise ValueError('Unrecognized dataset key: ' + dataset_key)
    return dict(zip(LABELS, [run] + parts + [regime, comparison_support]))


def process_source(source, input_root, output_root, resamples, seed):
    relative = source.relative_to(input_root)
    dest = output_root / relative.parent
    labels = source_labels(source, input_root)
    signature = dict(
        schema='measurement_strata_v1', source_sha256=sha(source),
        module_sha256={p.name: sha(p) for p in [Path(__file__), HERE/'postprocess.py', HERE/'fast_resampling.py']},
        resamples=resamples, seed=seed, strata=list(STRATA), metrics=list(METRICS),
        endpoint_pool='observed endpoints within each stratum/domain; same counts for all paired arms',
    )
    receipt_path = dest/'measurement_strata_receipt.json'
    if receipt_path.exists():
        old = json.loads(receipt_path.read_text())
        outputs = old.get('output_sha256', {})
        if old.get('signature') == signature and outputs and all((dest/name).is_file() and sha(dest/name) == value for name, value in outputs.items()):
            print(json.dumps(dict(event='measurement_strata_resume_skip', source=str(relative))), flush=True)
            return 'resume_skip'
    frame = pd.read_csv(source, float_precision='round_trip')
    outputs = analyse(frame, resamples, seed)
    names = ['paired_estimates.csv', 'paired_resampling_draws.csv.gz', 'paired_uncertainty.csv', 'stratum_support.csv', 'stratum_status.csv']
    dest.mkdir(parents=True, exist_ok=True)
    for name, result in zip(names, outputs):
        for key, value in reversed(list(labels.items())):
            result.insert(0, key, value)
        result.to_csv(dest/name, index=False)
    receipt = dict(
        status='complete', source=str(source), source_labels=labels, signature=signature,
        new_model_training=False, scientific_thresholds=None,
        primary_all_identified_outputs_modified=False,
        input_rows=len(frame), input_identified_D_rows=int(frame.observed_D_direction.ne(0).sum()),
        omitted_censor_involving_unidentified_D_rows=int((~frame.observed_D_exact & frame.observed_D_direction.eq(0)).sum()),
        output_sha256={name: sha(dest/name) for name in names},
        versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
    )
    receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(dict(event='measurement_strata_complete', source=str(relative), rows=len(frame))), flush=True)
    return 'complete'


def run_from_config(config_path, project_root=None):
    root = Path(project_root).resolve() if project_root else HERE.parent
    config_path = Path(config_path)
    config = json.loads(config_path.read_text())
    spec = config.get('measurement_strata', config)
    def resolve(value):
        path = Path(value)
        return path.resolve() if path.is_absolute() else (root/path).resolve()
    input_root, output_root = resolve(spec['input_dir']), resolve(spec['output_dir'])
    if output_root == input_root or input_root in output_root.parents:
        raise ValueError('Measurement output must be separate from the immutable primary replay tree')
    resamples, seed = int(spec.get('resamples', 400)), int(spec.get('seed', 20261022))
    if resamples < 0:
        raise ValueError('resamples must be nonnegative')
    sources = [input_root/name for name in spec['files']] if spec.get('files') else sorted(input_root.rglob('continuous_rectangle_data.csv.gz'))
    if not sources:
        raise ValueError('No final continuous rectangle replay tables found')
    states = [process_source(source, input_root, output_root, resamples, seed) for source in sources]
    output_root.mkdir(parents=True, exist_ok=True)
    receipt = dict(status='complete', config_sha256=sha(config_path), sources=len(sources), completed=states.count('complete'), resumed=states.count('resume_skip'), resamples=resamples, seed=seed, scientific_thresholds=None, new_model_training=False, source_sha256={str(p.relative_to(input_root)):sha(p) for p in sources})
    (output_root/'measurement_strata_run_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    return receipt


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--project-root', type=Path)
    args = p.parse_args()
    run_from_config(args.config, args.project_root)


if __name__ == '__main__':
    main()
