"""Package completed runs without fitting or editing their source artifacts.

Run from any directory. Source is the directory containing core_run02 etc.
Predictions, split JSON and compact metrics remain byte identical. Decoder
inheritance path metadata is an explicit derivative with both manifests retained. Omitted fitted coefficients and record-level prediction tables are
still listed in the original completion manifests as external audit artifacts.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
RUNS = ['core_run02', 'fingerprint_run01', 'both_cold_run01', 'growth_run02', 'decoder_sensitivity_run01']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')


def assert_portable(path):
    if path.suffix == '.json' or path.suffix == '.csv':
        text = path.read_text()
    elif path.name.endswith('.csv.gz'):
        with gzip.open(path, 'rt') as stream:
            text = stream.read()
    else:
        return
    if any('/' + folder + '/' in text for folder in ['Users', 'home']):
        raise ValueError('Review personal paths before packaging: ' + path.name)


def copy_exact(source, target, records):
    assert_portable(source)
    digest = sha(source)
    if target.exists() and sha(target) != digest:
        raise ValueError('Refusing to replace a different staged artifact: ' + str(target))
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        shutil.copyfile(source, target)
    if sha(target) != digest:
        raise ValueError('Copy hash mismatch: ' + str(target))
    records.append(dict(path=str(target.relative_to(ROOT)), sha256=digest, bytes=target.stat().st_size))


def stage_run(source, destination, name, include_audit_artifacts=False):
    run = source / name
    receipt = run / 'run_receipt.json'
    if not receipt.is_file() or json.loads(receipt.read_text()).get('status') != 'complete':
        raise ValueError('Run is not terminal complete: ' + name)
    signature = json.loads((run / 'analysis_signature.json').read_text())
    target = destination / name
    config = json.loads((run / 'effective_config.json').read_text())
    config['input_dir'] = 'data'
    config['output_dir'] = 'outputs/' + name
    if 'sensitivity' in config and 'frozen_results' in config['sensitivity']:
        config['sensitivity']['frozen_results'] = 'artifacts/frozen_predictions/core_run02'
    # Paths change, so the sanitized metadata is not labelled an original signature body.
    write_json(target / 'effective_config.json', config)
    write_json(target / 'frozen_source_provenance.json', dict(
        source_alias='results/' + name,
        original_analysis_signature=signature['signature'],
        original_signature_file_sha256=sha(run / 'analysis_signature.json'),
        original_effective_config_sha256=sha(run / 'effective_config.json'),
        input_sha256=signature['input_sha256'],
        module_sha256=signature['module_sha256'],
        sanitization='Original absolute paths omitted. Predictions/splits/compact metrics are unchanged. Decoder baseline_inheritance metadata is an explicit portable derivative; original complete.json retained as source_complete.original.json and both metadata hashes recorded.',
        original_gpu_inference_replayed=False,
    ))
    # The exact decoder source consumes these three historical identity fields.
    # This reduced portable record is explicitly not the original signed JSON body.
    write_json(target / 'analysis_signature.json', dict(
        schema='historical_signature_reference_v1',
        signature=signature['signature'],
        input_sha256=signature['input_sha256'],
        module_sha256=signature['module_sha256'],
        original_signature_file_sha256=sha(run / 'analysis_signature.json'),
        note='Historical signature identifier and byte manifests only; not a recomputable hash of this reduced portable file. See frozen_source_provenance.json.',
    ))
    records = []
    metadata_derivatives = []
    for marker in sorted(run.rglob('complete.json')):
        manifest = json.loads(marker.read_text())
        if manifest['analysis_signature'] != signature['signature']:
            raise ValueError('Split signature mismatch: ' + str(marker.relative_to(run)))
        for filename, expected in manifest['artifact_sha256'].items():
            artifact = marker.parent / filename
            if not artifact.is_file() or sha(artifact) != expected:
                raise ValueError('Source completion manifest mismatch: ' + str(artifact.relative_to(run)))
        folder_derivatives = []
        for artifact in sorted(marker.parent.iterdir()):
            if not artifact.is_file() or artifact.name == 'complete.json':
                continue
            if not include_audit_artifacts and (artifact.name.endswith('_model.npz') or artifact.name == 'record_predictions.csv.gz'):
                continue
            destination_file = target / artifact.relative_to(run)
            if artifact.name == 'baseline_inheritance.json':
                value = json.loads(artifact.read_text())
                old_path = str(value['frozen_directory'])
                relative = old_path.split('/results/core_run02/', 1)[-1]
                if relative == old_path:
                    raise ValueError('Unrecognized frozen baseline path layout')
                value['frozen_directory'] = 'artifacts/frozen_predictions/core_run02/' + relative
                write_json(destination_file, value)
                assert_portable(destination_file)
                derivative = dict(path=str(destination_file.relative_to(ROOT)), source_alias='results/' + name + '/' + str(artifact.relative_to(run)), original_sha256=sha(artifact), portable_sha256=sha(destination_file), transform='Only frozen_directory replaced by a repository-relative asset path; no predictions, arrays, rows, tuning or scalar results changed.')
                folder_derivatives.append(derivative)
                metadata_derivatives.append(derivative)
                records.append(dict(path=str(destination_file.relative_to(ROOT)), sha256=sha(destination_file), bytes=destination_file.stat().st_size))
            else:
                copy_exact(artifact, destination_file, records)
        public_marker = target / marker.relative_to(run)
        if folder_derivatives:
            copy_exact(marker, public_marker.with_name('source_complete.original.json'), records)
            public_manifest = dict(manifest)
            public_manifest['artifact_sha256'] = dict(manifest['artifact_sha256'])
            for item in folder_derivatives:
                public_manifest['artifact_sha256'][Path(item['path']).name] = item['portable_sha256']
            public_manifest['release_metadata_derivative'] = dict(original_complete_sha256=sha(marker), original_complete_file='source_complete.original.json', changes=folder_derivatives)
            write_json(public_marker, public_manifest)
            records.append(dict(path=str(public_marker.relative_to(ROOT)), sha256=sha(public_marker), bytes=public_marker.stat().st_size))
        else:
            copy_exact(marker, public_marker, records)
    # Include dataset-level aliases, descriptions and available aggregate summaries.
    for artifact in sorted(run.rglob('*')):
        if not artifact.is_file() or artifact.name in ['analysis_signature.json', 'effective_config.json']:
            continue
        if (artifact.parent / 'complete.json').is_file():
            continue
        if artifact.name in ['run_receipt.json', 'identical_score_aliases.json', 'identical_score_alias.json', 'input_description.json', 'summary_assays.csv', 'summary_preferences.csv', 'query_score_comparator_summary.csv', 'order.json']:
            copy_exact(artifact, target / artifact.relative_to(run), records)
    return dict(run=name, source_terminal_complete=True, include_audit_artifacts=include_audit_artifacts, metadata_derivatives=metadata_derivatives, original_analysis_signature=signature['signature'], files=len(records), bytes=sum(x['bytes'] for x in records), artifacts=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True, help='Source results directory; read-only')
    parser.add_argument('--run', choices=RUNS, action='append', required=True)
    parser.add_argument('--destination', type=Path, default=ROOT / 'artifacts/frozen_predictions')
    parser.add_argument('--include-audit-artifacts', action='store_true', help='Include every manifest-listed fit artifact; required for exact frozen decoder refits')
    args = parser.parse_args()
    source = args.source.resolve()
    destination = args.destination.resolve()
    if ROOT not in destination.parents:
        parser.error('Destination must be inside the release repository')
    if source == destination or source in destination.parents:
        parser.error('Destination must be outside the immutable source result tree')
    summary = [stage_run(source, destination, name, args.include_audit_artifacts) for name in args.run]
    manifest_path = ROOT / 'provenance/frozen_prediction_assets.json'
    prior = json.loads(manifest_path.read_text()).get('included_runs', []) if manifest_path.exists() else []
    combined = {item['run']: item for item in prior}
    combined.update({item['run']: item for item in summary})
    write_json(manifest_path, dict(
        status='selected completed runs copied and hash-verified; final postprocessor/figure integration is separate',
        destination=str(destination.relative_to(ROOT)),
        included_runs=[combined[name] for name in sorted(combined)],
        omitted='Without --include-audit-artifacts, fitted coefficients and per-record tables are omitted; exact decoder refits require the full core audit pack. Original inference weights and poses are never included.',
    ))
    print(json.dumps([dict(run=x['run'], files=x['files'], bytes=x['bytes']) for x in summary], indent=2))


if __name__ == '__main__':
    main()
