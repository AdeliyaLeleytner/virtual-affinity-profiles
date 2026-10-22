"""Verify public frozen artifacts, exact fit sources and explicit metadata derivatives."""
from pathlib import Path
import hashlib
import argparse
import json

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = {'core_run02': 'src', 'fingerprint_run01': 'src', 'both_cold_run01': 'src',
            'growth_run02': 'growth_src', 'decoder_sensitivity_run01': 'decoder_src'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', choices=sorted(PACKAGES), action='append', help='Verify only one selected restored pack; repeat if needed')
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'provenance/frozen_prediction_assets.json').read_text())
    summaries = []
    for run in manifest['included_runs']:
        if args.run and run['run'] not in args.run:
            continue
        folder = ROOT / manifest['destination'] / run['run']
        source = json.loads((folder / 'frozen_source_provenance.json').read_text())
        for name, expected in source['input_sha256'].items():
            if sha(ROOT / 'data' / name) != expected:
                raise AssertionError('Prepared input mismatch: ' + name)
        for name, expected in source['module_sha256'].items():
            if sha(ROOT / PACKAGES[run['run']] / 'profile_information' / name) != expected:
                raise AssertionError('Frozen fit-module mismatch: ' + run['run'] + ':' + name)
        for artifact in run['artifacts']:
            if sha(ROOT / artifact['path']) != artifact['sha256']:
                raise AssertionError('Public artifact mismatch: ' + artifact['path'])
        for change in run.get('metadata_derivatives', []):
            portable = ROOT / change['path']
            public_complete = json.loads((portable.parent / 'complete.json').read_text())
            original_path = portable.parent / public_complete['release_metadata_derivative']['original_complete_file']
            if sha(original_path) != public_complete['release_metadata_derivative']['original_complete_sha256']:
                raise AssertionError('Original complete manifest changed')
            original_complete = json.loads(original_path.read_text())
            if original_complete['artifact_sha256'][portable.name] != change['original_sha256']:
                raise AssertionError('Original metadata hash differs from original completion')
            if sha(portable) != change['portable_sha256']:
                raise AssertionError('Portable metadata hash changed')
            if public_complete['artifact_sha256'][portable.name] != change['portable_sha256']:
                raise AssertionError('Public derivative completion mismatch')
            for name in ['predictions.npz', 'split.json', 'preference_context_metrics.csv.gz']:
                if sha(portable.parent / name) != original_complete['artifact_sha256'][name]:
                    raise AssertionError('Scientific artifact changed while removing metadata path: ' + name)
        summaries.append(dict(run=run['run'], verified_artifacts=len(run['artifacts']),
                              explicit_metadata_derivatives=len(run.get('metadata_derivatives', [])),
                              full_audit_pack=run['include_audit_artifacts']))
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    main()
