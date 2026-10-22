"""Pin restored portable assets after validating their original scientific hashes."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def files_in(folder):
    if not folder.is_dir():
        raise FileNotFoundError(folder)
    return {path.resolve() for path in folder.rglob('*') if path.is_file() and '__pycache__' not in path.parts}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['all', 'postprocessing', 'growth'], default='all')
    parser.add_argument('--config', type=Path, help='Alternate postprocessing config, for a bounded restored-source selection')
    args = parser.parse_args()
    if args.config and args.stage != 'postprocessing':
        parser.error('--config requires --stage postprocessing')
    post_config = args.config.resolve() if args.config else ROOT / 'configs/postprocessing.json'
    subprocess.run([sys.executable, str(ROOT / 'scripts/verify_bundle.py')], check=True)
    runs = ['growth_run02'] if args.stage == 'growth' else ['core_run02', 'fingerprint_run01', 'both_cold_run01', 'decoder_sensitivity_run01']
    if args.stage == 'postprocessing' and args.config:
        selected=json.loads(post_config.read_text())['postprocessing']['results']
        runs=[Path(name).name for name in selected]
    if args.stage == 'all':
        runs.append('growth_run02')
    command = [sys.executable, str(ROOT / 'tools/verify_frozen_assets.py')]
    for run in runs:
        command.extend(['--run', run])
    subprocess.run(command, check=True)
    if args.stage in ['all', 'postprocessing']:
        config = json.loads(post_config.read_text())['postprocessing']
        paths = files_in(ROOT / config['input_dir'])
        for name in config['results']:
            paths |= files_in(ROOT / name)
        body = dict(schema='frozen_postprocessing_sources_v1',
                    note='Portable metadata derivatives are documented in frozen_prediction_assets.json; original scientific bytes were checked before pinning.',
                    files={str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)})
        (ROOT / config['source_manifest']).write_text(json.dumps(body, indent=2) + '\n')
        print('Pinned postprocessing:', len(paths), 'portable files')
    if args.stage in ['all', 'growth']:
        config_path = ROOT / 'configs/growth_analysis.json'
        config = json.loads(config_path.read_text())['growth_analysis']
        paths = {config_path}
        for folder in [ROOT / config['input_dir'], ROOT / config['frozen_code_dir'],
                       ROOT / config['postprocessing_code_dir'], ROOT / 'growth_analysis']:
            paths |= files_in(folder)
        for spec in config['datasets']:
            paths |= files_in((ROOT / spec['growth_dir']).parent)
        body = dict(scope='Portable frozen growth replay source pin; exact scientific arrays and source modules verified before pinning. No source inference or predictive fitting.',
                    files=[dict(path=str(p.relative_to(ROOT)), absolute=False, bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(paths)])
        (ROOT / config['source_manifest']).write_text(json.dumps(body, indent=2) + '\n')
        print('Pinned growth analysis:', len(paths), 'portable files')


if __name__ == '__main__':
    main()
