"""Run existing software tests with isolated stage-specific import paths."""
from pathlib import Path
import argparse
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def run(paths, arguments):
    environment = os.environ.copy()
    environment['PYTHONPATH'] = os.pathsep.join(str(ROOT / p) for p in paths)
    subprocess.run([sys.executable, *arguments], cwd=ROOT, env=environment, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--with-assets', action='store_true')
    args = parser.parse_args()
    for source, tests in [('src', 'tests'), ('growth_src', 'growth_tests'), ('decoder_src', 'decoder_tests')]:
        run([source], ['-m', 'unittest', 'discover', '-s', tests, '-p', 'test_*.py', '-v'])
    run(['replay_src', 'postprocessing'], ['-m', 'unittest', 'discover', '-s', 'postprocessing', '-p', 'test_*.py', '-v'])
    run(['growth_src', 'growth_analysis_support', '.'], ['-m', 'unittest', 'growth_analysis.test_technical', 'growth_analysis.test_stage.StageTests.test_manifest_rejects_changed_source', '-v'])
    run([], ['scripts/verify_bundle.py'])
    if args.with_assets:
        run([], ['tools/verify_frozen_assets.py'])
        run([], ['scripts/replay_geometry.py'])


if __name__ == '__main__':
    main()
