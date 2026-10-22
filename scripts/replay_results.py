"""Run published prediction-only analysis stages with their own frozen packages."""
from pathlib import Path
import argparse
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['postprocessing', 'measurement-strata', 'growth', 'all'])
    parser.add_argument('--config', type=Path, help='Committed alternate configuration for a single replay stage')
    args = parser.parse_args()
    if args.stage == 'all' and args.config:
        parser.error('--config is for one stage at a time')
    if args.stage == 'all':
        # Fresh interpreters keep the distinct stage-specific packages/imports separate.
        for stage in ['postprocessing', 'measurement-strata', 'growth']:
            subprocess.run([sys.executable, str(Path(__file__).resolve()), stage], check=True)
        return
    source = 'growth_src' if args.stage == 'growth' else 'replay_src'
    sys.path.insert(0, str(ROOT / source))
    sys.path.insert(0, str(ROOT))
    if args.stage == 'postprocessing':
        from postprocessing.dispatch import main as run
        run(args.config.resolve() if args.config else ROOT / 'configs/postprocessing.json')
    elif args.stage == 'measurement-strata':
        from postprocessing.measurement_strata import run_from_config
        run_from_config(args.config.resolve() if args.config else ROOT / 'configs/measurement_strata.json', project_root=ROOT)
    else:
        from growth_analysis.analyze import main as run
        run(args.config.resolve() if args.config else ROOT / 'configs/growth_analysis.json')


if __name__ == '__main__':
    main()
