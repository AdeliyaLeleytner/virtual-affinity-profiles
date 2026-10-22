"""Pin and verify frozen inputs before and after a replay. Never reads outcomes numerically."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]


def digest(path):
    hasher=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):hasher.update(chunk)
    return hasher.hexdigest()


def verify_manifest(path,root=ROOT):
    manifest=json.loads(Path(path).read_text())
    for row in manifest['files']:
        file=Path(row['path']) if row['absolute'] else Path(root)/row['path']
        if not file.is_file() or file.stat().st_size!=row['bytes'] or digest(file)!=row['sha256']:
            raise ValueError(f'Frozen source verification failed: {file}')
    return dict(verified_files=len(manifest['files']),manifest_sha256=digest(path))


def build(config_path,output,root=ROOT):
    config=json.loads(Path(config_path).read_text())['growth_analysis']
    resolve=lambda p:(Path(root)/p).resolve()
    files={Path(config_path).resolve()}
    files.update(p for p in resolve(config['input_dir']).rglob('*') if p.is_file())
    for directory in ['growth_analysis','postprocessing','profile_information','study']:
        files.update(p for p in (Path(root)/directory).glob('*.py') if p.is_file())
    files.add(Path(root)/'environment.lock.json')
    for spec in config['datasets']:
        growth=resolve(spec['growth_dir'])
        for name in ['summary_preferences.csv','summary_assays.csv','query_score_comparator_summary.csv']:
            files.add(growth/name)
        files.add(growth.parent/'input_description.json')
        for fold in range(config['expected_folds']):
            folder=growth/f'mol_fold_{fold}'
            for seed in config['seeds']:
                for policy in ['random','family_stratified']:
                    order=folder/f'{policy}_{seed}'
                    files.add(order/'order.json')
                    for cell in order.glob('k*'):
                        if not cell.is_dir():continue
                        for name in ['heldout_predictions.npz','split.json','preference_context_metrics.csv.gz','assay_metrics.csv.gz','complete.json']:
                            files.add(cell/name)
            query=folder/'query_score_comparators'
            for name in ['heldout_predictions.npz','split.json','preference_context_metrics.csv.gz','assay_metrics.csv.gz','complete.json','fit_diagnostics.json','score_acquisition.json']:
                files.add(query/name)
    rows=[]
    for path in sorted(files):
        try:label=str(path.relative_to(root));absolute=False
        except ValueError:label=str(path);absolute=True
        rows.append(dict(path=label,absolute=absolute,bytes=path.stat().st_size,sha256=digest(path)))
    manifest=dict(scope='Byte-level frozen source pin for registered analysis; all saved growth predictions, splits, compact statistics and required acquisition metadata; no coefficient files needed.',files=rows)
    Path(output).write_text(json.dumps(manifest,indent=2)+'\n')
    return verify_manifest(output,root)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path,default=ROOT/'config.json');parser.add_argument('--write',type=Path);parser.add_argument('--verify',type=Path)
    args=parser.parse_args()
    if args.write:result=build(args.config,args.write)
    elif args.verify:result=verify_manifest(args.verify)
    else:parser.error('Choose --write or --verify')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
