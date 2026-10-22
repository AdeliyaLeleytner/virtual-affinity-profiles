"""Hash the staged portable inputs/results before committing a replay config."""
from pathlib import Path
import argparse,hashlib,json


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    root=args.root.resolve();config=json.loads(args.config.read_text());spec=config.get('postprocessing',config)
    def resolve(name):
        path=Path(name);return path if path.is_absolute() else root/path
    dirs=[resolve(spec['input_dir'])]+[resolve(name) for name in spec['results']]
    output=args.output.resolve()
    files={}
    for folder in dirs:
        if not folder.is_dir():raise FileNotFoundError(folder)
        for path in sorted(folder.rglob('*')):
            if path.is_file() and path.resolve()!=output:
                try:name=str(path.resolve().relative_to(root))
                except ValueError:name=str(path.resolve())
                files[name]=hashlib.sha256(path.read_bytes()).hexdigest()
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(dict(schema='frozen_postprocessing_sources_v1',files=files),indent=2)+'\n')
    print(json.dumps(dict(event='snapshot_manifest_written',files=len(files),path=str(output))))


if __name__=='__main__':main()
