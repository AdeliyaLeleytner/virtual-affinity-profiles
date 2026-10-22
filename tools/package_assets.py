"""Create deterministic public asset archives or restore them without local paths."""
from pathlib import Path, PurePosixPath
import argparse
import gzip
import hashlib
import json
import shutil
import tarfile

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / 'provenance/archive_assets.json'


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def clean_header(info):
    info.uid = info.gid = 0
    info.uname = info.gname = ''
    info.mtime = 0
    info.pax_headers = {}
    return info


def create(names):
    source_manifest = json.loads((ROOT / 'provenance/frozen_prediction_assets.json').read_text())
    runs = {r['run']: r for r in source_manifest['included_runs']}
    selected = names or list(runs)
    output = ROOT / 'release_assets'
    output.mkdir(exist_ok=True)
    entries = {r['run']: r for r in json.loads(CATALOG.read_text()).get('archives', [])} if CATALOG.exists() else {}
    for name in selected:
        if name not in runs:
            raise ValueError('Run has not been staged: ' + name)
        folder = ROOT / source_manifest['destination'] / name
        archive = output / (name + '.tar.gz')
        temporary = archive.with_suffix('.partial')
        with temporary.open('wb') as file_stream:
            with gzip.GzipFile(filename='', mode='wb', fileobj=file_stream, mtime=0, compresslevel=6) as compressed:
                with tarfile.open(fileobj=compressed, mode='w|', format=tarfile.PAX_FORMAT) as tar:
                    tar.add(folder, arcname=name, recursive=True, filter=clean_header)
        temporary.replace(archive)
        files = [p for p in folder.rglob('*') if p.is_file()]
        entry = dict(run=name, filename=archive.name, sha256=sha(archive), bytes=archive.stat().st_size,
                     files=len(files), uncompressed_bytes=sum(p.stat().st_size for p in files),
                     full_audit_pack=runs[name]['include_audit_artifacts'],
                     restoration='artifacts/frozen_predictions/' + name)
        entries[name] = entry
        print(json.dumps(entry), flush=True)
    CATALOG.write_text(json.dumps(dict(schema='public_frozen_assets_v1',
        status='local archives; download URLs will be assigned after publication',
        archives=[entries[k] for k in sorted(entries)]), indent=2) + '\n')


def restore(archive):
    archive = archive.resolve()
    entries = json.loads(CATALOG.read_text())['archives']
    matches = [entry for entry in entries if entry['filename'] == archive.name]
    if len(matches) != 1 or sha(archive) != matches[0]['sha256']:
        raise ValueError('Archive does not match the repository checksum catalog')
    entry = matches[0]
    destination = ROOT / 'artifacts/frozen_predictions'
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            path = PurePosixPath(member.name)
            if path.is_absolute() or '..' in path.parts or '\\' in member.name or not path.parts or path.parts[0] != entry['run']:
                raise ValueError('Unsafe or unexpected archive path')
            if not (member.isfile() or member.isdir()):
                raise ValueError('Archive must contain only regular files and directories')
            target = destination.joinpath(*path.parts)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            stream = tar.extractfile(member)
            if stream is None:
                raise ValueError('Missing member stream')
            temporary = target.with_name(target.name + '.restore-partial')
            with temporary.open('wb') as output:
                shutil.copyfileobj(stream, output)
            if target.exists() and sha(target) != sha(temporary):
                temporary.unlink()
                raise ValueError('Refusing to overwrite a different local artifact: ' + str(path))
            temporary.replace(target)
    print('Restored', entry['run'], 'without using original machine paths. Run python tools/verify_frozen_assets.py.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    build = sub.add_parser('create')
    build.add_argument('--run', action='append')
    extract = sub.add_parser('restore')
    extract.add_argument('archive', type=Path)
    args = parser.parse_args()
    if args.command == 'create':
        create(args.run)
    else:
        restore(args.archive)


if __name__ == '__main__':
    main()
