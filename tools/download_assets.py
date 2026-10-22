"""Download versioned frozen-source archives and verify pinned SHA-256 values."""
from pathlib import Path
import argparse,hashlib,json,urllib.request

ROOT=Path(__file__).resolve().parents[1]


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',action='append');p.add_argument('--verify-only',action='store_true');args=p.parse_args()
    entries=json.loads((ROOT/'provenance/archive_assets.json').read_text())['archives']
    if args.run:
        missing=set(args.run)-{e['run'] for e in entries}
        if missing:raise ValueError('Unknown archive names: '+str(sorted(missing)))
        entries=[e for e in entries if e['run'] in args.run]
    destination=ROOT/'release_assets';destination.mkdir(exist_ok=True)
    for entry in entries:
        path=destination/entry['filename']
        if path.is_file() and digest(path)==entry['sha256']:
            print('Verified existing',entry['filename']);continue
        if args.verify_only:raise ValueError('Archive missing or hash mismatch: '+entry['filename'])
        temporary=path.with_suffix(path.suffix+'.download-partial')
        request=urllib.request.Request(entry['download_url'],headers={'User-Agent':'virtual-affinity-profiles/1.0.0'})
        try:
            with urllib.request.urlopen(request) as incoming,temporary.open('wb') as output:
                for block in iter(lambda:incoming.read(1024*1024),b''):output.write(block)
            if digest(temporary)!=entry['sha256']:raise ValueError('Downloaded archive hash mismatch: '+entry['filename'])
            temporary.replace(path)
        finally:
            if temporary.exists():temporary.unlink()
        print('Downloaded and verified',entry['filename'])


if __name__=='__main__':main()
