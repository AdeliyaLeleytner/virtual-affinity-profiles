"""Check the bundled compact output subset, including sanitized metadata hashes."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'provenance/completed_output_assets.json').read_text())
for row in manifest['files']:
    path=ROOT/row['path']
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
        raise SystemExit('Compact output mismatch: '+row['path'])
print('Verified',len(manifest['files']),'compact completed outputs; no source fits or replay.')
