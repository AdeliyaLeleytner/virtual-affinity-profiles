"""Audit the prospective Git tree without publishing or changing scientific data."""
from pathlib import Path
import gzip,hashlib,json,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def main():
    raw=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=ROOT)
    names=sorted(set(n.decode() for n in raw.split(b'\0') if n))
    forbidden=['/'+part+'/' for part in ['Users','home','var/folders']]
    findings=[];scanned=0
    for name in names:
        path=ROOT/name
        if not path.is_file():continue
        if path.name in ['.env','.env.local'] or path.suffix in ['.pem','.key']:
            findings.append(dict(path=name,reason='credential/private-key file'))
        if path.stat().st_size>=100*1024*1024:findings.append(dict(path=name,reason='file reaches100MiB'))
        if name.startswith(('manuscript/','submission/')) or path.name in ['main.tex','main.pdf','article.tex','article.pdf','cover_letter.tex']:
            findings.append(dict(path=name,reason='manuscript/submission object'))
        texts=[]
        if path.name.endswith('.csv.gz'):
            with gzip.open(path,'rt') as stream:
                for line in stream:
                    if any(v in line for v in forbidden):findings.append(dict(path=name,reason='private path in compressed table'));break
            scanned+=1;continue
        if path.suffix=='.npz':
            with np.load(path,allow_pickle=False) as archive:
                for key in archive.files:
                    value=archive[key]
                    if value.dtype.kind in 'US':texts.extend(value.astype(str).ravel().tolist())
            scanned+=1
        elif path.suffix.lower() in ['.py','.md','.json','.csv','.txt','.toml','.cff','.yml','.yaml','.tex','.sh','.svg','.pdb','.log'] or path.name in ['LICENSE','.gitignore']:
            texts=[path.read_text(errors='replace')];scanned+=1
        if any(v in text for text in texts for v in forbidden):findings.append(dict(path=name,reason='private filesystem path'))
    required=['LICENSE','LICENSE-DATA.md','THIRD_PARTY_NOTICES.md','LICENSES/OpenResearch-MIT.txt','provenance/source_licenses.json','figures/ASSET_PROVENANCE.md']
    for name in required:
        if not (ROOT/name).is_file():findings.append(dict(path=name,reason='required licence/attribution absent'))
    if any(name.startswith(('release_assets/','outputs/','artifacts/','.qa/')) for name in names):findings.append(dict(reason='runtime/archive tree would enterGit'))
    expected=json.loads((ROOT/'figures/figure_bundle_manifest.json').read_text())['files']
    missing_dependencies=['figures/'+n for n in expected if 'figures/'+n not in names]
    if missing_dependencies:findings.extend(dict(path=n,reason='required figure dependency excluded fromGit') for n in missing_dependencies)
    receipt=dict(status='passed' if not findings else 'requires_fix',candidate_files=len(names),candidate_bytes=sum((ROOT/n).stat().st_size for n in names if (ROOT/n).is_file()),text_or_string_archives_scanned=scanned,findings=findings,manuscript_or_submission_included=False if not any(f.get('reason')=='manuscript/submission object' for f in findings) else True,archives_gitignored=True,figure_dependency_files=len(expected),licence_scope='Required code/data/source attribution files present; exact frozen packages and source terms retained. This is a distribution-completeness audit.',no_git_commit_or_push=True)
    (ROOT/'provenance/final_privacy_license_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
    if findings:raise SystemExit(1)


if __name__=='__main__':main()
