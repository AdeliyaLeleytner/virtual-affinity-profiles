"""Stage the verified minimal cached figure bundle with explicit path derivatives."""
from pathlib import Path
import argparse,hashlib,json,re
ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);args=p.parse_args();project=args.project_root.resolve();source=project/'manuscript/figures';dest=ROOT/'figures';dest.mkdir(exist_ok=True)
    original=json.loads((source/'figure_bundle_manifest.json').read_text())
    names=set(original['files'])
    extras=['ASSET_PROVENANCE.md','FIGURE_PORTABILITY.md','figure_environment.json','requirements-figures-validated.txt','requirements-figure-assets-optional.txt','regenerate_illustration_assets.py','graphical_abstract_render_receipt.json']
    names.update(name for name in extras if (source/name).is_file())
    stems=['figure_geometry_readouts_support','figure_geometry_shared_shape','figure3_broad_paired_value','figure4_spd_paired_value','figure5_panel_acquisition','figure6_measured_choice_examples','figure_supplement_residual_correlation','figure_supplement_growth_level_chemistry','graphical_abstract_final']
    for stem in stems:
        for suffix in ['.pdf','.svg','.png','_caption.md','_render_receipt.json']:
            if (source/(stem+suffix)).is_file():names.add(stem+suffix)
    records=[]
    for name in sorted(names):
        path=source/name
        if name in original['files'] and sha(path)!=original['files'][name]:raise AssertionError('Authoring bundle changed: '+name)
        target=dest/name;target.parent.mkdir(parents=True,exist_ok=True);before=path.read_bytes();after=before
        if path.suffix in ['.json','.csv','.md','.log','.txt']:
            text=before.decode('utf-8').replace(str(project/'manuscript/figures'),'figures').replace(str(project/'artifacts/profile_postprocessing_run02'),'results/profile_postprocessing_run02').replace(str(project/'artifacts/panel_growth_analysis'),'results/panel_growth_analysis').replace(str(project),'source-evidence')
            for folder in ['Users','home']:text=re.sub('/'+folder+r'/[^"\n]+?/([0-9a-f-]{36})/repo/',r'orx-source/\1/repo/',text)
            if name=='figure_environment.json':
                environment=json.loads(text);environment['validated_interpreter']='python3.12';text=json.dumps(environment,indent=2)+'\n'
            if any('/'+part+'/' in text for part in ['Users','home','var/folders']):raise ValueError('Private path remains: '+name)
            after=text.encode('utf-8')
        target.write_bytes(after)
        records.append(dict(source_id='manuscript/figures/'+name,path='figures/'+name,source_sha256=hashlib.sha256(before).hexdigest(),sha256=sha(target),bytes=len(after),byte_identical=before==after,path_metadata_derivative=before!=after))
    public=dict(scope=original['scope'],validated_environment=original['validated_environment'],files={name:sha(dest/name) for name in original['files']},original_manifest_sha256=sha(source/'figure_bundle_manifest.json'),metadata_derivative_provenance='../provenance/figure_bundle_integration.json')
    (dest/'figure_bundle_manifest.json').write_text(json.dumps(public,indent=2)+'\n')
    (ROOT/'provenance/figure_bundle_integration.json').write_text(json.dumps(dict(scope='All six main figures, two supplements and final graphical abstract, cached renderer dependencies and attribution; no manuscript or scientific fits.',source_manifest_sha256=sha(source/'figure_bundle_manifest.json'),public_minimal_files=len(public['files']),files=records),indent=2)+'\n')
    (ROOT/'requirements-figures.txt').write_text((source/'requirements-figures-validated.txt').read_text())
    print(json.dumps(dict(files=len(records),minimal_dependencies=len(public['files']),bytes=sum(r['bytes'] for r in records))))


if __name__=='__main__':main()
