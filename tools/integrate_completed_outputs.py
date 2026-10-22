"""Stage compact completed analysis outputs; preserve raw source/archive trees."""
from pathlib import Path
import argparse
import hashlib
import json
import re

ROOT=Path(__file__).resolve().parents[1]
PROFILE={'paired_estimates.csv','paired_uncertainty.csv','heterogeneity.csv','continuous_trends.csv','illustrative_choices.csv','illustrative_source_variant_decisions.csv','within_assay_summary.csv','metric_defined_support.csv','source_dataset_inventory.csv','paired_association.csv','paired_association_uncertainty.csv'}
STRATA={'paired_estimates.csv','paired_uncertainty.csv','stratum_support.csv','stratum_status.csv'}


def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--project-root',type=Path,required=True);args=parser.parse_args()
    project=args.project_root.resolve()
    maps=[('artifacts/profile_postprocessing_run02','outputs/profile_postprocessing_run02'),('artifacts/profile_measurement_strata_run02','outputs/profile_measurement_strata_run02'),('artifacts/panel_growth_analysis','outputs/panel_growth_analysis'),('run-postprocessing-fixed/postprocessing','postprocessing'),('run-postprocessing-fixed/profile_information','replay_src/profile_information'),('run-growth-analysis/postprocessing','growth_analysis_support'),('run-growth-analysis/growth_analysis','growth_analysis'),('run-growth-canonical/profile_information','growth_src/profile_information'),('run-converged/profile_inputs','data'),('run-postprocessing-fixed/profile_inputs','data'),('preparation/data','data'),('results','artifacts/frozen_predictions')]
    def sanitize(text):
        for source,target in maps:text=text.replace(str(project/source),target)
        text=text.replace(str(project),'source-evidence')
        for folder in ['Users','home']:
            text=re.sub('/'+folder+r'/[^"\n]+?/([0-9a-f-]{36})/repo/',r'orx-source/\1/repo/',text)
        forbidden=['/'+part+'/' for part in ['Users','home','var/folders']]
        if any(value in text for value in forbidden):raise ValueError('Unmapped private path remains in output metadata')
        return text
    selections=[('profile_postprocessing_run02',PROFILE),('profile_measurement_strata_run02',STRATA),('panel_growth_analysis',None)]
    records=[]
    for name,allowed in selections:
        source_root=project/'artifacts'/name
        if not source_root.is_dir():raise FileNotFoundError(source_root)
        for source in sorted(source_root.rglob('*')):
            if not source.is_file() or 'replay' in source.relative_to(source_root).parts:continue
            relative=source.relative_to(source_root)
            keep=source.name in allowed if allowed is not None else source.suffix=='.csv'
            keep=keep or source.suffix=='.json'
            if not keep:continue
            # Draw/rectangle censuses are regenerated from the five prediction archives.
            if source.name.endswith('.csv.gz') or 'continuous_rectangle_data' in source.name:continue
            target=ROOT/'results'/name/relative;target.parent.mkdir(parents=True,exist_ok=True)
            original=source.read_bytes();text=original.decode('utf-8');public=sanitize(text).encode('utf-8')
            target.write_bytes(public)
            records.append(dict(source_id=str(source.relative_to(project)),path=str(target.relative_to(ROOT)),source_sha256=hashlib.sha256(original).hexdigest(),sha256=sha(target),bytes=len(public),byte_identical=public==original,metadata_path_derivative=public!=original))
    support=project/'manuscript/tables'
    for source in sorted(support.glob('decision_support*')):
        if source.suffix not in ['.csv','.json']:continue
        target=ROOT/'results/decision_support'/source.name;target.parent.mkdir(parents=True,exist_ok=True)
        original=source.read_bytes();public=sanitize(original.decode()).encode();target.write_bytes(public)
        records.append(dict(source_id=str(source.relative_to(project)),path=str(target.relative_to(ROOT)),source_sha256=hashlib.sha256(original).hexdigest(),sha256=sha(target),bytes=len(public),byte_identical=public==original,metadata_path_derivative=public!=original))
    note=support/'support_methods_note.tex'
    if note.exists():
        target=ROOT/'results/decision_support'/note.name;target.write_bytes(note.read_bytes())
        records.append(dict(source_id=str(note.relative_to(project)),path=str(target.relative_to(ROOT)),source_sha256=sha(note),sha256=sha(target),bytes=target.stat().st_size,byte_identical=True,metadata_path_derivative=False))
    manifest=dict(schema='compact_completed_outputs_v1',scope='Compact numeric summaries and sanitized provenance from completed authoritative analyses; raw rectangles/draws remain regenerable from frozen prediction archives.',postprocessing_run_id='c133ae0c-9f6a-4da6-a321-7ed9f79e585a',postprocessing_commit='5807705b7a17ea4fb438880f887cbfde967cee3a',source_datasets=24,measurement_source_tables=72,metadata_policy='Only original machine path strings are mapped to public aliases; every source and public hash is recorded. Whole-run receipts may reference additional omitted replay outputs.',files=records)
    (ROOT/'provenance/completed_output_assets.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(files=len(records),bytes=sum(r['bytes'] for r in records),metadata_derivatives=sum(r['metadata_path_derivative'] for r in records))),flush=True)


if __name__=='__main__':main()
