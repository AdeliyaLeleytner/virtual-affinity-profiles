"""Fixed-command OpenResearch dispatch for committed postprocessing config."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import pandas as pd

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import postprocess
from fast_resampling import summary
from common_context_association import main as association_main
from metric_support import main as support_main


def resolve(root,value):
    path=Path(value)
    return path if path.is_absolute() else root/path


def main(config_path=None):
    root=HERE.parent
    config_path=Path(config_path) if config_path else root/'config.json'
    config=json.loads(config_path.read_text())
    spec=config.get('postprocessing',config)
    input_dir=resolve(root,spec['input_dir'])
    results=[resolve(root,p) for p in spec['results']]
    output=resolve(root,spec['output_dir'])
    manifest_path=resolve(root,spec['source_manifest'])
    manifest=json.loads(manifest_path.read_text())
    def verify_sources():
        expected_paths={resolve(root,name).resolve() for name in manifest['files']}
        actual_paths={p.resolve() for folder in [input_dir]+results for p in folder.rglob('*') if p.is_file() and p.resolve()!=manifest_path.resolve()}
        if actual_paths!=expected_paths:raise RuntimeError('Frozen source file membership changed')
        for name,expected in manifest['files'].items():
            path=resolve(root,name)
            actual=hashlib.sha256(path.read_bytes()).hexdigest()
            if actual!=expected:raise RuntimeError('Frozen source hash mismatch: '+name)
    verify_sources()
    for source in results:
        receipt=json.loads((source/'run_receipt.json').read_text())
        if receipt.get('status')!='complete':raise RuntimeError('Complete source run required: '+str(source))
        effective=json.loads((source/'effective_config.json').read_text())
        for dataset in effective['datasets']:
            name='__'.join([dataset['dataset'],dataset.get('channel',dataset['dataset']),dataset.get('scope','native_genes'),dataset.get('policy','all_records')])
            if spec.get('datasets') and name not in spec['datasets']:continue
            folder=source/name
            if not folder.is_dir():raise RuntimeError('Expected dataset snapshot missing: '+str(folder))
    args=argparse.Namespace(input_dir=input_dir,output_dir=output,dataset=spec.get('datasets'),allow_partial=False,resamples=int(spec.get('resamples',400)),seed=int(spec.get('seed',20261022)))
    postprocess.summary=summary
    inventory=[]
    for source in results:
        for folder in sorted(source.glob('*__all_records')):
            if args.dataset and folder.name not in args.dataset:continue
            if (folder/'identical_score_alias.json').exists():
                alias=json.loads((folder/'identical_score_alias.json').read_text())
                inventory.append(dict(source_run=source.name,dataset=folder.name,status='identical_score_alias',alias_of=alias['alias_of'],independent_prediction_evidence=False))
                continue
            if list(folder.glob('core/*/*/complete.json')):
                postprocess.process_dataset(source,folder,args)
                state=json.loads((output/source.name/folder.name/'dataset_state.json').read_text())
                if not state['dataset_complete']:raise RuntimeError('Incomplete expected split coverage: '+str(folder))
                dataset_output=output/source.name/folder.name
                assay_frames=[pd.read_csv(path) for path in sorted((dataset_output/'replay').glob('*/*/assay_metrics.csv.gz'))]
                assay_frames=[f for f in assay_frames if len(f)]
                if assay_frames:
                    assay_frame=pd.concat(assay_frames,ignore_index=True)
                    postprocess.save(assay_frame,dataset_output/'within_assay_metrics.csv.gz')
                    postprocess.save(postprocess.evaluation.summarize_assays(assay_frame),dataset_output/'within_assay_summary.csv')
                inventory.append(dict(source_run=source.name,dataset=folder.name,status='replayed_complete',alias_of='',independent_prediction_evidence=True,expected_splits=state.get('expected_splits'),completed_splits=state.get('completed_splits')))
            else:raise RuntimeError('No completed snapshot splits: '+str(folder))
    # Sidecars use the same saved rectangle census and technical draw settings.
    previous=list(sys.argv)
    try:
        sys.argv=['metric_support','--output-dir',str(output),'--save',str(output/'metric_defined_support.csv')]
        support_main()
        sys.argv=['common_context_association','--output-dir',str(output),'--save-dir',str(output/'common_association'),'--resamples',str(args.resamples),'--seed',str(args.seed)]
        association_main()
    finally:sys.argv=previous
    verify_sources()
    module_files=sorted(HERE.glob('*.py'))
    output.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(inventory).to_csv(output/'source_dataset_inventory.csv',index=False)
    postprocess.write_json(output/'dispatch_receipt.json',dict(status='complete',new_model_training=False,scientific_success=None,config_sha256=postprocess.sha(config_path),source_manifest_sha256=postprocess.sha(manifest_path),source_files_checked=len(manifest['files']),source_datasets=len(inventory),independent_replayed_datasets=sum(r['status']=='replayed_complete' for r in inventory),identical_score_aliases=sum(r['status']=='identical_score_alias' for r in inventory),module_sha256={p.name:postprocess.sha(p) for p in module_files},resamples=args.resamples,seed=args.seed,source_data_unchanged_after_replay=True))
    print(json.dumps(dict(event='postprocessing_complete',output_dir=str(output),scientific_success=None,new_model_training=False)),flush=True)


if __name__=='__main__':
    main()
