"""Deterministic one-block saved-decoder replay check; no scientific fitting."""
from pathlib import Path
import argparse,json,sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
sys.path.insert(0,str(HERE))
import postprocess


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,default=HERE.parent/'config.json');p.add_argument('--output-dir',type=Path,default=HERE.parent/'artifacts/nested_name_validation/smoke');args=p.parse_args()
    root=args.config.resolve().parent;config=json.loads(args.config.read_text())['postprocessing']
    resolve=lambda value:Path(value) if Path(value).is_absolute() else root/value
    sources=[resolve(value) for value in config['results']]
    source=next(value for value in sources if value.name=='decoder_sensitivity_run01')
    # Selection uses lexical block identity, not an outcome or model metric.
    marker=sorted(source.glob('*__all_records/core/*/*/complete.json'))[0]
    folder=marker.parent;dataset_key=folder.parents[2].name
    dataset,channel,scope,policy=dataset_key.split('__')
    d=postprocess.data.load(resolve(config['input_dir']),dataset,channel if dataset=='spd' else None,scope,policy)
    saved=np.load(folder/'predictions.npz',allow_pickle=False)
    models=sorted(k for k in saved.files if k not in ['test_rows','test_assays'])
    nested=sorted(name for name in models if '__' in name)
    if not nested:raise AssertionError('Chosen deterministic block has no nested model names')
    detail=postprocess.replay_one(d,folder,args.output_dir/'replay')
    contexts=postprocess.evaluation.detail_contexts(detail)
    summary=postprocess.evaluation.summarize_preferences(detail)
    if set(contexts.model)!=set(models) or set(summary.model)!=set(models):raise AssertionError('Full model names did not survive helper summaries')
    replay=json.loads((args.output_dir/'replay/replay_receipt.json').read_text())
    receipt=dict(status='passed',selection='first lexical completed decoder block',source=str(folder),source_prediction_sha256=postprocess.sha(folder/'predictions.npz'),source_compact_sha256=postprocess.sha(folder/'preference_context_metrics.csv.gz'),models=models,nested_models=nested,n_rectangle_rows=len(detail),n_context_rows=len(contexts),n_preference_summary_rows=len(summary),detail_compact_numeric_max_abs_difference=replay['compact_numeric_max_abs_difference'],evaluator_sha256=postprocess.sha(Path(postprocess.evaluation.__file__)),new_model_fit=False,whole_scientific_replay=False)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    postprocess.write_json(args.output_dir/'smoke_receipt.json',receipt)
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()
