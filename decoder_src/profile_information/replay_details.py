"""Regenerate every evaluated preference change from saved held-out predictions."""
from pathlib import Path
import argparse,json
import numpy as np
from . import data,evaluation,run


def main():
    p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--split-dir',type=Path,required=True);p.add_argument('--dataset',required=True);p.add_argument('--channel');p.add_argument('--scope',default='native_genes');p.add_argument('--policy',default='all_records');p.add_argument('--output-dir',type=Path,required=True);args=p.parse_args()
    d=data.load(args.input_dir,args.dataset,args.channel,args.scope,args.policy)
    split=json.loads((args.split_dir/'split.json').read_text());path=args.split_dir/'predictions.npz'
    if path.exists():
        saved=np.load(path,allow_pickle=False);pred={k:saved[k] for k in saved.files if k not in ['test_rows','test_assays']}
    else:
        saved=np.load(args.split_dir/'heldout_predictions.npz',allow_pickle=False);rr=saved['test_rows'];aa=saved['test_assays'];pred={}
        for name in saved.files:
            if name in ['test_rows','test_assays']:continue
            a=np.full(d['lower'].shape,np.nan);a[np.ix_(rr,aa)]=saved[name];pred[name]=a
    rr,aa=saved['test_rows'],saved['test_assays'];meta=dict(dataset=args.dataset,channel=args.channel or args.dataset,scope=args.scope,policy=args.policy,regime=split['regime'],split=split['id'],variant_set='detail_replay')
    frames=evaluation.evaluate(d,pred,rr,aa,meta,compact=False)
    args.output_dir.mkdir(parents=True,exist_ok=True);run.save_evaluation(args.output_dir,frames)
    run.write_json(args.output_dir/'replay_receipt.json',dict(source_predictions=str(path if path.exists() else args.split_dir/'heldout_predictions.npz'),n_preference_changes=len(frames[2]),source_outcomes='unchanged prepared records',new_model_fit=False))


if __name__=='__main__':main()
