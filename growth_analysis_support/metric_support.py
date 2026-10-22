"""Metric-specific defined denominators, without changing any estimate."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def audit_frame(frame):
    rows=[]
    models=sorted(c[:-3] for c in frame if c.endswith('__D') and frame[c].notna().any())
    for domain,g in [('all_source_operational',frame),('metadata_matched_operational',frame[frame.quantitative_context_matched])]:
        if g.empty:continue
        for model in models:
            pair=g[model+'__D'].notna() & g.global_chemistry__D.notna()
            f=g[pair]
            for metric,key in [('D_direction_concordance','D_sign_credit'),('D_rmse','D_mse'),('pareto_correct_fraction','pareto_correct'),('symmetrized_D_correlation','D')]:
                if model.startswith('native_') and metric=='D_rmse':continue
                if metric=='symmetrized_D_correlation':
                    ex=f[f.observed_D_exact].copy();ex['_y2']=ex.observed_D_lower**2
                    by=ex.groupby(['assay_id_a','assay_id_b'])
                    yy=by._y2.sum();a=(yy>0)&(by[model+'__D2'].sum()>0);b=(yy>0)&(by.global_chemistry__D2.sum()>0)
                    narm,nbase,ncommon=int(a.sum()),int(b.sum()),int((a&b).sum())
                    defined=int(len(ex))
                else:
                    a=f[model+'__'+key].notna();b=f['global_chemistry__'+key].notna()
                    def count(mask):return f.loc[mask,['assay_id_a','assay_id_b']].drop_duplicates().shape[0]
                    narm,nbase,ncommon=count(a),count(b),count(a&b);defined=int((a&b).sum())
                rows.append(dict(domain=domain,model=model,metric=metric,n_total_rectangle_rows=len(g),n_paired_prediction_rows=len(f),n_metric_defined_rectangle_rows=defined,n_defined_assay_pairs_arm=narm,n_defined_assay_pairs_baseline=nbase,n_defined_assay_pairs_common=ncommon,defined_context_sets_agree=narm==nbase==ncommon))
    return pd.DataFrame(rows)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--save',type=Path,required=True);args=p.parse_args()
    frames=[]
    for path in sorted(args.output_dir.rglob('continuous_rectangle_data.csv.gz')):
        f=audit_frame(pd.read_csv(path))
        f.insert(0,'source',str(path.relative_to(args.output_dir)));frames.append(f)
    args.save.parent.mkdir(parents=True,exist_ok=True)
    (pd.concat(frames,ignore_index=True) if frames else pd.DataFrame()).to_csv(args.save,index=False)


if __name__=='__main__':main()
