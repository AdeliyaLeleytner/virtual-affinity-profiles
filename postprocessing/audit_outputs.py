"""Independent algebraic checks over generated tables, without fitting."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd


def main():
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);args=p.parse_args()
    checks=[]
    for path in sorted(args.output_dir.rglob('paired_estimates.csv')):
        f=pd.read_csv(path)
        native=f.model.str.startswith('native_')
        assert not (native&f.metric.isin(['D_rmse','D_mse'])).any(),path
        base=f[f.model.eq('global_chemistry')].arm_minus_baseline.dropna()
        assert np.allclose(base,0,rtol=0,atol=1e-12),path
        valid=f.estimate.notna()&f.baseline_estimate.notna()
        assert np.allclose(f.loc[valid,'arm_minus_baseline'],f.loc[valid,'estimate']-f.loc[valid,'baseline_estimate'],rtol=1e-10,atol=1e-12),path
        checks.append(dict(file=str(path),check='no native quantitative loss; paired delta algebra; baseline zero',rows=len(f),undefined_model_rows=int(f.estimate.isna().sum())))
    for path in sorted(args.output_dir.rglob('paired_resampling_draws.csv.gz')):
        try:f=pd.read_csv(path)
        except pd.errors.EmptyDataError:continue
        if not len(f):continue
        valid=f.arm_estimate.notna()&f.baseline_estimate.notna()
        assert np.allclose(f.loc[valid,'arm_minus_baseline'],f.loc[valid,'arm_estimate']-f.loc[valid,'baseline_estimate'],rtol=1e-10,atol=1e-12),path
        checks.append(dict(file=str(path),check='every stored draw is paired algebraically',rows=len(f)))
    for path in sorted(args.output_dir.rglob('illustrative_source_variant_decisions.csv')):
        try:v=pd.read_csv(path);f=pd.read_csv(path.parent/'illustrative_choices.csv')
        except pd.errors.EmptyDataError:continue
        for (idx,model),g in v.groupby(['illustration_index','model']):
            w=g.conditional_variant_pair_weight.to_numpy(float)
            assert np.isclose(w.sum(),1),path
            for local,stored in [('pareto_correct','pareto_correct'),('pareto_reversed','pareto_reverse'),('predicted_tradeoff','pareto_tradeoff'),('predicted_tie','pareto_tie')]:
                value=np.sum(w*g[local].to_numpy(float))
                assert np.isclose(value,float(f.iloc[idx][model+'__'+stored]),atol=1e-12), (path,idx,model,local)
        checks.append(dict(file=str(path),check='all-model source-variant decisions reproduce illustration fractions; weights total one',rows=len(v)))
    for path in sorted(args.output_dir.rglob('continuous_rectangle_data.csv.gz')):
        f=pd.read_csv(path)
        models=[c[:-3] for c in f if c.endswith('__D') and f[c].notna().any()]
        for model in models:
            sub=f[f.pareto_identified]
            columns=[model+'__'+s for s in ['pareto_correct','pareto_reverse','pareto_tradeoff','pareto_tie']]
            if len(sub):assert np.allclose(sub[columns].sum(axis=1),1,atol=1e-12),path
            credit=f[model+'__D_sign_credit'].dropna()
            assert credit.between(0,1).all(),path
            assert f.loc[f.observed_D_direction.eq(0),model+'__D_sign_credit'].isna().all(),path
            if model.startswith('native_'):assert f[model+'__D_mse'].isna().all(),path
        checks.append(dict(file=str(path),check='choice census partitions one; unknown direction excluded; native squared error undefined',rows=len(f)))
    receipt=dict(status='passed',scope='existing generated tables only; complete-dataset receipts are still required',checks=checks,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.receipt.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(status='passed',checked_files=len(checks))))


if __name__=='__main__':main()
