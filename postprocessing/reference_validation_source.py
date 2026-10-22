"""Replay frozen predictions; paired, dependency-aware descriptive analysis.

This module never imports model/run modules and never fits a predictive model.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'implementation'))
from profile_information import data, evaluation, splits

BASELINE = 'global_chemistry'
METRICS = {
    'D_direction_concordance': 'D_sign_credit',
    'D_predicted_tie_fraction': 'D_tie',
    'D_mse': 'D_mse',
    'pareto_correct_fraction': 'pareto_correct',
    'pareto_reversed_fraction': 'pareto_reverse',
    'pareto_tradeoff_fraction': 'pareto_tradeoff',
    'pareto_tie_fraction': 'pareto_tie',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=str) + '\n')


def save(frame, path):
    frame.to_csv(path, index=False, compression='gzip' if str(path).endswith('.gz') else None)


def node_factor(counts, left, right):
    """One multiplicity per distinct endpoint cluster, not per cell."""
    return counts[left] * np.where(left == right, 1, counts[right])


def pooled_mean(values, mass, contexts, context_mass, n_contexts, weighting):
    valid = np.isfinite(values)
    sums = np.bincount(contexts[valid], weights=mass[valid] * values[valid], minlength=n_contexts)
    den = np.bincount(contexts[valid], weights=mass[valid], minlength=n_contexts)
    good = (den > 0) & (context_mass > 0)
    if not good.any():
        return np.nan
    if weighting == 'actual_assay_pair_equal':
        return np.average(sums[good] / den[good], weights=context_mass[good])
    return np.sum(sums[good] * context_mass[good]) / np.sum(den[good] * context_mass[good])


def metric_estimate(frame, model, metric, mass=None, context_mass=None, weighting='actual_assay_pair_equal'):
    contexts = frame['_context'].to_numpy(int)
    n = int(contexts.max()) + 1
    if mass is None:
        mass = np.ones(len(frame))
    if context_mass is None:
        context_mass = np.ones(n)
    if metric == 'symmetrized_D_correlation':
        exact = frame.observed_D_exact.to_numpy(bool)
        y = np.where(exact, frame.observed_D_lower.to_numpy(float), np.nan)
        q = frame[model + '__D'].to_numpy(float)
        q2 = frame[model + '__D2'].to_numpy(float)
        if weighting == 'raw_rectangle_equal':
            cross = pooled_mean(y*q, mass, contexts, context_mass, n, weighting)
            yy = pooled_mean(y*y, mass, contexts, context_mass, n, weighting)
            qq = pooled_mean(np.where(exact, q2, np.nan), mass, contexts, context_mass, n, weighting)
            return cross/np.sqrt(yy*qq) if yy > 0 and qq > 0 else np.nan
        valid = exact & np.isfinite(q) & np.isfinite(q2)
        numer = np.bincount(contexts[valid], weights=mass[valid]*y[valid]*q[valid], minlength=n)
        yy = np.bincount(contexts[valid], weights=mass[valid]*y[valid]**2, minlength=n)
        qq = np.bincount(contexts[valid], weights=mass[valid]*q2[valid], minlength=n)
        good = (yy > 0) & (qq > 0) & (context_mass > 0)
        return np.average(numer[good]/np.sqrt(yy[good]*qq[good]), weights=context_mass[good]) if good.any() else np.nan
    key = 'D_mse' if metric == 'D_rmse' else METRICS[metric]
    result = pooled_mean(frame[model+'__'+key].to_numpy(float), mass, contexts, context_mass, n, weighting)
    return np.sqrt(result) if metric == 'D_rmse' else result


def prepare_context(frame):
    frame = frame.copy()
    frame['_context'] = pd.factorize(pd.MultiIndex.from_frame(frame[['assay_id_a','assay_id_b']]))[0]
    return frame


def summary(frame, resamples=400, seed=20261022):
    """Arms and baseline always evaluated on the same replayed rectangle rows."""
    output, draws = [], []
    models = sorted(c[:-3] for c in frame if c.endswith('__D'))
    metrics = list(METRICS) + ['D_rmse','symmetrized_D_correlation']
    for domain, sub in [('all_source_operational', frame), ('metadata_matched_operational', frame[frame.quantitative_context_matched])]:
        if sub.empty:
            continue
        sub = prepare_context(sub)
        ctx = sub.drop_duplicates('_context').sort_values('_context')
        chem = sorted(set(sub.chemical_group_a.astype(str)) | set(sub.chemical_group_b.astype(str)))
        genes = sorted(set(sub.gene_a) | set(sub.gene_b))
        cm, gm = dict(zip(chem, range(len(chem)))), dict(zip(genes, range(len(genes))))
        ca, cb = [sub[c].astype(str).map(cm).to_numpy(int) for c in ['chemical_group_a','chemical_group_b']]
        ga, gb = [ctx[c].map(gm).to_numpy(int) for c in ['gene_a','gene_b']]
        support = dict(domain=domain,n_rectangles=len(sub),n_actual_assay_pairs=len(ctx),n_genes=len(genes),n_chemical_groups=len(chem),n_molecules=len(set(sub.molecule_a)|set(sub.molecule_b)),n_exact=int(sub.observed_D_exact.sum()),n_identified_D=int(sub.observed_D_direction.ne(0).sum()),n_pareto=int(sub.pareto_identified.sum()))
        estimates = {}
        for weighting in ['actual_assay_pair_equal','raw_rectangle_equal']:
            for model in models:
                pair_mass=(sub[model+'__D'].notna() & sub[BASELINE+'__D'].notna()).to_numpy(float) if BASELINE in models else np.ones(len(sub))
                for metric in metrics:
                    if model.startswith('native_') and metric in ['D_mse','D_rmse']:
                        continue
                    estimate = metric_estimate(sub,model,metric,mass=pair_mass,weighting=weighting)
                    estimates[weighting,model,metric] = estimate
                    base = metric_estimate(sub,BASELINE,metric,mass=pair_mass,weighting=weighting) if BASELINE in models else np.nan
                    output.append(dict(**support,n_paired_rectangles=int(pair_mass.sum()),weighting=weighting,model=model,metric=metric,estimate=estimate,baseline_estimate=base,arm_minus_baseline=estimate-base))
        rng = np.random.default_rng(seed)
        for i in range(resamples):
            cc = rng.multinomial(len(chem),np.full(len(chem),1/len(chem)))
            gc = rng.multinomial(len(genes),np.full(len(genes),1/len(genes)))
            cw, gw = node_factor(cc,ca,cb), node_factor(gc,ga,gb)
            for uncertainty,mass,context_mass in [('chemical_group_fixed_panel',cw,np.ones(len(ctx))),('gene_context_fixed_chemistry',np.ones(len(sub)),gw),('joint_chemical_gene_sensitivity',cw,gw)]:
                for weighting in ['actual_assay_pair_equal','raw_rectangle_equal']:
                    for model in models:
                        if model == BASELINE:
                            continue
                        pair_mass=mass*((sub[model+'__D'].notna() & sub[BASELINE+'__D'].notna()).to_numpy(float)) if BASELINE in models else mass
                        for metric in metrics:
                            if model.startswith('native_') and metric in ['D_mse','D_rmse']:
                                continue
                            arm = metric_estimate(sub,model,metric,pair_mass,context_mass,weighting)
                            base = metric_estimate(sub,BASELINE,metric,pair_mass,context_mass,weighting) if BASELINE in models else np.nan
                            draws.append(dict(domain=domain,weighting=weighting,model=model,metric=metric,resampling=uncertainty,replicate=i,arm_estimate=arm,baseline_estimate=base,arm_minus_baseline=arm-base))
    points = pd.DataFrame(output)
    replicates = pd.DataFrame(draws)
    if len(replicates):
        intervals=[]
        keys=['domain','weighting','model','metric','resampling']
        for key,g in replicates.groupby(keys,dropna=False):
            finite=g.arm_minus_baseline.dropna()
            intervals.append(dict(zip(keys,key),n_resamples=resamples,n_defined_resamples=len(finite),paired_delta_q025=finite.quantile(.025),paired_delta_q50=finite.quantile(.5),paired_delta_q975=finite.quantile(.975),paired_delta_resampling_sd=finite.std(ddof=1)))
        intervals=pd.DataFrame(intervals).merge(points,on=['domain','weighting','model','metric'],how='left')
    else:
        intervals=pd.DataFrame()
    return points,replicates,intervals


def descriptive_slope(x,y,w):
    valid=np.isfinite(x)&np.isfinite(y)&np.isfinite(w)&(w>0)
    x,y,w=x[valid],y[valid],w[valid]
    if not len(x):return np.nan,np.nan,0
    mx,my=np.average(x,weights=w),np.average(y,weights=w)
    den=np.sum(w*(x-mx)**2)
    slope=np.sum(w*(x-mx)*(y-my))/den if den>0 else np.nan
    return slope,my-slope*mx,len(x)


def heterogeneity(frame,d):
    models=sorted(c[:-3] for c in frame if c.endswith('__D'))
    rows=[]
    frame=prepare_context(frame)
    for kind,columns in [('actual_assay_pair',['assay_id_a','assay_id_b']),('gene_pair',['gene_pair_key']),('family_pair',['family_pair'])]:
        for key,g in frame.groupby(columns,dropna=False):
            g=prepare_context(g)
            metadata=dict(zip(columns,key if isinstance(key,tuple) else (key,)))
            for model in models:
                for metric in ['D_direction_concordance','D_rmse','pareto_correct_fraction','pareto_reversed_fraction']:
                    if model.startswith('native_') and metric=='D_rmse':continue
                    estimate=metric_estimate(g,model,metric)
                    baseline=metric_estimate(g,BASELINE,metric) if BASELINE in models else np.nan
                    rows.append(dict(grouping=kind,**metadata,model=model,metric=metric,estimate=estimate,baseline_estimate=baseline,arm_minus_baseline=estimate-baseline,n_rectangles=len(g),n_molecules=len(set(g.molecule_a)|set(g.molecule_b)),n_assay_pairs=g._context.nunique(),n_exact=int(g.observed_D_exact.sum()),n_identified_D=int(g.observed_D_direction.ne(0).sum()),n_pareto=int(g.pareto_identified.sum()),n_metadata_matched=int(g.quantitative_context_matched.sum())))
    trends=[]
    for predictor in ['observed_abs_D','chemical_tanimoto']:
        for model in models:
            if model==BASELINE:continue
            for metric in ['D_direction_concordance','pareto_correct_fraction','pareto_reversed_fraction','D_mse']:
                if model.startswith('native_') and metric=='D_mse':continue
                value=frame[model+'__'+METRICS[metric]].to_numpy()-frame[BASELINE+'__'+METRICS[metric]].to_numpy() if BASELINE in models else np.full(len(frame),np.nan)
                x=frame[predictor].to_numpy(float)
                good=np.isfinite(x)&np.isfinite(value)
                for weighting in ['actual_assay_pair_equal','raw_rectangle_equal']:
                    counts=frame.loc[good].groupby('_context')._context.transform('size')
                    weights=np.zeros(len(frame));weights[good]=1/counts.to_numpy() if weighting=='actual_assay_pair_equal' else 1
                    slope,intercept,n=descriptive_slope(x,value,weights)
                    trends.append(dict(predictor=predictor,model=model,metric=metric,weighting=weighting,slope=slope,intercept=intercept,n_defined=n,interpretation='descriptive linear trend; no mechanistic or causal identification'))
    return pd.DataFrame(rows),pd.DataFrame(trends)


def practical(frame,d):
    """A deterministic illustration per choice pattern; entire census is retained."""
    rows=[];variants=[]
    models=sorted(c[:-3] for c in frame if c.endswith('__D'))
    sorted_frame=frame.sort_values(['split','assay_id_a','assay_id_b','molecule_id_a','molecule_id_b'])
    records=d['records'].set_index(['molecule_index','local_assay_index'])
    for model in models:
        for category in ['correct','reverse','tradeoff','tie']:
            candidates=sorted_frame[sorted_frame.pareto_identified & sorted_frame[model+'__pareto_'+category].gt(0)]
            if candidates.empty:continue
            r=candidates.iloc[0]
            record=r.drop(labels=[c for c in r.index if c.startswith('_')]).to_dict()
            record.update(illustration_model=model,illustration_category=category,selection_rule='first lexical split, assay IDs, molecule IDs with positive variant-integrated mass in this category; one per category, not effect-ranked')
            for end,mol in [('a',int(r.molecule_a)),('b',int(r.molecule_b))]:
                mr=d['molecules'].iloc[mol]
                record['name_'+end]=mr['name'];record['smiles_'+end]=mr.reference_smiles
                for target,j in [('a',int(r.assay_index_a)),('b',int(r.assay_index_b))]:
                    rec=records.loc[(mol,j)]
                    for field in ['ic50_uM','source_relation','source_row_id','lower','upper','lower_open','upper_open']:
                        record['observed_'+end+'_'+target+'_'+field]=rec[field]
            rows.append(record)
    return pd.DataFrame(rows)


def practical_variants(illustrations,d,dataset_dir,regime):
    """Expose all source-row decisions underlying each selected illustration."""
    from itertools import product
    output=[];cache={}
    for illustration,r in illustrations.iterrows():
        split=str(r['split'])
        if split not in cache:
            z=np.load(dataset_dir/'core'/regime/split/'predictions.npz',allow_pickle=False)
            cache[split]={k:z[k] for k in z.files}
        pred=cache[split];test=set(map(int,pred['test_rows']))
        ma,mb=int(r.molecule_a),int(r.molecule_b)
        aa,ab=int(r.assay_index_a),int(r.assay_index_b)
        names=[k for k in pred if k not in ['test_rows','test_assays'] and (r.comparison_support=='native_matched_common' or not k.startswith('native_'))]
        for va,vb in product(d['identity_rows'][ma],d['identity_rows'][mb]):
            if int(va) not in test or int(vb) not in test:continue
            valid=all(np.isfinite(pred[n][np.ix_([va,vb],[aa,ab])]).all() for n in names)
            if not valid:continue
            for model in names:
                pa=pred[model][vb,aa]-pred[model][va,aa]
                pb=pred[model][vb,ab]-pred[model][va,ab]
                output.append(dict(illustration_index=illustration,illustration_model=r.illustration_model,split=split,comparison_support=r.comparison_support,model=model,illustration_category=r.illustration_category,molecule_id_a=r.molecule_id_a,molecule_id_b=r.molecule_id_b,assay_id_a=r.assay_id_a,assay_id_b=r.assay_id_b,source_variant_row_a=int(va),source_variant_row_b=int(vb),source_variant_weight_a=d['variant_weights'][va],source_variant_weight_b=d['variant_weights'][vb],conditional_variant_pair_weight=d['variant_weights'][va]*d['variant_weights'][vb]/r.variant_pair_mass,predicted_A_t=pred[model][va,aa],predicted_A_u=pred[model][va,ab],predicted_B_t=pred[model][vb,aa],predicted_B_u=pred[model][vb,ab],predicted_delta_t=pa,predicted_delta_u=pb,predicted_D=pa-pb,pareto_correct=(np.sign(pa)==r.delta_a_direction)&(np.sign(pb)==r.delta_b_direction),pareto_reversed=(np.sign(pa)==-r.delta_a_direction)&(np.sign(pb)==-r.delta_b_direction),predicted_tradeoff=np.sign(pa)*np.sign(pb)==1,predicted_tie=(pa==0)|(pb==0),prediction_units='native source score units' if model.startswith('native_') else 'pIC50'))
    return pd.DataFrame(output)


def replay_one(d,source,dest):
    files={name:sha(source/name) for name in ['predictions.npz','split.json','complete.json']}
    completed=json.loads((source/'complete.json').read_text())
    for name in ['predictions.npz','split.json']:
        if completed['artifact_sha256'].get(name)!=files[name]:raise AssertionError('Frozen artifact manifest mismatch: '+str(source/name))
    inputs={name:sha(d['input_dir']/name) for name in ['binding_records.csv.gz','molecules.csv','assays.csv','molecule_features.npz','datasets/'+d['dataset']+'.npz']}
    receipt_path=dest/'replay_receipt.json'
    if receipt_path.exists():
        receipt=json.loads(receipt_path.read_text())
        valid_outputs=all((dest/name).is_file() and sha(dest/name)==value for name,value in receipt.get('output_sha256',{}).items())
        if receipt.get('source_sha256')==files and receipt.get('prepared_input_sha256')==inputs and receipt.get('evaluator_sha256')==sha(ROOT/'implementation/profile_information/evaluation.py') and valid_outputs:
            return pd.read_csv(dest/'preference_changes.csv.gz',float_precision='round_trip')
    saved=np.load(source/'predictions.npz',allow_pickle=False)
    pred={k:saved[k] for k in saved.files if k not in ['test_rows','test_assays']}
    split=json.loads((source/'split.json').read_text())
    meta=dict(dataset=d['dataset'],channel=d['channel'],scope=d['scope'],policy=d['policy'],regime=split['regime'],split=split['id'],variant_set='core')
    frames=evaluation.evaluate(d,pred,saved['test_rows'],saved['test_assays'],meta,compact=False,include_records=False)
    dest.mkdir(parents=True,exist_ok=True)
    names=['record_predictions','assay_metrics','preference_changes','assay_pair_support']
    for name,frame in zip(names,frames):
        if name!='record_predictions':save(frame if len(frame.columns) else pd.DataFrame(columns=['undefined_reason']),dest/(name+'.csv.gz'))
    # Independent replay consistency with the saved compact evaluator outputs.
    regenerated=evaluation.detail_contexts(frames[2]) if len(frames[2]) else pd.DataFrame()
    saved_context=source/'preference_context_metrics.csv.gz'
    max_difference=0.0
    if saved_context.exists() and len(regenerated):
        prior=pd.read_csv(saved_context,float_precision='round_trip')
        keys=['comparison_support','model','assay_id_a','assay_id_b']
        for f in [regenerated,prior]:
            for key in keys:f[key]=f[key].astype(str)
        a=regenerated.set_index(keys).sort_index();b=prior.set_index(keys).sort_index()
        if not a.index.equals(b.index):raise AssertionError('Replay context membership differs: '+str(source))
        for c in a.select_dtypes(include=np.number):
            if c in b:
                if not np.allclose(a[c],b[c],rtol=2e-10,atol=2e-10,equal_nan=True):raise AssertionError('Replay context differs: '+c+' '+str(source))
                delta=np.abs(a[c].to_numpy(float)-b[c].to_numpy(float));max_difference=max(max_difference,float(np.nanmax(delta)) if len(delta) else 0)
    write_json(receipt_path,dict(source=str(source),source_sha256=files,prepared_input_sha256=inputs,evaluator_sha256=sha(ROOT/'implementation/profile_information/evaluation.py'),new_model_fit=False,n_rectangle_rows=len(frames[2]),compact_numeric_max_abs_difference=max_difference,output_sha256={p.name:sha(p) for p in dest.glob('*.csv.gz')}))
    return frames[2]


def process_dataset(source_root,dataset_dir,args):
    spec=dict(zip(['dataset','channel','scope','policy'],dataset_dir.name.split('__')))
    d=data.load(args.input_dir,spec['dataset'],spec['channel'] if spec['dataset']=='spd' else None,spec['scope'],spec['policy'])
    cfg=json.loads((source_root/'effective_config.json').read_text())
    scheme='gene_folds' if spec['scope']=='all_binding_genes' else cfg['core'].get('target_scheme',{}).get(spec['dataset'],'pairs')
    expected=list(splits.outer_splits(d,cfg['core']['regimes'],scheme,cfg['core'].get('target_folds',5)))
    completed=list(dataset_dir.glob('core/*/*/complete.json'))
    expected_ids={(s['regime'],s['id']) for s in expected}
    completed_ids={(p.parent.parent.name,p.parent.name) for p in completed}
    out=args.output_dir/source_root.name/dataset_dir.name
    out.mkdir(parents=True,exist_ok=True)
    write_json(out/'dataset_state.json',dict(source=str(dataset_dir),expected_splits=len(expected_ids),completed_splits=len(completed_ids),dataset_complete=expected_ids<=completed_ids,missing_splits=sorted(expected_ids-completed_ids),summary_scope='completed saved splits only; refresh after new splits finish'))
    if not (expected_ids<=completed_ids) and not args.allow_partial:
        print('SKIP incomplete',dataset_dir.name,len(completed_ids),'/',len(expected_ids),flush=True);return
    by_regime={}
    for marker in sorted(completed):
        split_source=marker.parent
        f=replay_one(d,split_source,out/'replay'/split_source.parent.name/split_source.name)
        if len(f):by_regime.setdefault(split_source.parent.name,[]).append(f)
        print('REPLAY',source_root.name,dataset_dir.name,split_source.parent.name,split_source.name,len(f),flush=True)
    for regime,frames in by_regime.items():
        frame=pd.concat(frames,ignore_index=True)
        frame['chemical_tanimoto']=d['tanimoto'][frame.molecule_a.to_numpy(int),frame.molecule_b.to_numpy(int)]
        frame['observed_abs_D']=np.where(frame.observed_D_exact,np.abs(frame.observed_D_lower),np.nan)
        families=d['assays'].set_index('assay_id').uniprot_family_level.astype(str).to_dict()
        frame['family_pair']=['|'.join(sorted([families[a],families[b]])) for a,b in zip(frame.assay_id_a,frame.assay_id_b)]
        dest=out/regime;dest.mkdir(exist_ok=True)
        for support,sub in frame.groupby('comparison_support'):
            sd=dest/support;sd.mkdir(exist_ok=True)
            points,replicates,intervals=summary(sub,args.resamples,args.seed)
            save(points,sd/'paired_estimates.csv');save(replicates,sd/'paired_resampling_draws.csv.gz');save(intervals,sd/'paired_uncertainty.csv')
            contexts,trends=heterogeneity(sub,d)
            illustrations=practical(sub,d)
            save(contexts,sd/'heterogeneity.csv');save(trends,sd/'continuous_trends.csv');save(illustrations,sd/'illustrative_choices.csv')
            save(practical_variants(illustrations,d,dataset_dir,regime),sd/'illustrative_source_variant_decisions.csv')
            # Full continuous rows permit plotting without statistical bins or effect gates.
            cols=['split','assay_id_a','assay_id_b','gene_a','gene_b','gene_pair_key','family_pair','molecule_a','molecule_b','molecule_id_a','molecule_id_b','chemical_group_a','chemical_group_b','chemical_tanimoto','observed_abs_D','observed_D_lower','observed_D_upper','observed_D_exact','observed_D_direction','pareto_identified','variant_pair_mass','quantitative_context_matched']
            cols += [c for c in sub if '__' in c]
            save(sub[cols],sd/'continuous_rectangle_data.csv.gz')
            print('SUMMARIZED',source_root.name,dataset_dir.name,regime,support,len(sub),flush=True)
    write_json(out/'postprocessing_receipt.json',dict(source=str(dataset_dir),source_state=json.loads((out/'dataset_state.json').read_text()),new_model_training=False,resamples=args.resamples,seed=args.seed,script_sha256=sha(__file__),output_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='postprocessing_receipt.json'}))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--results',type=Path,action='append',required=True)
    parser.add_argument('--input-dir',type=Path,default=ROOT/'preparation/data')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'postprocessing/outputs')
    parser.add_argument('--dataset',action='append',help='Exact dataset directory name; repeat to filter.')
    parser.add_argument('--allow-partial',action='store_true')
    parser.add_argument('--resamples',type=int,default=400)
    parser.add_argument('--seed',type=int,default=20261022)
    args=parser.parse_args()
    if args.resamples<0:parser.error('resamples must be nonnegative')
    for root in args.results:
        for folder in sorted(root.glob('*__all_records')):
            if args.dataset and folder.name not in args.dataset:continue
            if (folder/'identical_score_alias.json').exists():continue
            if list(folder.glob('core/*/*/complete.json')):process_dataset(root,folder,args)


if __name__=='__main__':
    main()
