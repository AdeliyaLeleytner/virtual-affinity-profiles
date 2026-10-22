"""Bounded secondary profile decoder/objective sensitivity; never an information bound.

Frozen baseline and best-column predictions are copied, never silently retuned.
The two inner criteria select from the SAME cached pointwise-trained candidates.
"""
from pathlib import Path
import argparse, hashlib, json, shutil
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from . import data as input_data, features, model, splits, evaluation
from .run import clean, write_json, log, tensor_slice, save_evaluation, score_signature, completed_matches, mark_complete, analysis_signature


def weighted_kernel_map(k_train, k_all, weights):
    """Exact finite RBF map centered under identity-preserving variant mass."""
    k=np.asarray(k_train,float);cross=np.asarray(k_all,float)
    w=np.asarray(weights,float);w=w/w.sum()
    mean=k@w;grand=float(w@mean)
    center=k-mean[:,None]-mean[None,:]+grand
    sw=np.sqrt(w);weighted=sw[:,None]*center*sw[None,:]
    ev,u=np.linalg.eigh((weighted+weighted.T)/2)
    tol=np.finfo(float).eps*max(k.shape)*max(1.,float(np.max(np.abs(ev))))
    if ev.min() < -100*tol:raise ValueError('Profile kernel not numerically PSD')
    keep=ev>tol;basis=sw[:,None]*u[:,keep]/np.sqrt(ev[keep])[None,:]
    crosscenter=cross-mean[None,:]-(cross@w)[:,None]+grand
    return crosscenter@basis,dict(rank=int(keep.sum()),numerical_tolerance=float(tol),minimum_eigenvalue=float(ev.min()),training_weights=w.tolist(),centering='variant-weighted training kernel mean only')


def profile_blocks(data, train_rows, columns):
    rr=np.unique(np.asarray(train_rows,int))
    if len(rr)!=len(train_rows):raise ValueError('Duplicated training source rows')
    _,linear,_,_,diag=features.score_features(data,rr,columns)
    raw=data['scores'][:,columns]
    av,adiag=features.standardize(np.isfinite(raw).astype(float),rr,data['variant_weights'][rr])
    # Same residual scaling as frozen linear arm. Availability is explicit;
    # its RBF interaction is a decoder change, though baseline already sees it.
    inputs=np.c_[linear*np.sqrt(max(1,len(columns))),av]
    kt=features.rbf(inputs[rr],inputs[rr]);ka=features.rbf(inputs,inputs[rr])
    nonlinear,kdiag=weighted_kernel_map(kt,ka,data['variant_weights'][rr])
    offdiag=kt[np.triu_indices(len(kt),1)]
    kdiag['training_offdiagonal_kernel']=dict(n_pairs=len(offdiag),minimum=float(offdiag.min()) if len(offdiag) else None,quantiles=np.quantile(offdiag,[0,.1,.25,.5,.75,.9,1]).tolist() if len(offdiag) else [],quantile_probabilities=[0,.1,.25,.5,.75,.9,1],mean=float(offdiag.mean()) if len(offdiag) else None,weighting='unweighted unique source-row pairs; kernel centering separately uses identity-preserving variant mass')
    return dict(linear=linear,rbf=nonlinear),dict(training_source_rows=rr.tolist(),training_row_molecule_indices=data['row_molecule_indices'][rr].tolist(),standardization=diag,availability=adiag,kernel=kdiag,
        rbf_definition='exp(-squared_distance/(2*d)); d=number of residual+availability inputs; no labels or held-out rows set moments, centering or basis',rbf_input_scaling='unscaled train-standardized residual and availability; dimension divisor applied only inside RBF',rbf_input_dimensions=inputs.shape[1],unique_source_rows=len(rr))


def assemble_base(data, split, columns, weight):
    known=split['known'] if 'known' in split else split['regime']=='new_molecules_known_proteins'
    x,z,q,_,diag=features.assemble(data,split['train_rows'],split['train_assays'],columns,known,'global_chemistry',descriptor_weight=weight)
    return x,z,q,diag


def fit_candidate(data, split, base, block, strength, alpha, algorithm):
    x,z,q,bdiag=base
    # Dropping a zero block, rather than adding zero coefficients, preserves
    # the baseline optimization problem and coefficients exactly.
    x=x if strength==0 else np.c_[x,strength*block]
    rr,aa=split['train_rows'],split['train_assays']
    params,diag=model.fit(x[rr],z[aa],tensor_slice(q,rr,aa),data['lower'][np.ix_(rr,aa)],data['upper'][np.ix_(rr,aa)],data['weights'][np.ix_(rr,aa)],alpha,algorithm.get('optimizer'))
    return model.prediction(params,x,z,q),params,dict(**diag,base_features=bdiag,profile_strength=float(strength),effective_profile_penalty=None if strength==0 else float(alpha/strength**2))


def candidate_objectives(data, predictions, split):
    """Call the outer evaluator verbatim, including exact decimal/strict bounds.

    Criterion: equal weight per actual assay pair, all-source operational
    identified-D direction; ties receive 1/2. Exact-only D is a saved diagnostic.
    No Pareto blend, activity threshold, or test-outcome choice.
    """
    rr,aa=split['test_rows'],split['test_assays']
    meta=dict(dataset=data['dataset'],channel=data['channel'],scope=data['scope'],policy=data['policy'],regime='inner',split='inner',variant_set='decoder_candidates')
    frames=evaluation.evaluate(data,predictions,rr,aa,meta,compact=True,include_records=False)
    contexts=frames[2]
    rows=[]
    for name,pred in predictions.items():
        c=contexts[contexts.model==name] if len(contexts) else pd.DataFrame()
        direction=np.nan;count=0;rectangles=0
        if len(c):
            grouped=c.groupby(['assay_id_a','assay_id_b'])[['D_sign_credit_sum','D_sign_credit_count']].sum()
            valid=grouped.D_sign_credit_count>0
            direction=float((grouped.loc[valid,'D_sign_credit_sum']/grouped.loc[valid,'D_sign_credit_count']).mean()) if valid.any() else np.nan
            count=int(valid.sum());rectangles=int(grouped.D_sign_credit_count.sum())
        loss=model.interval_loss(pred[np.ix_(rr,aa)],data['lower'][np.ix_(rr,aa)],data['upper'][np.ix_(rr,aa)],data['weights'][np.ix_(rr,aa)])
        rows.append(dict(candidate=name,validation_interval_mse=loss,validation_D_direction_concordance=direction,n_defined_assay_pairs=count,n_identified_D=rectangles))
    return rows,frames


def choose_candidates(trials):
    choices={}
    for decoder in ['linear','rbf']:
        for objective,field,sign in [('pointwise','validation_interval_mse',1),('direction','validation_D_direction_concordance',-1)]:
            candidates=[t for t in trials if t['decoder']==decoder and np.isfinite(t[field])]
            # Equal objectives choose the smallest block (including baseline 0).
            choices[decoder+'__'+objective]=min(candidates,key=lambda t:(sign*t[field],t['strength'])) if candidates else None
    return choices


def verify_frozen(data, split, folder, config):
    original=json.loads((folder/'split.json').read_text())
    for field in ['train_rows','test_rows','train_assays','test_assays']:
        np.testing.assert_array_equal(original[field],split[field],err_msg='Frozen split mismatch: '+field)
    completion=json.loads((folder/'complete.json').read_text())
    for name,digest in completion['artifact_sha256'].items():
        if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Frozen artifact changed: '+str(folder/name))
    tuning=json.loads((folder/'tuning.json').read_text())
    inner=splits.inner_split(data,split,config['algorithm'].get('inner_validation_fold',0));inner['known']=split['regime']=='new_molecules_known_proteins'
    for field in ['train_rows','test_rows','train_assays','test_assays']:
        np.testing.assert_array_equal(tuning['inner_split'][field],inner[field],err_msg='Frozen inner split mismatch')
    return inner,tuning


def run_split(data, split, config, frozen, folder):
    folder.mkdir(parents=True,exist_ok=True)
    if completed_matches(folder,config):return
    inner,tuning=verify_frozen(data,split,frozen,config)
    columns=np.arange(data['scores'].shape[1]);weight=tuning['chosen_descriptor_weight'];alpha=tuning['chosen_alpha']['global_chemistry']
    write_json(folder/'split.json',split);write_json(folder/'inner_split.json',inner)
    shutil.copy2(frozen/'tuning.json',folder/'frozen_tuning.json')
    write_json(folder/'baseline_inheritance.json',dict(frozen_directory=frozen,source_score_signature=score_signature(data),alpha=alpha,descriptor_weight=weight,policy='frozen pointwise-selected baseline alpha/chemistry; same alpha for every profile strength; fixed best-column comparator uses original independent chosen alpha',source_complete_sha256=hashlib.sha256((frozen/'complete.json').read_bytes()).hexdigest()))
    base=assemble_base(data,inner,columns,weight);blocks,pdiag=profile_blocks(data,inner['train_rows'],columns)
    cached={};params={};diagnostics={};trials=[]
    pred,par,diag=fit_candidate(data,inner,base,None,0.,alpha,config['algorithm'])
    cached['baseline0']=pred;params['baseline0']=par;diagnostics['baseline0']=diag
    for decoder in ['linear','rbf']:
        for strength in config['sensitivity']['strength_grid']:
            if strength==0:continue
            name=decoder+'_s'+str(strength)
            log('decoder_inner_candidate_begin',dataset=data['dataset'],regime=split['regime'],split=split['id'],candidate=name)
            pred,par,diag=fit_candidate(data,inner,base,blocks[decoder],strength,alpha,config['algorithm'])
            cached[name]=pred;params[name]=par;diagnostics[name]=diag
    objectives,frames=candidate_objectives(data,cached,inner);byname={t['candidate']:t for t in objectives}
    for decoder in ['linear','rbf']:
        for strength in config['sensitivity']['strength_grid']:
            name='baseline0' if strength==0 else decoder+'_s'+str(strength)
            trials.append(dict(**byname[name],decoder=decoder,strength=strength))
    choices=choose_candidates(trials)
    write_json(folder/'inner_tuning.json',dict(trials=trials,choices=choices,baseline_alpha=alpha,descriptor_weight=weight,profile_features=pdiag,fit_diagnostics=diagnostics,policy='same cached candidate predictions for both criteria; undefined direction leaves direction arm undefined; exact baseline selected for ties at strength 0'))
    innerfolder=folder/'inner';innerfolder.mkdir(exist_ok=True);save_evaluation(innerfolder,frames)
    np.savez_compressed(innerfolder/'candidate_coefficients.npz',**params)
    np.savez_compressed(innerfolder/'candidate_predictions.npz',**{n:p[np.ix_(inner['test_rows'],inner['test_assays'])] for n,p in cached.items()},test_rows=inner['test_rows'],test_assays=inner['test_assays'])
    del cached,params,diagnostics,blocks,base
    saved=np.load(frozen/'predictions.npz',allow_pickle=False)
    predictions={n:saved[n].copy() for n in ['global_chemistry','global_plus_selected_column','full_profile']}
    saved.close()
    # The original linear full_profile is retained as historical comparator;
    # it used its independently selected alpha, unlike this bounded grid.
    for name in ['global_chemistry','global_plus_selected_column','full_profile']:
        shutil.copy2(frozen/(name+'_model.npz'),folder/(name+'_model.npz'))
    base=assemble_base(data,split,columns,weight);blocks,pdiag=profile_blocks(data,split['train_rows'],columns)
    refits={};outerdiag={};selected={}
    for arm,choice in choices.items():
        if choice is None:
            selected[arm]=dict(undefined_reason='no defined inner direction comparisons');continue
        strength=choice['strength'];key=(choice['decoder'],strength)
        if strength==0:
            predictions[arm]=predictions['global_chemistry'].copy()
            selected[arm]=dict(**choice,reused_frozen_baseline=True)
            continue
        if key not in refits:
            pred,par,diag=fit_candidate(data,split,base,blocks[choice['decoder']],strength,alpha,config['algorithm'])
            refits[key]=pred;outerdiag[str(key)]=diag
            np.savez_compressed(folder/(choice['candidate']+'_model.npz'),coefficients=par,train_rows=split['train_rows'],train_assays=split['train_assays'],alpha=alpha,descriptor_weight=weight,profile_strength=strength)
        predictions[arm]=refits[key];selected[arm]=dict(**choice,reused_frozen_baseline=False)
    write_json(folder/'outer_diagnostics.json',dict(selected=selected,profile_features=pdiag,fits=outerdiag))
    np.savez_compressed(folder/'predictions.npz',**predictions,test_rows=split['test_rows'],test_assays=split['test_assays'])
    meta=dict(dataset=data['dataset'],channel=data['channel'],scope=data['scope'],policy=data['policy'],regime=split['regime'],split=split['id'],variant_set='decoder_sensitivity')
    frames=evaluation.evaluate(data,predictions,split['test_rows'],split['test_assays'],meta,compact=True)
    save_evaluation(folder,frames)
    mark_complete(folder,dict(status='complete',scientific_success=None,selected=selected),config)
    log('decoder_sensitivity_split_complete',**meta,selected=selected)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,default=Path(__file__).resolve().parents[1]/'config.json');parser.add_argument('--cost-only',action='store_true');args=parser.parse_args()
    config=json.loads(args.config.read_text());root=args.config.resolve().parent
    inputs=Path(config['input_dir']);inputs=inputs if inputs.is_absolute() else root/inputs
    output=Path(config['output_dir']);output=output if output.is_absolute() else root/output
    frozen=Path(config['sensitivity']['frozen_results']);output.mkdir(parents=True,exist_ok=True)
    plan=[];loaded=[]
    frozen_signature=json.loads((frozen/'analysis_signature.json').read_text())
    # Exact prepared bytes, frozen module bytes and optimizer are inherited.
    for rel,digest in frozen_signature['input_sha256'].items():
        if hashlib.sha256((inputs/rel).read_bytes()).hexdigest()!=digest:raise RuntimeError('Input bytes differ from original core: '+rel)
    package=Path(__file__).resolve().parent
    for name,digest in frozen_signature['module_sha256'].items():
        if hashlib.sha256((package/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Frozen inherited module differs: '+name)
    inherited=json.loads((frozen/'effective_config.json').read_text())
    if config['algorithm']!=inherited['algorithm']:raise RuntimeError('Baseline algorithm contract differs')
    for spec in config['datasets']:
        data=input_data.load(inputs,**spec);name='__'.join([data['dataset'],data['channel'],data['scope'],data['policy']])
        ss=list(splits.outer_splits(data,config['core']['regimes'],config['core']['target_scheme'][data['dataset']],config['core']['target_folds']))
        for split in ss:verify_frozen(data,split,frozen/name/'core'/split['regime']/split['id'],config)
        plan.append(dict(dataset=name,splits=len(ss),inner_fits_per_split=1+2*(len(config['sensitivity']['strength_grid'])-1),outer_fits_max_per_split=4,baseline_outer_fits=0))
        loaded.append((data,name,ss))
    write_json(output/'cost_plan.json',plan);log('decoder_sensitivity_cost_plan',plan=plan)
    if args.cost_only:return
    sig=analysis_signature(config,inputs)
    if (output/'analysis_signature.json').exists() and json.loads((output/'analysis_signature.json').read_text())['signature']!=sig['signature']:raise RuntimeError('Changed analysis; choose a new output directory')
    config['_analysis_signature']=sig['signature'];write_json(output/'analysis_signature.json',sig);write_json(output/'effective_config.json',config)
    with threadpool_limits(limits=config['algorithm']['numerical_threads']):
        for data,name,ss in loaded:
            for split in ss:
                log('decoder_sensitivity_split_begin',dataset=name,regime=split['regime'],split=split['id'])
                run_split(data,split,config,frozen/name/'core'/split['regime']/split['id'],output/name/'core'/split['regime']/split['id'])
            folder=output/name/'core'
            assays=[pd.read_csv(p) for p in folder.glob('*/*/assay_metrics.csv.gz')]
            prefs=[pd.read_csv(p) for p in folder.glob('*/*/preference_context_metrics.csv.gz')]
            if assays:evaluation.summarize_assays(pd.concat(assays,ignore_index=True)).to_csv(folder/'summary_assays.csv',index=False)
            if prefs:evaluation.summarize_contexts(pd.concat(prefs,ignore_index=True)).to_csv(folder/'summary_preferences.csv',index=False)
    write_json(output/'run_receipt.json',dict(status='complete',scientific_success=None,analysis_signature=sig['signature'],plan=plan))
    log('decoder_sensitivity_complete',output=output,scientific_success=None)


if __name__=='__main__':main()
