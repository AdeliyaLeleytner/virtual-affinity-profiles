"""Staged, config-driven study runner. Never picks scientific success cutoffs."""
from pathlib import Path
import argparse,gzip,hashlib,json,os,platform,time
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from . import data as input_data
from . import model,features,splits,evaluation


def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.ndarray):return clean(value.tolist())
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    if isinstance(value,Path):return str(value)
    return value


def write_json(path,value):
    Path(path).write_text(json.dumps(clean(value),indent=2,ensure_ascii=False,allow_nan=False)+'\n')


def log(event,**kw):print(json.dumps(clean(dict(event=event,**kw)),allow_nan=False),flush=True)


def tensor_slice(q,rr,aa):return q[rr][:,aa,:]


def fit_predict(data,split,cols,variant,alpha,descriptor_weight,algorithm,cache=None,initial=None):
    rr,aa=split['train_rows'],split['train_assays'];known=split.get('regime')=='new_molecules_known_proteins'
    if 'known' in split:known=split['known']
    if not len(rr) or not len(aa) or not (data['weights'][np.ix_(rr,aa)]>0).any():
        return None,dict(undefined_reason='no observed training intervals'),None
    observed=data['weights'][np.ix_(rr,aa)]>0
    lo=data['lower'][np.ix_(rr,aa)][observed];hi=data['upper'][np.ix_(rr,aa)][observed]
    if np.isneginf(lo).all() or np.isposinf(hi).all():
        return None,dict(undefined_reason='all training observations bound the same side; unpenalized global intercept not identified',observed_cells=int(observed.sum())),None
    selection=None;selection_report=[]
    if variant in ('global_plus_selected_column','chemistry_plus_selected_score'):selection,selection_report=features.select_columns(data,rr,aa,cols,known)
    x,z,q,native,fdiag=features.assemble(data,rr,aa,cols,known,variant,cache,descriptor_weight,selection)
    params,diag=model.fit(x[rr],z[aa],tensor_slice(q,rr,aa),data['lower'][np.ix_(rr,aa)],
                          data['upper'][np.ix_(rr,aa)],data['weights'][np.ix_(rr,aa)],alpha,
                          algorithm.get('optimizer'),initial)
    prediction=model.prediction(params,x,z,q)
    return prediction,dict(**diag,features=fdiag,column_selection=selection_report),params


def tuning(data,outer,cols,variants,algorithm,folder):
    inner=splits.inner_split(data,outer,algorithm.get('inner_validation_fold',0))
    inner['known']=outer['regime']=='new_molecules_known_proteins'
    trials=[];alphas=algorithm['alpha_grid'];weights=algorithm.get('descriptor_weights',[0.0,0.5])
    if any(not len(inner[k]) for k in ['train_rows','test_rows','train_assays','test_assays']):
        write_json(folder/'tuning.json',dict(undefined_reason='empty eligible inner training or validation axis',inner_split=inner))
        return None,None
    cache_by_weight={}
    for weight in weights:
        xc,cd=features.chemistry(data,inner['train_rows'],weight)
        z,zd=features.proteins(data,inner['train_assays'],inner['known']);cache=(xc,cd,z,zd);cache_by_weight[weight]=cache
        for alpha in alphas:
            pred,diag,_=fit_predict(data,inner,cols,'global_chemistry',alpha,weight,algorithm,cache)
            loss=np.nan if pred is None else model.interval_loss(pred[np.ix_(inner['test_rows'],inner['test_assays'])],data['lower'][np.ix_(inner['test_rows'],inner['test_assays'])],data['upper'][np.ix_(inner['test_rows'],inner['test_assays'])],data['weights'][np.ix_(inner['test_rows'],inner['test_assays'])])
            trials.append(dict(variant='global_chemistry',alpha=alpha,descriptor_weight=weight,validation_interval_mse=loss,fit=diag))
    finite=[t for t in trials if np.isfinite(t['validation_interval_mse'])]
    if not finite:
        write_json(folder/'tuning.json',dict(trials=trials,undefined_reason='no finite inner-validation interval loss'))
        return None,None
    best=min(finite,key=lambda t:(t['validation_interval_mse'],t['descriptor_weight'],t['alpha']))
    weight=best['descriptor_weight'];choices={'global_chemistry':best['alpha']}
    for variant in variants:
        if variant=='global_chemistry':continue
        for alpha in alphas:
            pred,diag,_=fit_predict(data,inner,cols,variant,alpha,weight,algorithm,cache_by_weight[weight])
            loss=np.nan if pred is None else model.interval_loss(pred[np.ix_(inner['test_rows'],inner['test_assays'])],data['lower'][np.ix_(inner['test_rows'],inner['test_assays'])],data['upper'][np.ix_(inner['test_rows'],inner['test_assays'])],data['weights'][np.ix_(inner['test_rows'],inner['test_assays'])])
            trials.append(dict(variant=variant,alpha=alpha,descriptor_weight=weight,validation_interval_mse=loss,fit=diag))
        eligible=[t for t in trials if t['variant']==variant and np.isfinite(t['validation_interval_mse'])]
        if eligible:choices[variant]=min(eligible,key=lambda t:(t['validation_interval_mse'],t['alpha']))['alpha']
    write_json(folder/'tuning.json',dict(trials=trials,chosen_alpha=choices,chosen_descriptor_weight=weight,
               inner_split=inner,policy='one inherited inner molecular fold or a simultaneous held-out inner gene block; no outer test outcomes'))
    return choices,weight


def save_evaluation(folder,frames):
    preference_name='preference_context_metrics' if 'preference_context_schema' in frames[2] else 'preference_changes'
    for name,frame in zip(['record_predictions','assay_metrics',preference_name,'assay_pair_support'],frames):
        if not len(frame.columns):frame=pd.DataFrame(columns=['undefined_reason'])
        frame.to_csv(folder/(name+'.csv.gz'),index=False,compression='gzip')


def native_scores(data):
    native=np.full(data['lower'].shape,np.nan);lookup={g:i for i,g in enumerate(data['column_genes'])}
    for j,g in enumerate(data['assays'].gene):
        if g in lookup:native[:,j]=data['scores'][:,lookup[g]]
    return native


def core_dataset(data,config,output):
    algorithm=config['algorithm'];variants=config['core']['variants'];cols=np.arange(data['scores'].shape[1])
    scheme=config['core'].get('target_scheme',{}).get(data['dataset'],'pairs')
    if data['scope']=='all_binding_genes':scheme='gene_folds'
    all_assays=[];fit_counts=0
    for split in splits.outer_splits(data,config['core']['regimes'],scheme,config['core'].get('target_folds',5)):
        folder=output/split['regime']/split['id'];folder.mkdir(parents=True,exist_ok=True)
        if completed_matches(folder,config):
            log('resume_skip',folder=folder);continue
        write_json(folder/'split.json',split)
        if any(not len(split[k]) for k in ['train_rows','test_rows','train_assays','test_assays']):
            write_json(folder/'undefined.json',dict(reason='empty training or evaluation axis',split=split));continue
        log('core_split_begin',dataset=data['dataset'],channel=data['channel'],scope=data['scope'],regime=split['regime'],split=split['id'])
        choices,weight=tuning(data,split,cols,variants,algorithm,folder)
        if choices is None:
            native=native_scores(data)
            meta=dict(dataset=data['dataset'],channel=data['channel'],scope=data['scope'],policy=data['policy'],regime=split['regime'],split=split['id'],variant_set='core')
            frames=evaluation.evaluate(data,{'native_score':native},split['test_rows'],split['test_assays'],meta,compact=True)
            save_evaluation(folder,frames)
            np.savez_compressed(folder/'predictions.npz',native_score=native,test_rows=split['test_rows'],test_assays=split['test_assays'])
            mark_complete(folder,dict(status='complete_native_only',undefined_learned_reason='see tuning.json',scientific_success=None),config)
            continue
        xc,cd=features.chemistry(data,split['train_rows'],weight);z,zd=features.proteins(data,split['train_assays'],split['regime']=='new_molecules_known_proteins')
        cache=(xc,cd,z,zd);predictions={};diagnostics={}
        for variant in variants:
            if variant not in choices:continue
            pred,diag,params=fit_predict(data,split,cols,variant,choices[variant],weight,algorithm,cache)
            diagnostics[variant]=diag
            if pred is not None:
                predictions[variant]=pred
                np.savez_compressed(folder/(variant+'_model.npz'),coefficients=params,train_rows=split['train_rows'],train_assays=split['train_assays'],alpha=choices[variant],descriptor_weight=weight)
                fit_counts+=1
        # Native score is never calibrated/fitted or named as trained output.
        native=native_scores(data)
        predictions['native_score']=native
        np.savez_compressed(folder/'predictions.npz',**predictions,test_rows=split['test_rows'],test_assays=split['test_assays'])
        meta=dict(dataset=data['dataset'],channel=data['channel'],scope=data['scope'],policy=data['policy'],regime=split['regime'],split=split['id'],variant_set='core')
        frames=evaluation.evaluate(data,predictions,split['test_rows'],split['test_assays'],meta,compact=not config['core'].get('save_preference_details',False))
        save_evaluation(folder,frames);all_assays.append(frames[1]);write_json(folder/'fit_diagnostics.json',diagnostics)
        mark_complete(folder,dict(status='complete',fitted_variants=list(choices),scientific_success=None),config)
        log('core_split_complete',**meta,preference_changes=len(frames[2]),evaluated_assay_rows=len(frames[1]))
    # Include prior completed splits on a resumed run.
    tables=[]
    for path in output.glob('*/*/assay_metrics.csv.gz'):
        tables.append(pd.read_csv(path))
    if tables:
        table=pd.concat(tables,ignore_index=True);evaluation.summarize_assays(table).to_csv(output/'summary_assays.csv',index=False)
    preferences=[pd.read_csv(p) for p in output.glob('*/*/preference_context_metrics.csv.gz')]
    if preferences:evaluation.summarize_contexts(pd.concat(preferences,ignore_index=True)).to_csv(output/'summary_preferences.csv',index=False)
    return fit_counts


def panel_orders(data,seeds):
    # Class annotations are used only for ordering columns, never labels.
    path=data['input_dir']/'datasets'/(data['dataset']+'_columns.csv')
    if path.exists():
        frame=pd.read_csv(path);name='broad_class' if 'broad_class' in frame else 'mechanistic_class'
        classes=frame[name].fillna('unannotated').astype(str).to_numpy() if name in frame else np.repeat('unannotated',len(data['column_genes']))
    else:classes=np.repeat('unannotated',len(data['column_genes']))
    if len(classes)!=len(data['column_genes']):raise ValueError('Column class metadata mismatch')
    for seed in seeds:
        rng=np.random.default_rng(seed);random=rng.permutation(len(classes))
        yield 'random_'+str(seed),random,dict(policy='uniform permutation',seed=seed)
        buckets={c:list(rng.permutation(np.flatnonzero(classes==c))) for c in sorted(set(classes))}
        order=[];class_order=list(rng.permutation(list(buckets)))
        while any(buckets.values()):
            for c in class_order:
                if buckets[c]:order.append(buckets[c].pop())
        yield 'family_stratified_'+str(seed),np.array(order,int),dict(policy='round-robin broad class; random within class and class order',seed=seed)


def growth_dataset(data,config,output):
    algorithm=config['algorithm'];growth=config['growth'];rows=[]
    for split in splits.outer_splits(data,['new_molecules_known_proteins']):
        basefolder=output/split['id'];basefolder.mkdir(parents=True,exist_ok=True)
        # Hyperparameters are selected WITHOUT ANY panel-score column.
        choices,weight=tuning(data,split,np.empty(0,int),['global_chemistry'],algorithm,basefolder)
        if choices is None:continue
        alpha=choices['global_chemistry'];xc,cd=features.chemistry(data,split['train_rows'],weight)
        z,zd=features.proteins(data,split['train_assays'],True);cache=(xc,cd,z,zd)
        comparator_folder=basefolder/'query_score_comparators';comparator_folder.mkdir(exist_ok=True)
        if not completed_matches(comparator_folder,config):
            predictions={};diagnostics={};all_columns=np.arange(data['scores'].shape[1])
            for variant in ['chemistry_plus_target_score','chemistry_plus_selected_score']:
                pred,diag,_=fit_predict(data,split,all_columns,variant,alpha,weight,algorithm,cache)
                if pred is not None:predictions[variant]=pred
                diagnostics[variant]=diag
            meta=dict(dataset=data['dataset'],channel=data['channel'],scope=data['scope'],policy=data['policy'],regime=split['regime'],split=split['id'],variant_set='query_score_comparators')
            frames=evaluation.evaluate(data,predictions,split['test_rows'],split['test_assays'],meta,compact=True,include_records=False)
            save_evaluation(comparator_folder,frames)
            rr,aa=split['test_rows'],split['test_assays']
            np.savez_compressed(comparator_folder/'heldout_predictions.npz',**{name:pred[np.ix_(rr,aa)] for name,pred in predictions.items()},test_rows=rr,test_assays=aa)
            write_json(comparator_folder/'split.json',split);write_json(comparator_folder/'fit_diagnostics.json',diagnostics)
            write_json(comparator_folder/'score_acquisition.json',dict(
                chemistry_plus_target_score=dict(training='matching native column for each labelled assay',deployment_per_single_target=1,deployment_per_distinct_target_pair=2,deployment_whole_profile='union of queried genes; not one shared panel column'),
                chemistry_plus_selected_score=dict(training_candidate_columns=len(all_columns),selection='training actual-assay intervals only',deployment_per_single_target=1,deployment_per_distinct_target_pair='one or two selected columns',deployment_whole_profile='union of assay-selected columns'),
                full_panel_mean_used=False,plotting='separate comparator; never place at k=1 of a shared fixed panel'))
            mark_complete(comparator_folder,dict(status='complete'),config)
        for name,order,orderdiag in panel_orders(data,growth['seeds']):
            folder=basefolder/name;folder.mkdir(exist_ok=True)
            write_json(folder/'order.json',dict(columns=order,column_genes=data['column_genes'][order],**orderdiag))
            for k in range(1,len(order)+1):
                cell=folder/('k'+str(k));cell.mkdir(exist_ok=True)
                if completed_matches(cell,config):continue
                cols=order[:k];predictions={};diags={}
                for variant in growth['variants']:
                    pred,diag,_=fit_predict(data,split,cols,variant,alpha,weight,algorithm,cache)
                    if pred is not None:predictions[variant]=pred
                    diags[variant]=diag
                meta=dict(dataset=data['dataset'],channel=data['channel'],scope=data['scope'],policy=data['policy'],regime=split['regime'],split=split['id'],variant_set='growth',order=name,k=k)
                frames=evaluation.evaluate(data,predictions,split['test_rows'],split['test_assays'],meta,compact=True,include_records=False)
                save_evaluation(cell,frames);write_json(cell/'fit_diagnostics.json',diags)
                rr,aa=split['test_rows'],split['test_assays']
                np.savez_compressed(cell/'heldout_predictions.npz',**{name:pred[np.ix_(rr,aa)] for name,pred in predictions.items()},test_rows=rr,test_assays=aa)
                write_json(cell/'split.json',split)
                mark_complete(cell,dict(status='complete',alpha=alpha,descriptor_weight=weight,purchased_columns=cols,full_panel_mean_used=False,hyperparameters_from='chemical-only inner training'),config)
                log('growth_k_complete',**meta)
    paths=sorted(output.glob('*/*/k*/assay_metrics.csv.gz'))
    if paths:
        with gzip.open(output/'all_k_assay_metrics.csv.gz','wt') as stream:
            for i,path in enumerate(paths):pd.read_csv(path).to_csv(stream,index=False,header=i==0)
    curve_keys=sorted({(p.parent.parent.name,p.parent.name) for p in output.glob('*/*/k*/preference_context_metrics.csv.gz')})
    summaries=[];assay_summaries=[]
    for order,k in curve_keys:
        pp=sorted(output.glob('*/'+order+'/'+k+'/preference_context_metrics.csv.gz'))
        frames=[pd.read_csv(p) for p in pp];frames=[f for f in frames if len(f)]
        if frames:summaries.append(evaluation.summarize_contexts(pd.concat(frames,ignore_index=True)))
        ap=sorted(output.glob('*/'+order+'/'+k+'/assay_metrics.csv.gz'))
        ff=[pd.read_csv(p) for p in ap];ff=[f for f in ff if len(f)]
        if ff:assay_summaries.append(evaluation.summarize_assays(pd.concat(ff,ignore_index=True)))
    if summaries:pd.concat(summaries,ignore_index=True).to_csv(output/'summary_preferences.csv',index=False)
    if assay_summaries:pd.concat(assay_summaries,ignore_index=True).to_csv(output/'summary_assays.csv',index=False)
    comparator=[pd.read_csv(p) for p in output.glob('*/query_score_comparators/preference_context_metrics.csv.gz')]
    comparator=[f for f in comparator if len(f)]
    if comparator:evaluation.summarize_contexts(pd.concat(comparator,ignore_index=True)).to_csv(output/'query_score_comparator_summary.csv',index=False)



def cost_plan(config,input_dir):
    rows=[]
    for spec in config['datasets']:
        data=input_data.load(input_dir,spec['dataset'],spec.get('channel'),spec.get('scope','native_genes'),spec.get('policy','all_records'))
        scheme=config['core'].get('target_scheme',{}).get(data['dataset'],'pairs')
        if data['scope']=='all_binding_genes':scheme='gene_folds'
        ns=sum(1 for _ in splits.outer_splits(data,config['core']['regimes'],scheme,config['core'].get('target_folds',5)))
        na=len(config['algorithm']['alpha_grid']);nw=len(config['algorithm']['descriptor_weights']);nv=len(config['core']['variants'])
        growth=bool(spec.get('growth',False) and config.get('growth',{}).get('enabled',False))
        nk=data['scores'].shape[1];ng=5*2*len(config.get('growth',{}).get('seeds',[]))*nk*len(config.get('growth',{}).get('variants',[])) if growth else 0
        rows.append(dict(**input_data.describe(data),core_outer_splits=ns,core_fit_count_upper_bound=ns*(na*nw+na*(nv-1)+nv),growth_fit_count_without_base_tuning=ng,target_scheme=scheme))
    return rows


def score_signature(data):
    digest=hashlib.sha256()
    score=np.asarray(data['scores'],dtype='<f8').copy()
    mask=np.isfinite(score);score[~mask]=0
    for value in [score,mask,data['row_molecule_indices'],data['variant_weights']]:
        digest.update(np.ascontiguousarray(value).tobytes())
    digest.update(json.dumps(data['column_genes'].tolist()).encode())
    digest.update((data['scope']+'|'+data['policy']).encode())
    return digest.hexdigest()


def completed_matches(folder,config):
    path=folder/'complete.json'
    if not path.exists():return False
    previous=json.loads(path.read_text())
    if previous.get('analysis_signature')!=config['_analysis_signature']:
        raise RuntimeError('Completed output signature mismatch: '+str(folder))
    artifacts=previous.get('artifact_sha256')
    if not artifacts:raise RuntimeError('Completion has no artifact manifest: '+str(folder))
    for name,expected in artifacts.items():
        p=folder/name
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
            raise RuntimeError('Completed artifact missing/changed: '+str(p))
    return config.get('resume',True)


def mark_complete(folder,payload,config):
    artifacts={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.iterdir()) if p.is_file() and p.name!='complete.json'}
    write_json(folder/'complete.json',dict(**payload,analysis_signature=config['_analysis_signature'],artifact_sha256=artifacts))


def analysis_signature(config,inputs):
    files={str(p.relative_to(inputs)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(inputs.rglob('*')) if p.is_file()}
    package=Path(__file__).resolve().parent
    modules={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(package.glob('*.py'))}
    body=dict(config={k:v for k,v in config.items() if not k.startswith('_')},input_sha256=files,module_sha256=modules)
    signature=hashlib.sha256(json.dumps(body,sort_keys=True).encode()).hexdigest()
    return dict(signature=signature,**body)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,default=Path(__file__).resolve().parents[1]/'config.json');parser.add_argument('--cost-only',action='store_true');args=parser.parse_args()
    config=json.loads(args.config.read_text());root=args.config.resolve().parent
    inputs=Path(config['input_dir']);inputs=inputs if inputs.is_absolute() else root/inputs
    output=Path(config['output_dir']);output=output if output.is_absolute() else root/output;output.mkdir(parents=True,exist_ok=True)
    plan=cost_plan(config,inputs);write_json(output/'cost_plan.json',plan);log('cost_plan',datasets=plan)
    if args.cost_only:return
    signature=analysis_signature(config,inputs);signature_path=output/'analysis_signature.json'
    if signature_path.exists() and json.loads(signature_path.read_text())['signature']!=signature['signature']:
        raise RuntimeError('Analysis configuration/data/code changed. Use a new output directory; no mixed partial reuse.')
    write_json(signature_path,signature);config['_analysis_signature']=signature['signature'];write_json(output/'effective_config.json',config)
    completed={};aliases=[]
    with threadpool_limits(limits=config['algorithm'].get('numerical_threads',1)):
        for spec in config['datasets']:
            data=input_data.load(inputs,spec['dataset'],spec.get('channel'),spec.get('scope','native_genes'),spec.get('policy','all_records'))
            name='__'.join([data['dataset'],data['channel'],data['scope'],data['policy']]);folder=output/name;folder.mkdir(exist_ok=True)
            write_json(folder/'input_description.json',input_data.describe(data));score_key=score_signature(data)
            work=dict(dataset=data['dataset'],scope=data['scope'],policy=data['policy'],core=config['core'],growth_enabled=bool(config.get('growth',{}).get('enabled',False) and spec.get('growth',False)))
            work_key=hashlib.sha256(json.dumps(work,sort_keys=True).encode()).hexdigest();key=(score_key,work_key)
            if key in completed:
                alias=dict(dataset=name,alias_of=completed[key],score_signature=score_key,effective_work_signature=work_key,analysis_signature=config['_analysis_signature'],reason='identical finite score array, missingness, source identity rows, variant weights and target order; not independent prediction evidence')
                write_json(folder/'identical_score_alias.json',alias);aliases.append(alias);log('identical_score_alias',**alias);continue
            log('dataset_begin',**input_data.describe(data))
            if config['core'].get('enabled',True):core_dataset(data,config,folder/'core')
            if config.get('growth',{}).get('enabled',False) and spec.get('growth',False):growth_dataset(data,config,folder/'growth')
            completed[key]=name;log('dataset_complete',dataset=name)
    write_json(output/'identical_score_aliases.json',aliases)
    write_json(output/'run_receipt.json',dict(status='complete',scientific_success=None,python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,algorithm=config['algorithm'],analysis_signature=config['_analysis_signature'],unique_score_datasets=len(completed),aliased_score_datasets=len(aliases)))


if __name__=='__main__':main()
