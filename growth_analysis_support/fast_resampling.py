"""Batched sparse implementation of the reference endpoint-product resampler.

Only numerical aggregation is changed. The reference remains postprocess.summary.
"""
import numpy as np
import pandas as pd
from scipy import sparse
import postprocess as reference

REFERENCE_SUMMARY = reference.summary


def _ratio(numer,denom,context_mass,weighting):
    good=(denom>0)&(context_mass>0)
    if weighting=='actual_assay_pair_equal':
        context_value=np.divide(numer,denom,out=np.zeros_like(numer),where=good)
        top=np.sum(context_value*context_mass*good,axis=0)
        bottom=np.sum(context_mass*good,axis=0)
    else:
        top=np.sum(np.where(good,numer*context_mass,0),axis=0)
        bottom=np.sum(np.where(good,denom*context_mass,0),axis=0)
    return np.divide(top,bottom,out=np.full_like(top,np.nan),where=bottom>0)


def _correlation(qy,yy,qq,context_mass,weighting):
    if weighting=='raw_rectangle_equal':
        top=np.sum(qy*context_mass,axis=0)
        denominator=np.sqrt(np.sum(yy*context_mass,axis=0)*np.sum(qq*context_mass,axis=0))
        return np.divide(top,denominator,out=np.full_like(top,np.nan),where=denominator>0)
    den=np.sqrt(yy*qq)
    good=(den>0)&(context_mass>0)
    corr=np.divide(qy,den,out=np.zeros_like(qy),where=good)
    top=np.sum(corr*context_mass*good,axis=0)
    bottom=np.sum(context_mass*good,axis=0)
    return np.divide(top,bottom,out=np.full_like(top,np.nan),where=bottom>0)


def summary(frame,resamples=400,seed=20261022,batch_size=16):
    # Remove model columns that exist only because different support strata were
    # concatenated. Their all-undefined placeholder rows have no estimand.
    missing=[c[:-3] for c in frame if c.endswith('__D') and not frame[c].notna().any()]
    frame=frame.drop(columns=[c for c in frame if any(c.startswith(m+'__') for m in missing)])
    points,_,_=REFERENCE_SUMMARY(frame,resamples=0,seed=seed)
    if resamples==0:return points,pd.DataFrame(),pd.DataFrame()
    draws=[]
    for domain,sub in [('all_source_operational',frame),('metadata_matched_operational',frame[frame.quantitative_context_matched])]:
        if sub.empty:continue
        sub=reference.prepare_context(sub)
        models=sorted(c[:-3] for c in sub if c.endswith('__D'))
        ctx=sub.drop_duplicates('_context').sort_values('_context')
        nctx=len(ctx);codes=sub._context.to_numpy(int)
        chem=sorted(set(sub.chemical_group_a.astype(str))|set(sub.chemical_group_b.astype(str)))
        genes=sorted(set(sub.gene_a)|set(sub.gene_b))
        cm,gm=dict(zip(chem,range(len(chem)))),dict(zip(genes,range(len(genes))))
        ca,cb=[sub[c].astype(str).map(cm).to_numpy(int) for c in ['chemical_group_a','chemical_group_b']]
        # Canonical chemical dyads permit one common product vector to serve
        # every context, model and metric.
        dyads=np.column_stack([np.minimum(ca,cb),np.maximum(ca,cb)])
        unique,dyad_code=np.unique(dyads,axis=0,return_inverse=True)
        ga,gb=[ctx[c].map(gm).to_numpy(int) for c in ['gene_a','gene_b']]
        rng=np.random.default_rng(seed)
        chemical_counts=[];gene_counts=[]
        for _ in range(resamples):
            chemical_counts.append(rng.multinomial(len(chem),np.full(len(chem),1/len(chem))))
            gene_counts.append(rng.multinomial(len(genes),np.full(len(genes),1/len(genes))))
        chemical_counts=np.array(chemical_counts);gene_counts=np.array(gene_counts)
        gy=gene_counts[:,ga]*gene_counts[:,gb]
        for model in models:
            if model==reference.BASELINE:continue
            pair=(sub[model+'__D'].notna()&sub[reference.BASELINE+'__D'].notna()).to_numpy() if reference.BASELINE in models else np.ones(len(sub),bool)
            names=[model,reference.BASELINE] if reference.BASELINE in models else [model]
            matrices={};fixed={}
            def register(key,value):
                value=np.asarray(value,float)
                finite=np.isfinite(value)&pair
                matrix=sparse.coo_matrix((np.where(finite,value,0),(codes,dyad_code)),shape=(nctx,len(unique))).tocsr()
                matrix.eliminate_zeros()
                matrices[key]=matrix;fixed[key]=np.asarray(matrix.sum(axis=1))
            for name in names:
                for metric in reference.METRICS.values():
                    if model.startswith('native_') and metric=='D_mse':continue
                    value=sub[name+'__'+metric].to_numpy(float)
                    register((name,metric,'num'),value)
                    register((name,metric,'den'),np.isfinite(value).astype(float))
                exact=sub.observed_D_exact.to_numpy(bool)
                y=np.where(exact,sub.observed_D_lower.to_numpy(float),np.nan)
                q=sub[name+'__D'].to_numpy(float);q2=sub[name+'__D2'].to_numpy(float)
                register((name,'corr','qy'),y*q)
                register((name,'corr','yy'),y*y)
                register((name,'corr','qq'),np.where(exact,q2,np.nan))
            metrics=list(reference.METRICS)+['D_rmse','symmetrized_D_correlation']
            if model.startswith('native_'):metrics=[m for m in metrics if m not in ['D_mse','D_rmse']]
            for start in range(0,resamples,batch_size):
                stop=min(resamples,start+batch_size);nb=stop-start
                counts=chemical_counts[start:stop]
                cw=counts[:,unique[:,0]]*np.where(unique[:,0]==unique[:,1],1,counts[:,unique[:,1]])
                chemistry={key:matrix@cw.T for key,matrix in matrices.items()}
                fixed_batch={key:np.broadcast_to(value,(nctx,nb)) for key,value in fixed.items()}
                gw=gy[start:stop].T
                for uncertainty,stats,context_mass in [('chemical_group_fixed_panel',chemistry,np.ones((nctx,nb))),('gene_context_fixed_chemistry',fixed_batch,gw),('joint_chemical_gene_sensitivity',chemistry,gw)]:
                    for weighting in ['actual_assay_pair_equal','raw_rectangle_equal']:
                        for metric in metrics:
                            values={}
                            for name in names:
                                if metric=='symmetrized_D_correlation':
                                    val=_correlation(stats[name,'corr','qy'],stats[name,'corr','yy'],stats[name,'corr','qq'],context_mass,weighting)
                                else:
                                    key='D_mse' if metric=='D_rmse' else reference.METRICS[metric]
                                    val=_ratio(stats[name,key,'num'],stats[name,key,'den'],context_mass,weighting)
                                    if metric=='D_rmse':val=np.sqrt(val)
                                values[name]=val
                            arm=values[model];baseline=values.get(reference.BASELINE,np.full(nb,np.nan))
                            for i in range(nb):draws.append(dict(domain=domain,weighting=weighting,model=model,metric=metric,resampling=uncertainty,replicate=start+i,arm_estimate=arm[i],baseline_estimate=baseline[i],arm_minus_baseline=arm[i]-baseline[i]))
    replicates=pd.DataFrame(draws)
    intervals=[]
    keys=['domain','weighting','model','metric','resampling']
    for key,g in replicates.groupby(keys,dropna=False):
        finite=g.arm_minus_baseline.dropna()
        intervals.append(dict(zip(keys,key),n_resamples=resamples,n_defined_resamples=len(finite),paired_delta_q025=finite.quantile(.025),paired_delta_q50=finite.quantile(.5),paired_delta_q975=finite.quantile(.975),paired_delta_resampling_sd=finite.std(ddof=1)))
    intervals=pd.DataFrame(intervals).merge(points,on=['domain','weighting','model','metric'],how='left')
    return points,replicates,intervals
