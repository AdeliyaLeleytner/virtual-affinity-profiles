"""Paired association on common-defined contexts, including each resampling draw.

This sidecar does not modify the frozen evaluator or reference postprocessing.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd
from scipy import sparse
from fast_resampling import _correlation
from postprocess import prepare_context


def estimate_pair(a,b,context_mass,weighting):
    common=(a['yy']>0)&(a['qq']>0)&(b['yy']>0)&(b['qq']>0)
    # The raw estimand pools moments over the same finite rectangles. Individual
    # zero-norm assay pairs remain in this pooled estimand; global norms define it.
    mass=context_mass*common if weighting=='actual_assay_pair_equal' else context_mass
    av=_correlation(a['qy'],a['yy'],a['qq'],mass,weighting)
    bv=_correlation(b['qy'],b['yy'],b['qq'],mass,weighting)
    return av,bv,common


def analyse(frame,resamples=400,seed=20261022,batch_size=16):
    points=[];draws=[]
    models=sorted(c[:-3] for c in frame if c.endswith('__D') and frame[c].notna().any())
    if 'global_chemistry' not in models:return pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
    for domain,sub in [('all_source_operational',frame),('metadata_matched_operational',frame[frame.quantitative_context_matched])]:
        if sub.empty:continue
        sub=prepare_context(sub);ctx=sub.drop_duplicates('_context').sort_values('_context');nctx=len(ctx)
        chem=sorted(set(sub.chemical_group_a.astype(str))|set(sub.chemical_group_b.astype(str)));cm=dict(zip(chem,range(len(chem))))
        genes=sorted(set(sub.gene_a)|set(sub.gene_b));gm=dict(zip(genes,range(len(genes))))
        ca,cb=[sub[c].astype(str).map(cm).to_numpy(int) for c in ['chemical_group_a','chemical_group_b']]
        unique,dyad=np.unique(np.sort(np.column_stack([ca,cb]),axis=1),axis=0,return_inverse=True)
        ga,gb=[ctx[c].map(gm).to_numpy(int) for c in ['gene_a','gene_b']]
        rng=np.random.default_rng(seed);cc=[];gc=[]
        for _ in range(resamples):
            cc.append(rng.multinomial(len(chem),np.full(len(chem),1/len(chem))))
            gc.append(rng.multinomial(len(genes),np.full(len(genes),1/len(genes))))
        cc=np.asarray(cc);gc=np.asarray(gc)
        for model in models:
            pair=sub[model+'__D'].notna()&sub.global_chemistry__D.notna()
            exact=pair&sub.observed_D_exact
            y=sub.observed_D_lower.to_numpy(float)
            matrices={};fixed={}
            for name in [model,'global_chemistry']:
                q=sub[name+'__D'].to_numpy(float);q2=sub[name+'__D2'].to_numpy(float)
                for metric,values in [('qy',y*q),('yy',y*y),('qq',q2)]:
                    values=np.where(exact&np.isfinite(values),values,0)
                    mat=sparse.coo_matrix((values,(sub._context.to_numpy(int),dyad)),shape=(nctx,len(unique))).tocsr();mat.eliminate_zeros()
                    matrices[name,metric]=mat;fixed[name,metric]=np.asarray(mat.sum(axis=1))
            a={k:fixed[model,k] for k in ['qy','yy','qq']};b={k:fixed['global_chemistry',k] for k in ['qy','yy','qq']}
            for weighting in ['actual_assay_pair_equal','raw_rectangle_equal']:
                av,bv,common=estimate_pair(a,b,np.ones((nctx,1)),weighting)
                points.append(dict(domain=domain,weighting=weighting,model=model,metric='symmetrized_D_correlation',estimate=float(av[0]),baseline_estimate=float(bv[0]),arm_minus_baseline=float(av[0]-bv[0]),n_common_defined_assay_pairs=int(common.sum()),n_arm_defined_assay_pairs=int(((a['yy']>0)&(a['qq']>0)).sum()),n_baseline_defined_assay_pairs=int(((b['yy']>0)&(b['qq']>0)).sum()),n_paired_exact_rectangles=int(exact.sum()),definition='actual-pair macro uses common-defined contexts; raw pools the same finite rectangle moments'))
            for start in range(0,resamples,batch_size):
                stop=min(resamples,start+batch_size);nb=stop-start;c=cc[start:stop]
                cw=c[:,unique[:,0]]*np.where(unique[:,0]==unique[:,1],1,c[:,unique[:,1]])
                chem_stats={key:mat@cw.T for key,mat in matrices.items()}
                fixed_stats={key:np.broadcast_to(value,(nctx,nb)) for key,value in fixed.items()}
                gw=(gc[start:stop,ga]*gc[start:stop,gb]).T
                for scheme,stats,context_mass in [('chemical_group_fixed_panel',chem_stats,np.ones((nctx,nb))),('gene_context_fixed_chemistry',fixed_stats,gw),('joint_chemical_gene_sensitivity',chem_stats,gw)]:
                    a={k:stats[model,k] for k in ['qy','yy','qq']};b={k:stats['global_chemistry',k] for k in ['qy','yy','qq']}
                    for weighting in ['actual_assay_pair_equal','raw_rectangle_equal']:
                        av,bv,common=estimate_pair(a,b,context_mass,weighting)
                        for i in range(nb):draws.append(dict(domain=domain,weighting=weighting,model=model,metric='symmetrized_D_correlation',resampling=scheme,replicate=start+i,arm_estimate=av[i],baseline_estimate=bv[i],arm_minus_baseline=av[i]-bv[i],n_sampled_common_defined_assay_pairs=int((common[:,i]&(context_mass[:,i]>0)).sum())))
    points=pd.DataFrame(points);draws=pd.DataFrame(draws);intervals=[]
    if len(draws):
        keys=['domain','weighting','model','metric','resampling']
        for key,g in draws.groupby(keys):
            f=g.arm_minus_baseline.dropna()
            intervals.append(dict(zip(keys,key),n_resamples=resamples,n_defined_resamples=len(f),paired_delta_q025=f.quantile(.025),paired_delta_q50=f.quantile(.5),paired_delta_q975=f.quantile(.975),paired_delta_resampling_sd=f.std(ddof=1)))
    return points,draws,pd.DataFrame(intervals)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--save-dir',type=Path,required=True);p.add_argument('--resamples',type=int,default=400);p.add_argument('--seed',type=int,default=20261022);args=p.parse_args()
    for source in sorted(args.output_dir.rglob('continuous_rectangle_data.csv.gz')):
        f=pd.read_csv(source,float_precision='round_trip');point,draws,intervals=analyse(f,args.resamples,args.seed)
        dest=args.save_dir/source.parent.relative_to(args.output_dir);dest.mkdir(parents=True,exist_ok=True)
        point.to_csv(dest/'paired_association.csv',index=False);draws.to_csv(dest/'paired_association_draws.csv.gz',index=False);intervals.to_csv(dest/'paired_association_uncertainty.csv',index=False)
        receipt=dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seed=args.seed,resamples=args.resamples,new_training=False,policy='Common-defined assay contexts are intersected separately for each arm-baseline pair and each draw in primary macro association; raw estimand is pooled moments on common rectangle support.')
        (dest/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print('ASSOCIATION',str(source),flush=True)


if __name__=='__main__':main()
