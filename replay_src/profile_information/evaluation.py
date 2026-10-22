"""Actual-assay molecular ordering and four-cell preference changes.

No activity/selectivity cutoff. Zero is the algebraic no-change point. Distinct
assays remain distinct. Interval signs are used only when source bounds identify
an order; unknown comparisons never become negative labels.
"""
from itertools import combinations,product
from decimal import Decimal,localcontext
import numpy as np
import pandas as pd


def direction(lower,upper,lower_open=False,upper_open=False):
    lo,hi=np.asarray(lower),np.asarray(upper)
    pos=(lo>0)|((lo==0)&lower_open)
    neg=(hi<0)|((hi==0)&upper_open)
    return pos.astype(np.int8)-neg.astype(np.int8)


def difference(lo_a,hi_a,lo_b,hi_b,lop_a=False,uop_a=False,lop_b=False,uop_b=False):
    """B minus A."""
    return lo_b-hi_a,hi_b-lo_a,np.asarray(lop_b)|np.asarray(uop_a),np.asarray(uop_b)|np.asarray(lop_a)


def exact_delta_delta(concentrations):
    """D=(B,t-A,t)-(B,u-A,u); source cell order A,t;A,u;B,t;B,u."""
    with localcontext() as context:
        context.prec=60
        a,b,c,d=[Decimal(str(v)) for v in concentrations]
        num,den=a*d,b*c
        if num==den:return 0.0
        relative=(num-den)/den
        # Exact decimal product comparison preserves zero/sign; double log1p
        # is stable near one and avoids costly arbitrary-precision logarithms.
        if abs(relative)<Decimal('0.5'):
            return float(np.log1p(float(relative))/np.log(10.0))
        return float(np.log10(float(num/den)))


def quantitative_context(a,b):
    fields=('Format','Readout','Technology','endpoint_interpretation','assay_species')
    values=[]
    for key in fields:
        x,y=a.get(key),b.get(key)
        if pd.isna(x) or pd.isna(y) or not str(x).strip() or not str(y).strip():return False
        values.append(str(x).strip()==str(y).strip())
    return all(values)


def _weighted_variants(data,predictions,molecule_ids,assay_indices,testrows=None):
    """Return per-identity source rows. No predictions averaged before ranking."""
    rows={int(m):np.array([r for r in data['identity_rows'][int(m)] if testrows is None or int(r) in testrows],int) for m in molecule_ids}
    maxvariants=max((len(x) for x in rows.values()),default=0)
    index=np.full((len(molecule_ids),maxvariants),-1,int)
    weights=np.zeros(index.shape)
    for i,m in enumerate(molecule_ids):
        rr=rows[int(m)];index[i,:len(rr)]=rr;weights[i,:len(rr)]=data['variant_weights'][rr]
    return index,weights


def evaluate(data,predictions,test_rows,test_assays,metadata,compact=False,include_records=True):
    """Return record, assay, rectangle, and assay-pair support dataframes.

    `predictions` maps model names to full source-row x actual-assay arrays.
    Learned comparisons share all source variants. Native scores are reported in
    a second group on variants finite for native and all learned comparators.
    """
    testrows=set(map(int,test_rows));mols=set(map(int,data['row_molecule_indices'][test_rows]))
    assays=data['assays'];records=data['records']
    rec=records[records.molecule_index.isin(mols)&records.local_assay_index.isin(test_assays)].copy()
    native_names=[n for n in predictions if n.startswith('native_')]
    learned_names=[n for n in predictions if not n.startswith('native_')]
    groups=[('learned_common',learned_names)] if learned_names else []
    if native_names:groups.append(('native_matched_common',learned_names+native_names))
    record_rows=[];assay_rows=[];rectangles=[];support=[]
    by_assay={int(j):g.set_index('molecule_index') for j,g in rec.groupby('local_assay_index')}
    for j in test_assays:
        j=int(j);g=by_assay.get(j,pd.DataFrame());meta=assays.iloc[j]
        for name,pred in (predictions.items() if include_records else []):
            for row in g.itertuples() if len(g) else []:
                rr=data['identity_rows'][int(row.Index)];rr=np.array([r for r in rr if int(r) in testrows])
                values=pred[rr,j];finite=np.isfinite(values);rr=rr[finite];values=values[finite]
                if len(rr):
                    w=data['variant_weights'][rr];w=w/w.sum()
                    err=values-np.clip(values,row.lower,row.upper)
                    record_rows.append(dict(**metadata,model=name,molecule_index=int(row.Index),assay_index=j,
                        assay_id=meta.assay_id,gene=meta.gene,source_relation=row.source_relation,
                        lower=row.lower,upper=row.upper,lower_open=row.lower_open,upper_open=row.upper_open,
                        predicted_mean=float(w@values),predicted_variance=float(w@((values-w@values)**2)),
                        interval_mse=float(w@(err**2)) if not name.startswith('native_') else np.nan,
                        n_source_variants=len(rr),chemical_group=str(data['molecules'].iloc[int(row.Index)].chemical_group)))
        ids=np.array(sorted(g.index),int) if len(g) else np.empty(0,int)
        if len(ids):g=g.loc[ids]
        left,right=np.triu_indices(len(ids),1)
        if len(left):
            lo,hi,lop,uop=difference(g.lower.to_numpy()[left],g.upper.to_numpy()[left],g.lower.to_numpy()[right],g.upper.to_numpy()[right],g.lower_open.to_numpy()[left],g.upper_open.to_numpy()[left],g.lower_open.to_numpy()[right],g.upper_open.to_numpy()[right])
            truth=direction(lo,hi,lop,uop);is_exact=(g.source_relation.to_numpy()[left]=='=')&(g.source_relation.to_numpy()[right]=='=')
        else:truth=np.empty(0,int);is_exact=np.empty(0,bool)
        ind,vw=_weighted_variants(data,predictions,ids,[j],testrows)
        for comparison,names in groups:
            mass=np.zeros(len(left));numer={n:np.zeros(len(left)) for n in names}
            tie={n:np.zeros(len(left)) for n in names}
            for v1,v2 in product(range(ind.shape[1]),repeat=2):
                ra,rb=ind[left,v1],ind[right,v2];weight=vw[left,v1]*vw[right,v2]
                good=(ra>=0)&(rb>=0)&(truth!=0)
                aidx=np.maximum(ra,0);bidx=np.maximum(rb,0)
                for n in names:good&=np.isfinite(predictions[n][aidx,j])&np.isfinite(predictions[n][bidx,j])
                ww=np.where(good,weight,0);mass+=ww
                for n in names:
                    d=predictions[n][bidx,j]-predictions[n][aidx,j]
                    numer[n]+=np.where(good,ww*(np.sign(d)*truth+1)/2,0)
                    tie[n]+=np.where(good,ww*(d==0),0)
            for domain,mask in [('all_identified',np.ones(len(left),bool)),('exact_exact',is_exact),('exact_censored',~is_exact)]:
                valid=mask&(mass>0)
                for n in names:
                    assay_rows.append(dict(**metadata,comparison_support=comparison,model=n,assay_index=j,
                        assay_id=meta.assay_id,gene=meta.gene,mechanistic_class=meta.mechanistic_class,
                        domain=domain,n_observed_identities=len(ids),n_possible_pairs=len(left),
                        n_identified_pairs=int((mask&(truth!=0)).sum()),n_evaluated_pairs=int(valid.sum()),
                        concordance=float(np.mean(numer[n][valid]/mass[valid])) if valid.any() else np.nan,
                        prediction_tie_fraction=float(np.mean(tie[n][valid]/mass[valid])) if valid.any() else np.nan,
                        undefined_reason='' if valid.any() else 'no interval-identified pair with common predictions'))
    for a,b in combinations(map(int,test_assays),2):
        if assays.iloc[a].gene==assays.iloc[b].gene:continue
        ga=by_assay.get(a);gb=by_assay.get(b)
        ids=np.intersect1d(ga.index,gb.index) if ga is not None and gb is not None else np.empty(0,int)
        la,lb=np.triu_indices(len(ids),1)
        compatible=quantitative_context(assays.iloc[a],assays.iloc[b])
        base=dict(**metadata,assay_index_a=a,assay_index_b=b,assay_id_a=assays.iloc[a].assay_id,
                  assay_id_b=assays.iloc[b].assay_id,gene_a=assays.iloc[a].gene,gene_b=assays.iloc[b].gene,gene_pair_key='|'.join(sorted([str(assays.iloc[a].gene),str(assays.iloc[b].gene)])),
                  quantitative_context_matched=compatible,
                  quantitative_scope='metadata-matched operational IC50 contrast; not thermodynamic Kd selectivity' if compatible else 'all-source operational IC50 contrast; assay metadata not matched; not thermodynamic selectivity')
        sup=dict(**base,n_common_identities=len(ids),n_possible_rectangles=len(la))
        if not len(la):
            support.append(dict(**sup,n_exact_rectangles=0,n_identified_D=0,n_pareto_changes=0,n_saved_rectangles=0,n_exact_nochange_D=0,n_observed_tradeoffs=0,n_observed_any_exact_nochange=0,n_observed_direction_unknown=0));continue
        cache=data.setdefault('_contrast_truth_cache',{})
        cache_key=(a,b,tuple(ids.tolist()))
        if cache_key not in cache:
            ta,tb=ga.loc[ids],gb.loc[ids]
            da=difference(ta.lower.to_numpy()[la],ta.upper.to_numpy()[la],ta.lower.to_numpy()[lb],ta.upper.to_numpy()[lb],ta.lower_open.to_numpy()[la],ta.upper_open.to_numpy()[la],ta.lower_open.to_numpy()[lb],ta.upper_open.to_numpy()[lb])
            db=difference(tb.lower.to_numpy()[la],tb.upper.to_numpy()[la],tb.lower.to_numpy()[lb],tb.upper.to_numpy()[lb],tb.lower_open.to_numpy()[la],tb.upper_open.to_numpy()[la],tb.lower_open.to_numpy()[lb],tb.upper_open.to_numpy()[lb])
            dlo,dhi=da[0]-db[1],da[1]-db[0];dlop,duop=da[2]|db[3],da[3]|db[2]
            exact=(ta.source_relation.to_numpy()[la]=='=')&(ta.source_relation.to_numpy()[lb]=='=')&(tb.source_relation.to_numpy()[la]=='=')&(tb.source_relation.to_numpy()[lb]=='=')
            # Source decimal arithmetic removes fictitious exact ties from separately
            # rounded log concentrations. Strict bounds remain open as supplied.
            ca=ta.ic50_uM.to_numpy();cb=tb.ic50_uM.to_numpy()
            for k in np.flatnonzero(np.isfinite(dlo)|np.isfinite(dhi)):
                val=exact_delta_delta([ca[la[k]],cb[la[k]],ca[lb[k]],cb[lb[k]]])
                if np.isfinite(dlo[k]):dlo[k]=val
                if np.isfinite(dhi[k]):dhi[k]=val
            ds=direction(dlo,dhi,dlop,duop);sa=direction(*da);sb=direction(*db)
            pareto=(sa*sb==-1)
            keep=exact|(ds!=0)|pareto
            tie_a=(da[0]==0)&(da[1]==0)&~da[2]&~da[3];tie_b=(db[0]==0)&(db[1]==0)&~db[2]&~db[3]
            counts=dict(n_exact_rectangles=int(exact.sum()),n_identified_D=int((ds!=0).sum()),n_pareto_changes=int(pareto.sum()),n_saved_rectangles=int(keep.sum()),n_exact_nochange_D=int((exact&(ds==0)).sum()),n_observed_tradeoffs=int((sa*sb==1).sum()),n_observed_any_exact_nochange=int((tie_a|tie_b).sum()),n_observed_direction_unknown=int((((sa==0)&~tie_a)|((sb==0)&~tie_b)).sum()))
            la,lb=la[keep],lb[keep];dlo,dhi,ds,sa,sb,exact,pareto=[x[keep] for x in (dlo,dhi,ds,sa,sb,exact,pareto)]
            cache[cache_key]=(la,lb,dlo,dhi,ds,sa,sb,exact,pareto,counts)
        la,lb,dlo,dhi,ds,sa,sb,exact,pareto,counts=cache[cache_key]
        support.append(dict(**sup,**counts))
        if not len(la):continue
        ind,vw=_weighted_variants(data,predictions,ids,[a,b],testrows);n=len(la)
        for comparison,names in groups:
            mass=np.zeros(n);stats={name:{k:np.zeros(n) for k in ['D','D2','D_mse','D_sign_credit','D_tie','pareto_correct','pareto_reverse','pareto_tradeoff','pareto_tie','delta_a','delta_b']} for name in names}
            for v1,v2 in product(range(ind.shape[1]),repeat=2):
                ra,rb=ind[la,v1],ind[lb,v2];weight=vw[la,v1]*vw[lb,v2]
                good=(ra>=0)&(rb>=0);ia=np.maximum(ra,0);ib=np.maximum(rb,0)
                for name in names:
                    p=predictions[name];good&=np.isfinite(p[ia,a])&np.isfinite(p[ia,b])&np.isfinite(p[ib,a])&np.isfinite(p[ib,b])
                ww=np.where(good,weight,0);mass+=ww
                for name in names:
                    p=predictions[name];pa=p[ib,a]-p[ia,a];pb=p[ib,b]-p[ia,b];d=pa-pb
                    if name.startswith('native_'):
                        scale=np.abs(p[ib,a])+np.abs(p[ia,a])+np.abs(p[ib,b])+np.abs(p[ia,b])
                        ambiguous=good&(np.abs(d)<=8*np.finfo(float).eps*np.maximum(1,scale))
                        for k in np.flatnonzero(ambiguous):
                            d[k]=float(Decimal(str(p[ib[k],a]))-Decimal(str(p[ia[k],a]))-Decimal(str(p[ib[k],b]))+Decimal(str(p[ia[k],b])))
                    add=dict(D=d,D2=d*d,D_mse=(d-np.where(exact,dlo,0))**2,D_sign_credit=(np.sign(d)*ds+1)/2,D_tie=(d==0),
                             pareto_correct=(np.sign(pa)==sa)&(np.sign(pb)==sb),pareto_reverse=(np.sign(pa)==-sa)&(np.sign(pb)==-sb),
                             pareto_tradeoff=(np.sign(pa)*np.sign(pb)==1),pareto_tie=(pa==0)|(pb==0),delta_a=pa,delta_b=pb)
                    for key,val in add.items():stats[name][key]+=np.where(good,ww*val,0)
            valid=np.flatnonzero(mass>0)
            if not len(valid):continue
            ma,mb=ids[la[valid]].astype(int),ids[lb[valid]].astype(int)
            values_by_model={}
            for name in names:
                values={key:arr[valid]/mass[valid] for key,arr in stats[name].items()}
                if name.startswith('native_'):values['D_mse'][:]=np.nan
                else:values['D_mse'][~exact[valid]]=np.nan
                for key in ['D_sign_credit','D_tie']:values[key][ds[valid]==0]=np.nan
                for key in ['pareto_correct','pareto_reverse','pareto_tradeoff','pareto_tie']:values[key][~pareto[valid]]=np.nan
                values_by_model[name]=values
            if compact:
                # Only context sufficient statistics are retained. Every pair
                # still contributes, with no scientific subsampling. Saved
                # held-out predictions support exact detail regeneration.
                for name,values in values_by_model.items():
                    row=dict(**base,comparison_support=comparison,model=name,
                             preference_context_schema='sums_v1',n_rectangles=len(valid),
                             n_molecules=len(np.unique(np.r_[ma,mb])),n_exact=int(exact[valid].sum()),
                             n_identified_D=int((ds[valid]!=0).sum()),n_pareto=int(pareto[valid].sum()))
                    for key in ['D_mse','D_sign_credit','D_tie','pareto_correct','pareto_reverse','pareto_tradeoff','pareto_tie']:
                        v=values[key];finite=np.isfinite(v)
                        row[key+'_sum']=float(v[finite].sum());row[key+'_count']=int(finite.sum())
                    ex=exact[valid];y=dlo[valid][ex];q=values['D'][ex];q2=values['D2'][ex]
                    row.update(D_qy_sum=float(q@y),D_y2_sum=float(y@y),D_q2_sum=float(q2.sum()),D_exact_count=int(ex.sum()))
                    rectangles.append(row)
            else:
                molecular=data['molecules']
                chunk=dict(molecule_a=ma,molecule_b=mb,
                    molecule_id_a=molecular.standard_molecule_id.to_numpy()[ma],molecule_id_b=molecular.standard_molecule_id.to_numpy()[mb],
                    chemical_group_a=molecular.chemical_group.astype(str).to_numpy()[ma],chemical_group_b=molecular.chemical_group.astype(str).to_numpy()[mb],
                    same_chemical_group=molecular.chemical_group.astype(str).to_numpy()[ma]==molecular.chemical_group.astype(str).to_numpy()[mb],
                    observed_D_lower=dlo[valid],observed_D_upper=dhi[valid],observed_D_exact=exact[valid],observed_D_direction=ds[valid],
                    delta_a_direction=sa[valid],delta_b_direction=sb[valid],pareto_identified=pareto[valid],variant_pair_mass=mass[valid])
                for name,values in values_by_model.items():
                    for key,v in values.items():chunk[name+'__'+key]=v
                frame=pd.DataFrame(chunk)
                for key,value in dict(**base,comparison_support=comparison).items():frame[key]=value
                rectangles.append(frame)
    rectframe=pd.DataFrame(rectangles) if compact else (pd.concat(rectangles,ignore_index=True) if rectangles else pd.DataFrame())
    return [pd.DataFrame(record_rows),pd.DataFrame(assay_rows),rectframe,pd.DataFrame(support)]



def summarize_assays(frame):
    if not len(frame):return pd.DataFrame()
    keys=[k for k in ['dataset','channel','scope','policy','regime','variant_set','order','k','comparison_support','model','domain'] if k in frame]
    rows=[]
    for key,g in frame.groupby(keys,dropna=False):
        # Average split-specific contributions to one actual assay first.
        a=g.groupby(['assay_id','gene']).apply(lambda z:pd.Series({'pairs':z.n_evaluated_pairs.sum(),'credit':np.nansum(z.concordance*z.n_evaluated_pairs)}),include_groups=False)
        a['concordance']=a.credit/a.pairs.replace(0,np.nan)
        rows.append(dict(zip(keys,key if isinstance(key,tuple) else (key,)),
                         n_actual_assays=len(a),n_defined_assays=int(a.concordance.notna().sum()),n_pairs=int(a.pairs.sum()),
                         macro_assay_concordance=float(a.concordance.mean()),macro_gene_concordance=float(a.reset_index().groupby('gene').concordance.mean().mean()),
                         pooled_concordance=float(a.credit.sum()/a.pairs.sum()) if a.pairs.sum() else np.nan))
    return pd.DataFrame(rows)


def summarize_preferences(frame):
    """Continuous effect summaries with explicit alternative observation weights."""
    if not len(frame):return pd.DataFrame()
    names=sorted({c[:-len('__D')] for c in frame if c.endswith('__D')})
    meta_keys=[k for k in ['dataset','channel','scope','policy','regime','variant_set','order','k','comparison_support'] if k in frame]
    output=[]
    for meta,group in frame.groupby(meta_keys,dropna=False):
        metadata=dict(zip(meta_keys,meta if isinstance(meta,tuple) else (meta,)))
        for domain,g in [('all_source_operational',group),('metadata_matched_operational',group[group.quantitative_context_matched])]:
            if not len(g):continue
            for weighting in ['contrast_equal','assay_pair_equal_adjusted_alternatives','gene_pair_equal','chemical_pair_equal']:
                # Assay alternatives share the mass of one chemical/gene contrast.
                contrast_keys=['molecule_a','molecule_b','gene_pair_key']
                w=np.ones(len(g)) if weighting in ['raw_rectangle_equal','actual_assay_pair_equal'] else 1/g.groupby(contrast_keys,dropna=False).molecule_a.transform('size').to_numpy(float)
                key={'actual_assay_pair_equal':['assay_id_a','assay_id_b'],'assay_pair_equal_adjusted_alternatives':['assay_id_a','assay_id_b'],'gene_pair_equal':['gene_pair_key'],'chemical_pair_equal':['molecule_a','molecule_b']}.get(weighting)
                if key:
                    temp=g[key].copy();temp['_w']=w;w/=temp.groupby(key,dropna=False)._w.transform('sum').to_numpy()
                for name in names:
                    def avg(values):
                        values=np.asarray(values,float);valid=np.isfinite(values)
                        return float(w[valid]@values[valid]/w[valid].sum()) if valid.any() else np.nan
                    exact=g.observed_D_exact.to_numpy(bool)
                    y=np.where(exact,g.observed_D_lower.to_numpy(),np.nan)
                    q=g[name+'__D'].to_numpy(float);q2=g[name+'__D2'].to_numpy(float)
                    qy=avg(q*y);yy=avg(y*y);qq=avg(np.where(exact,q2,np.nan))
                    mse=avg(g[name+'__D_mse'])
                    output.append(dict(**metadata,domain=domain,weighting=weighting,model=name,
                        n_rectangles=len(g),n_molecules=len(set(g.molecule_a)|set(g.molecule_b)),
                        n_assay_pairs=g[['assay_id_a','assay_id_b']].drop_duplicates().shape[0],
                        n_gene_pairs=g.gene_pair_key.nunique(),
                        n_exact=int(exact.sum()),n_identified_D=int((g.observed_D_direction!=0).sum()),n_pareto=int(g.pareto_identified.sum()),
                        D_rmse=float(np.sqrt(mse)) if np.isfinite(mse) else np.nan,
                        D_direction_concordance=avg(g[name+'__D_sign_credit']),
                        D_predicted_tie_fraction=avg(g[name+'__D_tie']),
                        symmetrized_D_correlation=qy/np.sqrt(yy*qq) if yy>0 and qq>0 else np.nan,
                        pareto_correct_fraction=avg(g[name+'__pareto_correct']),pareto_reversed_fraction=avg(g[name+'__pareto_reverse']),
                        pareto_tradeoff_fraction=avg(g[name+'__pareto_tradeoff']),pareto_tie_fraction=avg(g[name+'__pareto_tie'])))
    primary=summarize_contexts(detail_contexts(frame))
    return pd.concat([primary,pd.DataFrame(output)],ignore_index=True)


def summarize_contexts(frame):
    """Aggregate streamed context sums: equal assay-pair and pooled contrasts."""
    if not len(frame):return pd.DataFrame()
    meta_keys=[k for k in ['dataset','channel','scope','policy','regime','variant_set','order','k','comparison_support','model'] if k in frame]
    out=[];metrics=['D_mse','D_sign_credit','D_tie','pareto_correct','pareto_reverse','pareto_tradeoff','pareto_tie']
    for key,group in frame.groupby(meta_keys,dropna=False):
        metadata=dict(zip(meta_keys,key if isinstance(key,tuple) else (key,)))
        for domain,g in [('all_source_operational',group),('metadata_matched_operational',group[group.quantitative_context_matched])]:
            if not len(g):continue
            columns=[m+suffix for m in metrics for suffix in ['_sum','_count']]+['D_qy_sum','D_y2_sum','D_q2_sum','D_exact_count','n_rectangles','n_exact','n_identified_D','n_pareto']
            contexts=g.groupby(['assay_id_a','assay_id_b','gene_a','gene_b'])[columns].sum()
            for weighting in ['actual_assay_pair_equal','raw_rectangle_equal']:
                row=dict(**metadata,domain=domain,weighting=weighting,n_assay_pairs=len(contexts),
                    n_gene_pairs=g.gene_pair_key.nunique(),
                    n_rectangles=int(contexts.n_rectangles.sum()),n_exact=int(contexts.n_exact.sum()),
                    n_identified_D=int(contexts.n_identified_D.sum()),n_pareto=int(contexts.n_pareto.sum()))
                for m in metrics:
                    count=contexts[m+'_count'];total=contexts[m+'_sum']
                    row[m]=float((total/count.replace(0,np.nan)).mean()) if weighting=='actual_assay_pair_equal' else (float(total.sum()/count.sum()) if count.sum() else np.nan)
                if weighting=='actual_assay_pair_equal':
                    denom=np.sqrt(contexts.D_y2_sum*contexts.D_q2_sum);row['symmetrized_D_correlation']=float((contexts.D_qy_sum/denom.replace(0,np.nan)).mean())
                else:
                    denom=np.sqrt(contexts.D_y2_sum.sum()*contexts.D_q2_sum.sum());row['symmetrized_D_correlation']=float(contexts.D_qy_sum.sum()/denom) if denom>0 else np.nan
                row['D_rmse']=float(np.sqrt(row['D_mse'])) if np.isfinite(row['D_mse']) else np.nan
                row.update(D_direction_concordance=row['D_sign_credit'],D_predicted_tie_fraction=row['D_tie'],pareto_correct_fraction=row['pareto_correct'],pareto_reversed_fraction=row['pareto_reverse'],pareto_tradeoff_fraction=row['pareto_tradeoff'],pareto_tie_fraction=row['pareto_tie'])
                for original in ['D_sign_credit','D_tie','pareto_correct','pareto_reverse','pareto_tradeoff','pareto_tie']:row.pop(original)
                out.append(row)
    return pd.DataFrame(out)


def detail_contexts(frame):
    """Convert regenerated detail to the exact streamed sufficient-statistic schema."""
    if not len(frame):return pd.DataFrame()
    names=sorted({c[:-len('__D')] for c in frame if c.endswith('__D')})
    group_keys=[k for k in ['dataset','channel','scope','policy','regime','split','variant_set','order','k','comparison_support','assay_index_a','assay_index_b','assay_id_a','assay_id_b','gene_a','gene_b','gene_pair_key','quantitative_context_matched','quantitative_scope'] if k in frame]
    out=[]
    for key,g in frame.groupby(group_keys,dropna=False):
        metadata=dict(zip(group_keys,key if isinstance(key,tuple) else (key,)))
        for name in names:
            valid=g[name+'__D'].notna();sub=g[valid]
            if not len(sub):continue
            exact=sub.observed_D_exact.to_numpy(bool);y=sub.observed_D_lower.to_numpy()[exact]
            q=sub[name+'__D'].to_numpy()[exact];q2=sub[name+'__D2'].to_numpy()[exact]
            row=dict(**metadata,model=name,preference_context_schema='sums_v1',n_rectangles=len(sub),n_molecules=len(set(sub.molecule_a)|set(sub.molecule_b)),
                     n_exact=int(exact.sum()),n_identified_D=int((sub.observed_D_direction!=0).sum()),n_pareto=int(sub.pareto_identified.sum()),
                     D_qy_sum=float(q@y),D_y2_sum=float(y@y),D_q2_sum=float(q2.sum()),D_exact_count=int(exact.sum()))
            for metric in ['D_mse','D_sign_credit','D_tie','pareto_correct','pareto_reverse','pareto_tradeoff','pareto_tie']:
                values=sub[name+'__'+metric];row[metric+'_sum']=float(values.sum());row[metric+'_count']=int(values.notna().sum())
            out.append(row)
    return pd.DataFrame(out)
