"""Train-only chemical/protein features and explicitly missing score encodings."""
import numpy as np


def kernel_map(k_train,k_all):
    k=np.asarray(k_train,float);cross=np.asarray(k_all,float)
    mean=k.mean(0);grand=float(mean.mean())
    centered=k-mean[None,:]-mean[:,None]+grand
    ev,u=np.linalg.eigh((centered+centered.T)/2)
    tolerance=np.finfo(float).eps*max(k.shape)*max(1.0,float(np.max(np.abs(ev))))
    if ev.min() < -100*tolerance:raise ValueError('Kernel not numerically PSD')
    keep=ev>tolerance
    basis=u[:,keep]/np.sqrt(ev[keep])[None,:]
    features=(cross-mean[None,:]-cross.mean(1)[:,None]+grand)@basis
    return features,dict(rank=int(keep.sum()),numerical_tolerance=float(tolerance),train_kernel_rows=len(k))


def standardize(values,train,weights=None):
    """Fit each feature's moments on observed training inputs only.

    Missing standardized feature entries are coded as zero with availability
    provided separately. This does not change or impute any raw score/outcome.
    """
    x=np.asarray(values,float)
    if x.ndim==1:x=x[:,None]
    train=np.asarray(train,int)
    w=np.ones(len(train)) if weights is None else np.asarray(weights,float)
    t=x[train];ok=np.isfinite(t)
    mass=(w[:,None]*ok).sum(0)
    mean=np.divide((w[:,None]*np.where(ok,t,0)).sum(0),mass,out=np.zeros(x.shape[1]),where=mass>0)
    var=np.divide((w[:,None]*np.where(ok,(t-mean)**2,0)).sum(0),mass,out=np.zeros(x.shape[1]),where=mass>0)
    sd=np.sqrt(var);sd=np.where(sd>0,sd,1.0)
    result=np.where(np.isfinite(x),(x-mean)/sd,0.0)
    return result,dict(mean=mean.tolist(),scale=sd.tolist(),training_observed_mass=mass.tolist(),missing_encoding='zero standardized input plus explicit observed mask')


def rbf(a,b):
    distance=np.maximum((a*a).sum(1)[:,None]+(b*b).sum(1)[None,:]-2*a@b.T,0)
    return np.exp(-distance/(2*max(1,a.shape[1])))


def chemistry(data,train_rows,descriptor_weight=0.5):
    train_mol=np.unique(data['row_molecule_indices'][train_rows])
    all_mol=data['row_molecule_indices']
    desc=data['descriptors']
    zdesc,diag=standardize(desc,train_mol)
    tani=data['tanimoto']
    kt=(1-descriptor_weight)*tani[np.ix_(train_mol,train_mol)]+descriptor_weight*rbf(zdesc[train_mol],zdesc[train_mol])
    ka=(1-descriptor_weight)*tani[np.ix_(all_mol,train_mol)]+descriptor_weight*rbf(zdesc[all_mol],zdesc[train_mol])
    features,kdiag=kernel_map(kt,ka)
    return features,dict(training_molecule_indices=train_mol.tolist(),descriptor=diag,kernel=kdiag,descriptor_weight=float(descriptor_weight),definition='(1-w) MorganR3 Tanimoto + w standardized descriptor RBF')


def proteins(data,train_assays,known):
    n=len(data['assays'])
    spectrum=data['protein_features'][data['assays'].protein_index.to_numpy(int)]
    k=spectrum@spectrum.T
    # Known assay IDs permit separate responses/slopes; protein similarity adds
    # a transferable prior. Neither feature contains experimental outcomes.
    if known:k=(k+np.eye(n))/2
    return kernel_map(k[np.ix_(train_assays,train_assays)],k[:,train_assays])


def score_features(data,train_rows,selected_columns):
    cols=np.asarray(selected_columns,int);s=data['scores'][:,cols]
    observed=np.isfinite(s);count=observed.sum(1)
    g=np.divide(np.where(observed,s,0).sum(1),count,out=np.full(len(s),np.nan),where=count>0)
    residual=s-g[:,None]
    weights=data['variant_weights'][train_rows]
    zg,gdiag=standardize(g,train_rows,weights)
    av,adiag=standardize(observed.astype(float),train_rows,weights)
    rz,rdiag=standardize(residual,train_rows,weights)
    coverage=count/max(1,len(cols))
    zc,cdiag=standardize(coverage,train_rows,weights)
    # Scores of unpurchased columns are absent from every input here.
    base=np.c_[zg,zc,av/np.sqrt(max(1,len(cols)))]
    full=rz/np.sqrt(max(1,len(cols)))
    genes=data['column_genes'][cols]
    lookup={str(gene):i for i,gene in enumerate(genes)}
    if len(lookup)!=len(genes):raise ValueError('Multiple contexts of one gene need separate source datasets')
    n_assays=len(data['assays']);native=np.full((len(s),n_assays),np.nan)
    target_residual=np.full_like(native,np.nan)
    for j,gene in enumerate(data['assays'].gene.astype(str)):
        if gene in lookup:
            native[:,j]=s[:,lookup[gene]]
            target_residual[:,j]=residual[:,lookup[gene]]
    # One global score-unit transform estimated only on training ligand/assay
    # entries will be fitted in assemble() once the actual training assays are known.
    return base,full,native,target_residual,dict(columns=cols.tolist(),global_mean=gdiag,availability=adiag,profile=rdiag,coverage=cdiag,definition='native observed-column mean and source-row residual; no raw score filling')


def assemble(data,train_rows,train_assays,selected_columns,known,variant,cache=None,descriptor_weight=0.5,column_selection=None):
    if cache is None:
        xc,cdiag=chemistry(data,train_rows,descriptor_weight);z,zdiag=proteins(data,train_assays,known)
    else:xc,cdiag,z,zdiag=cache
    base,profile,native,residual,sdiag=score_features(data,train_rows,selected_columns)
    x=np.c_[xc,base]
    if variant=='full_profile':x=np.c_[x,profile]
    elif variant not in ('global_chemistry','target_residual','global_plus_selected_column','chemistry_plus_selected_score','chemistry_plus_target_score'):raise ValueError(variant)
    if variant in ('chemistry_plus_target_score','chemistry_plus_selected_score'):x=xc
    if variant in ('global_plus_selected_column','chemistry_plus_selected_score'):
        if column_selection is None:raise ValueError('Training-only selected columns required')
        replacement=np.full_like(native,np.nan)
        for j,col in enumerate(column_selection):
            if col>=0:replacement[:,j]=data['scores'][:,col]
        # Compare the selected column to exactly the same purchased-panel mean.
        observed=np.isfinite(data['scores'][:,selected_columns]);count=observed.sum(1)
        g=np.divide(np.where(observed,data['scores'][:,selected_columns],0).sum(1),count,out=np.full(len(x),np.nan),where=count>0)
        native=replacement;residual=replacement-g[:,None]
    if variant in ('chemistry_plus_target_score','chemistry_plus_selected_score'):residual=native.copy()
    available=np.isfinite(native).astype(float)
    # Both arms have exactly the same matching-target score availability signal.
    q=available[:,:,None]
    if variant in ('target_residual','global_plus_selected_column','chemistry_plus_selected_score','chemistry_plus_target_score'):
        tr=residual[np.ix_(train_rows,train_assays)]
        w=data['variant_weights'][train_rows,None]*np.isfinite(tr)
        mass=w.sum();mean=float((w*np.where(np.isfinite(tr),tr,0)).sum()/mass) if mass else 0.0
        var=float((w*np.where(np.isfinite(tr),(tr-mean)**2,0)).sum()/mass) if mass else 0.0
        scale=np.sqrt(var) if var>0 else 1.0
        q=np.dstack([available,np.where(np.isfinite(residual),(residual-mean)/scale,0)])
        sdiag['matched_residual']=dict(mean=mean,scale=float(scale),training_observed_mass=float(mass))
    return x,z,q,native,dict(chemistry=cdiag,protein=zdiag,score=sdiag,variant=variant)


def nearest_columns(data,selected_columns):
    cols=np.asarray(selected_columns,int)
    ids=data['column_protein_indices'][cols]
    valid=ids>=0
    result=np.full(len(data['assays']),-1,int)
    if valid.any():
        pf=data['protein_features'];sim=pf[data['assays'].protein_index.to_numpy(int)]@pf[ids[valid]].T
        result=cols[valid][np.argmax(sim,axis=1)]
    return result


def select_columns(data,train_rows,train_assays,selected_columns,known):
    """Univariate column selection uses training assay labels only.

    Equal identity mass is preserved across source variants. Censored rows enter
    only interval-identified within-assay orders. Ranking includes both signs;
    the subsequent regularized regression estimates the sign. A label-unidentified
    assay falls back to its sequence-nearest candidate, recorded explicitly.
    """
    selected=np.asarray(selected_columns,int);choice=nearest_columns(data,selected)
    report=[]
    if not known:
        return choice,[dict(assay_index=int(j),column_index=int(choice[j]),rule='sequence_nearest_no_labels') for j in range(len(choice))]
    for j in train_assays:
        rows=np.asarray(train_rows)[data['weights'][train_rows,j]>0]
        left,right=np.triu_indices(len(rows),1)
        a,b=rows[left],rows[right]
        different=data['row_molecule_indices'][a]!=data['row_molecule_indices'][b]
        lo,hi=data['lower'],data['upper'];lop,uop=data['lower_open'],data['upper_open']
        pos=(lo[a,j]>hi[b,j])|((lo[a,j]==hi[b,j])&(lop[a,j]|uop[b,j]))
        neg=(hi[a,j]<lo[b,j])|((hi[a,j]==lo[b,j])&(uop[a,j]|lop[b,j]))
        truth=pos.astype(int)-neg.astype(int);use=different&(truth!=0)
        a,b,truth=a[use],b[use],truth[use]
        if len(a) and len(selected):
            x=data['scores'][np.ix_(a,selected)];y=data['scores'][np.ix_(b,selected)]
            ok=np.isfinite(x)&np.isfinite(y);w=(data['variant_weights'][a]*data['variant_weights'][b])[:,None]
            mass=(w*ok).sum(0)
            association=np.divide((w*np.where(ok,np.sign(x-y)*truth[:,None],0)).sum(0),mass,out=np.full(len(selected),np.nan),where=mass>0)
            if np.isfinite(association).any():
                k=int(np.nanargmax(np.abs(association)));choice[j]=selected[k]
                report.append(dict(assay_index=int(j),column_index=int(choice[j]),rule='training_interval_order_absolute_association',association=float(association[k]),identity_pair_mass=float(mass[k]),candidate_columns=len(selected)))
                continue
        report.append(dict(assay_index=int(j),column_index=int(choice[j]),rule='no_identified_training_order_sequence_fallback',candidate_columns=len(selected)))
    return choice,report
