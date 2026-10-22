"""Group-disjoint molecules and simultaneously held-out protein contexts."""
from itertools import combinations
import numpy as np


def _observed_train(data,rows,assays):
    rows=np.asarray(rows,int);assays=np.asarray(assays,int)
    if not len(rows) or not len(assays):return rows[:0],assays[:0]
    weight=data['weights'][np.ix_(rows,assays)]
    return rows[(weight>0).any(1)],assays[(weight>0).any(0)]


def _known_rows(data,rows,train_assays):
    rows=np.asarray(rows,int)
    if not len(train_assays):return rows[:0]
    return rows[(data['weights'][np.ix_(rows,train_assays)]>0).any(1)]


def target_blocks(data,scheme='pairs',n_folds=5):
    present=set(data['assays'].gene.astype(str))
    genes=[g for g in data['column_genes'] if g in present]
    genes+=sorted(present-set(genes))
    if scheme=='pairs':return [('pair_'+str(i),list(pair)) for i,pair in enumerate(combinations(genes,2))]
    if scheme=='gene_folds':return [('genes_fold_'+str(f),[g for i,g in enumerate(genes) if i%n_folds==f]) for f in range(n_folds)]
    raise ValueError(scheme)


def outer_splits(data,regimes,scheme='pairs',n_target_folds=5):
    rows=data['eligible_rows'];all_assays=np.arange(len(data['assays']));genes=data['assays'].gene.astype(str).to_numpy()
    if 'new_molecules_known_proteins' in regimes:
        for f in sorted(np.unique(data['outer_folds'][rows])):
            candidate_train=rows[data['outer_folds'][rows]!=f];test=rows[data['outer_folds'][rows]==f]
            train,tr_assays=_observed_train(data,candidate_train,all_assays)
            if set(data['chemical_groups'][train])&set(data['chemical_groups'][test]):raise ValueError('Chemical group leakage')
            known_genes=set(genes[tr_assays]);test_assays=all_assays[np.isin(genes,list(known_genes))]
            unknown_assays=all_assays[~np.isin(genes,list(known_genes))]
            yield dict(id='mol_fold_'+str(f),regime='new_molecules_known_proteins',outer_fold=int(f),
                       train_rows=train,test_rows=test,train_assays=tr_assays,test_assays=test_assays,
                       heldout_genes=[],target_scheme='known',
                       excluded_test_assays_no_training_gene=unknown_assays,
                       excluded_test_record_weight_no_training_gene=float(data['weights'][np.ix_(test,unknown_assays)].sum()))
    if 'new_proteins_known_molecules' in regimes or 'both_cold' in regimes:
        for name,heldout in target_blocks(data,scheme,n_target_folds):
            tr=all_assays[~np.isin(genes,heldout)];te=all_assays[np.isin(genes,heldout)]
            if 'new_proteins_known_molecules' in regimes:
                known=_known_rows(data,rows,tr);train,train_assays=_observed_train(data,known,tr)
                excluded=np.setdiff1d(rows,known)
                yield dict(id=name,regime='new_proteins_known_molecules',outer_fold=None,
                           train_rows=train,test_rows=known,train_assays=train_assays,test_assays=te,
                           heldout_genes=heldout,target_scheme=scheme,
                           excluded_test_rows_no_training_protein_label=excluded,
                           excluded_test_record_weight_no_training_protein_label=float(data['weights'][np.ix_(excluded,te)].sum()))
            if 'both_cold' in regimes:
                for f in sorted(np.unique(data['outer_folds'][rows])):
                    candidate_train=rows[data['outer_folds'][rows]!=f];test=rows[data['outer_folds'][rows]==f]
                    train,train_assays=_observed_train(data,candidate_train,tr)
                    yield dict(id=name+'_mol_fold_'+str(f),regime='both_cold',outer_fold=int(f),
                               train_rows=train,test_rows=test,train_assays=train_assays,test_assays=te,
                               heldout_genes=heldout,target_scheme=scheme)


def inner_split(data,outer,validation_fold=0):
    rows=outer['train_rows'];assays=outer['train_assays'];g=data['assays'].gene.astype(str).to_numpy()
    if outer['regime'] in ('new_molecules_known_proteins','both_cold'):
        f=outer['outer_fold'];field='inner_fold_outer'+str(f)
        labels=data['molecules'][field].to_numpy(int)[data['row_molecule_indices'][rows]]
        tr=rows[labels!=validation_fold];va=rows[labels==validation_fold]
        if set(data['chemical_groups'][tr])&set(data['chemical_groups'][va]):raise ValueError('Inner chemical group leakage')
        if outer['regime']=='both_cold':
            ug=sorted(set(g[assays]));val=set(ug[validation_fold%3::3])
            ta=assays[~np.isin(g[assays],list(val))];vs=assays[np.isin(g[assays],list(val))]
        else:ta=vs=assays
        tr,ta=_observed_train(data,tr,ta)
        if outer['regime']=='new_molecules_known_proteins':vs=vs[np.isin(g[vs],list(set(g[ta])))]
        return dict(train_rows=tr,test_rows=va,train_assays=ta,test_assays=vs)
    ug=sorted(set(g[assays]));val=set(ug[validation_fold%3::3])
    ta=assays[~np.isin(g[assays],list(val))];vs=assays[np.isin(g[assays],list(val))]
    known=_known_rows(data,rows,ta);tr,ta=_observed_train(data,known,ta)
    return dict(train_rows=tr,test_rows=known,train_assays=ta,test_assays=vs)
