"""Prepared, provenance-retaining inputs. No assay pooling or score averaging."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd


def load(input_dir,dataset,channel=None,scope='native_genes',policy='all_records'):
    root=Path(input_dir)
    molecules=pd.read_csv(root/'molecules.csv',float_precision='round_trip')
    assays=pd.read_csv(root/'assays.csv',float_precision='round_trip')
    records=pd.read_csv(root/'binding_records.csv.gz',float_precision='round_trip')
    proteins=pd.read_csv(root/'proteins.csv')
    mol=np.load(root/'molecule_features.npz',allow_pickle=False)
    protein=np.load(root/'protein_features.npz',allow_pickle=False)
    source=np.load(root/'datasets'/(dataset+'.npz'),allow_pickle=False)
    scores=source['scores'].copy()
    if scores.ndim==3:
        if channel is None:raise ValueError('SPD channel must be explicit')
        index=list(source['channel_names'].astype(str)).index(channel)
        scores=scores[index]
        genes=source['target_ids'].astype(str)
        column_ids=genes.copy()
    else:
        if channel is not None:raise ValueError('Unexpected channel')
        genes=source['column_genes'].astype(str);column_ids=source['column_ids'].astype(str)
    rows=source['row_molecule_indices'].astype(int)
    variant_weights=source['variant_weights'].astype(float)
    if len(set(genes))!=len(genes):raise ValueError('Keep receptor contexts in separate source datasets')
    if scores.shape!=(len(rows),len(genes)):raise ValueError('Scores shape')
    if policy=='representative':
        flag=records.representative_result_drug_assay_group_pair.astype(str).str.lower().eq('true')
        records=records[flag].copy()
    elif policy!='all_records':raise ValueError(policy)
    if records.duplicated(['molecule_index','assay_id']).any():raise ValueError('Actual-assay identity duplicate')
    records=records[records.molecule_index.isin(rows)].copy()
    if scope=='native_genes':records=records[records.gene.isin(genes)].copy()
    elif scope!='all_binding_genes':raise ValueError(scope)
    # Keep all actual assays in the designated gene universe, including zero-label
    # or all-censored contexts. Undefined evaluations are reported downstream.
    assays=assays[assays.gene.isin(genes) if scope=='native_genes' else np.ones(len(assays),bool)].copy()
    assays=assays.sort_values('assay_index').reset_index(drop=True)
    amap=dict(zip(assays.assay_id,range(len(assays))))
    records['local_assay_index']=records.assay_id.map(amap).astype(int)
    identity_rows={int(m):np.flatnonzero(rows==m) for m in np.unique(rows)}
    lower=np.full((len(rows),len(assays)),np.nan);upper=lower.copy();w=np.zeros_like(lower)
    lower_open=np.zeros(lower.shape,bool);upper_open=lower_open.copy()
    for rec in records.itertuples(index=False):
        ix=identity_rows[int(rec.molecule_index)];j=int(rec.local_assay_index)
        lower[ix,j]=rec.lower;upper[ix,j]=rec.upper;w[ix,j]=variant_weights[ix]
        lower_open[ix,j]=rec.lower_open;upper_open[ix,j]=rec.upper_open
    for ix in identity_rows.values():
        if not np.isclose(variant_weights[ix].sum(),1):raise ValueError('Variant weights do not total one')
    eligible=(w>0).any(1)&(molecules.outer_fold.to_numpy()[rows]>=0)
    pmap=dict(zip(proteins.gene.astype(str),proteins.protein_index.astype(int)))
    column_protein_indices=np.array([pmap.get(g,-1) for g in genes],int)
    return dict(input_dir=root,dataset=dataset,channel=channel or dataset,scope=scope,policy=policy,
                molecules=molecules,assays=assays,records=records,proteins=proteins,
                row_molecule_indices=rows,identity_rows=identity_rows,variant_weights=variant_weights,
                column_genes=genes,column_ids=column_ids,column_protein_indices=column_protein_indices,
                scores=scores,lower=lower,upper=upper,lower_open=lower_open,upper_open=upper_open,
                weights=w,eligible_rows=np.flatnonzero(eligible),
                descriptors=mol['descriptors'],tanimoto=mol['ligand_tanimoto'],
                protein_features=protein['protein_spectrum'],
                outer_folds=molecules.outer_fold.to_numpy(int)[rows],
                chemical_groups=molecules.chemical_group.astype(str).to_numpy()[rows])


def describe(data):
    use=data['weights']>0
    return dict(dataset=data['dataset'],channel=data['channel'],scope=data['scope'],policy=data['policy'],
                source_variant_rows=len(data['scores']),input_identities=len(data['identity_rows']),
                evaluated_identities=len(np.unique(data['row_molecule_indices'][data['eligible_rows']])),
                score_columns=data['scores'].shape[1],assays=len(data['assays']),genes=data['assays'].gene.nunique(),
                source_records=len(data['records']),exact_records=int(data['records'].source_relation.eq('=').sum()),
                censored_records=int(data['records'].source_relation.ne('=').sum()),
                missing_raw_scores=int(np.isnan(data['scores']).sum()),
                observed_variant_cells=int(use.sum()),identity_weighted_observations=float(data['weights'].sum()),
                numeric_interval='source IC50 uM -> pIC50; no endpoint conversion to Kd',
                native_context='one source-specific receptor context per nominal gene; canonical sequence does not certify construct')
