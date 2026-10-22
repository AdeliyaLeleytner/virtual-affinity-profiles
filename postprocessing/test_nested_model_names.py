"""Pure-prediction regression: nested decoder names must remain distinct."""
import unittest
import numpy as np
import pandas as pd
from profile_information import evaluation


def fixture():
    row_molecule=np.array([0,0,1,2]);weights=np.array([.25,.75,1.,1.])
    concentration=np.array([[1.,2.],[4.,1.],[2.,8.]])
    y=6-np.log10(concentration)
    assays=pd.DataFrame(dict(assay_id=['t','u'],gene=['T','U'],mechanistic_class=['x','x'],Format=['f','f'],Readout=['r','r'],Technology=['tech','tech'],endpoint_interpretation=['IC50','IC50'],assay_species=['human','human']))
    records=[]
    for i in range(3):
        for j in range(2):records.append(dict(molecule_index=i,local_assay_index=j,lower=y[i,j],upper=y[i,j],lower_open=False,upper_open=False,source_relation='=',ic50_uM=concentration[i,j]))
    data=dict(assays=assays,records=pd.DataFrame(records),molecules=pd.DataFrame(dict(standard_molecule_id=['A','B','C'],chemical_group=['a','b','c'])),row_molecule_indices=row_molecule,identity_rows={i:np.flatnonzero(row_molecule==i) for i in range(3)},variant_weights=weights)
    base=y[row_molecule]
    predictions={
        'global_chemistry':base.copy(),
        'full_profile':base+np.array([[0,.2],[.1,0],[.2,-.1],[0,.3]]),
        'linear__pointwise':base*.8,
        'linear__direction':base[:,::-1],
        'rbf__pointwise':base*1.2,
        'rbf__direction':-base,
        'native_score':np.array([[-8.,-3.],[-7.,-4.],[-6.,-9.],[-4.,-5.]]),
    }
    metadata=dict(dataset='fixture',channel='fixture',scope='native_genes',policy='all_records',regime='new_molecules_known_proteins',split='fixture',variant_set='core')
    return data,predictions,metadata


def sort(frame):
    keys=[k for k in ['comparison_support','model','domain','weighting','assay_id_a','assay_id_b'] if k in frame]
    return frame.sort_values(keys).reset_index(drop=True).sort_index(axis=1)


class NestedModelNames(unittest.TestCase):
    def test_nested_and_legacy_contexts_equal_direct_compact_evaluation(self):
        data,pred,meta=fixture();rows=np.arange(4);assays=np.arange(2)
        detail=evaluation.evaluate(data,pred,rows,assays,meta,compact=False,include_records=False)[2]
        compact=evaluation.evaluate(data,pred,rows,assays,meta,compact=True,include_records=False)[2]
        regenerated=evaluation.detail_contexts(detail)
        self.assertEqual(set(regenerated.model),set(pred))
        pd.testing.assert_frame_equal(sort(regenerated),sort(compact),check_dtype=False,rtol=1e-12,atol=1e-12)

    def test_full_preference_summary_preserves_values_after_name_relabelling(self):
        data,pred,meta=fixture();rows=np.arange(4);assays=np.arange(2)
        nested=evaluation.evaluate(data,pred,rows,assays,meta,compact=False,include_records=False)[2]
        rename={name:name.replace('__','_') for name in pred}
        legacy=evaluation.evaluate(data,{rename[name]:value for name,value in pred.items()},rows,assays,meta,compact=False,include_records=False)[2]
        actual=evaluation.summarize_preferences(nested)
        expected=evaluation.summarize_preferences(legacy)
        self.assertEqual(set(actual.model),set(pred))
        actual['model']=actual.model.map(rename)
        pd.testing.assert_frame_equal(sort(actual),sort(expected),check_dtype=False,rtol=1e-12,atol=1e-12)


if __name__=='__main__':unittest.main()
