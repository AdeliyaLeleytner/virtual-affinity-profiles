"""Synthetic-only technical checks; never executes the scientific replay."""
import unittest
import numpy as np
import pandas as pd
from growth_analysis.analyze import same_bits, validate_curve_summary, selected_costs


class TechnicalTests(unittest.TestCase):
    def test_bitwise_nan_and_signed_zero(self):
        self.assertTrue(same_bits(np.array([1., np.nan]), np.array([1., np.nan])))
        self.assertFalse(same_bits(np.array([0.]), np.array([-0.])))

    def test_curve_grid_and_support(self):
        rows=[]
        for model in ['global_chemistry','full_profile']:
            for order in ['random_1','family_stratified_1']:
                for k in range(3):
                    rows.append(dict(model=model,order=order,k=k,domain='all_source_operational',
                        weighting='actual_assay_pair_equal',comparison_support='learned_common',
                        n_assay_pairs=2,n_gene_pairs=2,n_rectangles=8,n_exact=4,n_identified_D=8,n_pareto=3,
                        D_direction_concordance=.5+k*.1,pareto_correct_fraction=.2+k*.1))
        frame=pd.DataFrame(rows)
        self.assertTrue(validate_curve_summary(frame,2,['random_1','family_stratified_1'])['support_constant'])
        frame.loc[0,'n_rectangles']=7
        with self.assertRaises(ValueError):validate_curve_summary(frame,2,['random_1','family_stratified_1'])

    def test_query_union_is_not_one_shared_column(self):
        data=dict(column_genes=np.array(['A','B','C']), assays=pd.DataFrame({'gene':['A','B'],'assay_id':['a','b']}),
                  weights=np.ones((2,2)))
        split=dict(id='fold0',train_rows=[0,1],train_assays=[0,1],test_assays=[0,1])
        diag={'chemistry_plus_selected_score':{'column_selection':[{'assay_index':0,'column_index':2},{'assay_index':1,'column_index':2}]}}
        support=pd.DataFrame([dict(assay_index_a=0,assay_index_b=1,assay_id_a='a',assay_id_b='b',gene_a='A',gene_b='B',n_saved_rectangles=1)])
        assays,pairs,unions=selected_costs(data,split,diag,{},support)
        pairs=pd.DataFrame(pairs).set_index('model')
        self.assertEqual(pairs.loc['chemistry_plus_target_score','deployment_query_column_count'],2)
        self.assertEqual(pairs.loc['chemistry_plus_selected_score','deployment_query_column_count'],1)
        self.assertEqual(pairs.loc['chemistry_plus_selected_score','training_score_candidate_columns'],3)
        self.assertTrue(pairs.loc['chemistry_plus_selected_score','training_search_not_free'])


if __name__=='__main__':unittest.main()
