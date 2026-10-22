"""Independent tiny cases for weights, dependency and paired uncertainty."""
import unittest
import numpy as np
import pandas as pd
from postprocess import node_factor,prepare_context,metric_estimate,summary,descriptive_slope


def fixture():
    f=pd.DataFrame(dict(assay_id_a=[1,1,2],assay_id_b=[3,3,3],gene_a=['X','X','Y'],gene_b=['Z','Z','Z'],chemical_group_a=['a','a','b'],chemical_group_b=['b','c','c'],molecule_a=[0,0,1],molecule_b=[1,2,2],observed_D_exact=[True,True,False],observed_D_lower=[1.,1.,np.nan],observed_D_direction=[1,1,0],pareto_identified=[True,True,False],quantitative_context_matched=[True,True,True]))
    for model in ['global_chemistry','full_profile','native_score']:
        for metric,v in dict(D=[1.,1.,2.],D2=[1.,1.,4.],D_mse=[0.,0.,np.nan],D_sign_credit=[1.,1.,0.],D_tie=[0.,0.,np.nan],pareto_correct=[1.,0.,np.nan],pareto_reverse=[0.,1.,np.nan],pareto_tradeoff=[0.,0.,np.nan],pareto_tie=[0.,0.,np.nan]).items():
            f[model+'__'+metric]=v
    return f


class WeightTests(unittest.TestCase):
    def test_actual_assay_pair_and_rectangle_are_distinct(self):
        f=prepare_context(fixture())
        self.assertEqual(metric_estimate(f,'full_profile','D_direction_concordance'),.5)
        self.assertAlmostEqual(metric_estimate(f,'full_profile','D_direction_concordance',weighting='raw_rectangle_equal'),2/3)

    def test_unknown_is_not_in_denominator(self):
        f=prepare_context(fixture());f.loc[2,'full_profile__D_sign_credit']=np.nan
        self.assertEqual(metric_estimate(f,'full_profile','D_direction_concordance'),1.)

    def test_dyadic_factors_do_not_treat_four_cells_as_independent(self):
        got=node_factor(np.array([2,3,0]),np.array([0,0,1,2]),np.array([1,0,2,2]))
        np.testing.assert_array_equal(got,[6,2,0,0])

    def test_resampling_stays_paired(self):
        points,draws,intervals=summary(fixture(),resamples=9,seed=12)
        finite=draws.arm_minus_baseline.dropna()
        self.assertTrue(np.all(finite==0))
        self.assertTrue(np.all(intervals.paired_delta_q025.dropna()==0))
        self.assertFalse(((points.model=='native_score') & points.metric.isin(['D_mse','D_rmse'])).any())
        _,again,_=summary(fixture(),resamples=9,seed=12)
        pd.testing.assert_frame_equal(draws,again)

    def test_exact_variant_second_moment_not_square_of_mean(self):
        f=prepare_context(fixture().iloc[:1].copy())
        f.loc[0,'full_profile__D']=1
        f.loc[0,'full_profile__D2']=2
        self.assertAlmostEqual(metric_estimate(f,'full_profile','symmetrized_D_correlation'),1/np.sqrt(2))

    def test_context_resampling_preserves_assay_pair_weights(self):
        f=prepare_context(fixture())
        # Context X-Z counted six times; Y-Z absent. It remains a mean within X-Z.
        got=metric_estimate(f,'full_profile','D_direction_concordance',context_mass=np.array([6,0]))
        self.assertEqual(got,1)

    def test_missing_arm_pairs_baseline_on_same_rows(self):
        f=fixture();f.loc[:1,"full_profile__D"]=np.nan
        points,_,_=summary(f,resamples=0)
        row=points[(points.model=="full_profile") & (points.metric=="D_direction_concordance") & (points.domain=="all_source_operational") & (points.weighting=="actual_assay_pair_equal")].iloc[0]
        self.assertEqual(row.n_paired_rectangles,1)
        self.assertEqual(row.estimate,0)
        self.assertEqual(row.baseline_estimate,0)

    def test_continuous_slope(self):
        slope,intercept,n=descriptive_slope(np.array([0.,1.,2.]),np.array([2.,5.,8.]),np.ones(3))
        self.assertEqual((slope,intercept,n),(3,2,3))


if __name__=='__main__':unittest.main()
