"""Numerical equivalence to the slow direct rectangle resampling oracle."""
import unittest
import numpy as np
import pandas as pd
from test_postprocess import fixture
from postprocess import summary as slow
from fast_resampling import summary as fast


class SparseEquivalence(unittest.TestCase):
    def compare(self,f):
        a,ad,ai=slow(f,resamples=23,seed=18)
        b,bd,bi=fast(f,resamples=23,seed=18,batch_size=7)
        pd.testing.assert_frame_equal(a,b,atol=1e-12,rtol=1e-12)
        keys=['domain','weighting','model','metric','resampling','replicate']
        pd.testing.assert_frame_equal(ad.sort_values(keys).reset_index(drop=True),bd.sort_values(keys).reset_index(drop=True),atol=1e-12,rtol=1e-12)
        keys=keys[:-1]
        pd.testing.assert_frame_equal(ai.sort_values(keys).reset_index(drop=True),bi.sort_values(keys).reset_index(drop=True),atol=1e-12,rtol=1e-12)

    def test_same_draws_same_estimands(self):self.compare(fixture())

    def test_same_group_dyads_missing_arm_and_source_moments(self):
        f=fixture()
        f.loc[0,'chemical_group_b']='a'
        f.loc[0,'full_profile__D2']=3
        f.loc[0,'full_profile__D_mse']=2
        f.loc[2,'full_profile__D']=np.nan
        self.compare(f)

    def test_all_undefined_model_is_omitted(self):
        f=fixture()
        for c in f:
            if c.startswith('native_score__'):f[c]=np.nan
        a,_,_=fast(f,resamples=2)
        self.assertFalse(a.model.eq('native_score').any())


if __name__=='__main__':unittest.main()
