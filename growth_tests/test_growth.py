import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from profile_information import features,run


class GrowthContract(unittest.TestCase):
    def test_zero_panel_does_not_read_score_arrays(self):
        class ForbiddenScores:
            def __getitem__(self,key):
                raise AssertionError('Zero-cost anchor read scores')
        data=dict(scores=ForbiddenScores(),row_molecule_indices=np.arange(3),assays=pd.DataFrame({'gene':['a','b']}))
        xc=np.arange(6).reshape(3,2);z=np.eye(2);cache=(xc,{},z,{})
        with patch.object(features,'score_features',side_effect=AssertionError('Score extraction called')):
            for arm in ['global_chemistry','full_profile']:
                x,p,q,native,diag=features.assemble(data,[0,1],[0,1],[],True,arm,cache=cache)
                np.testing.assert_array_equal(x,xc);np.testing.assert_array_equal(p,z)
                self.assertEqual(q.shape,(3,2,0));self.assertTrue(np.isnan(native).all())
                self.assertFalse(diag['score']['score_inputs_used'])

    def test_budgets_include_zero_and_preserve_default(self):
        self.assertEqual(list(run.growth_budgets({},3)),[1,2,3])
        self.assertEqual(list(run.growth_budgets({'start_k':0},3)),[0,1,2,3])
        for invalid in [-1,4,0.5,True]:
            with self.assertRaises(ValueError):run.growth_budgets({'start_k':invalid},3)

    def test_resume_checks_signature_and_artifact_bytes(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);(folder/'result.txt').write_text('original')
            config=dict(_analysis_signature='v1',resume=True)
            run.mark_complete(folder,dict(status='complete'),config)
            self.assertTrue(run.completed_matches(folder,config))
            with self.assertRaisesRegex(RuntimeError,'signature mismatch'):
                run.completed_matches(folder,dict(_analysis_signature='v2'))
            (folder/'result.txt').write_text('changed')
            with self.assertRaisesRegex(RuntimeError,'missing/changed'):run.completed_matches(folder,config)

    def test_analysis_signature_changes_for_start_k_and_input(self):
        with tempfile.TemporaryDirectory() as temp:
            inputs=Path(temp);(inputs/'source.txt').write_text('a')
            a=run.analysis_signature(dict(growth=dict(start_k=0)),inputs)['signature']
            self.assertNotEqual(a,run.analysis_signature(dict(growth=dict(start_k=1)),inputs)['signature'])
            (inputs/'source.txt').write_text('b')
            self.assertNotEqual(a,run.analysis_signature(dict(growth=dict(start_k=0)),inputs)['signature'])




class CanonicalGrowthContract(unittest.TestCase):
    def fixture(self):
        rng=np.random.default_rng(79);n=12
        data=dict(scores=rng.normal(size=(n,4)),variant_weights=np.ones(n),
                  row_molecule_indices=np.arange(n),assays=pd.DataFrame({'gene':['a','c']}),
                  column_genes=np.array(['a','b','c','d']))
        data['scores'][2,1]=np.nan;data['scores'][10,3]=np.nan
        xc=rng.normal(size=(n,3));z=rng.normal(size=(2,2));cache=(xc,{},z,{})
        return data,np.arange(8),np.arange(2),cache

    def fit(self,assembled):
        from threadpoolctl import threadpool_limits
        from profile_information import model
        x,z,q=assembled[:3];rng=np.random.default_rng(15)
        y=6+rng.normal(size=(8,len(z)));w=np.ones_like(y)
        with threadpool_limits(limits=1):
            theta,diag=model.fit(x[:8],z,q[:8],y,y,w,1.,dict(maxiter=5000,ftol=1e-10,gtol=1e-6,maxls=40))
            prediction=model.prediction(theta,x,z,q)
        return theta,prediction

    def test_same_purchased_set_has_bitwise_identical_features_fit_prediction(self):
        data,rr,aa,cache=self.fixture()
        for variant in ['global_chemistry','full_profile']:
            a=features.assemble(data,rr,aa,[3,0,2,1],True,variant,cache=cache)
            b=features.assemble(data,rr,aa,[1,2,3,0],True,variant,cache=cache)
            for x,y in zip(a[:4],b[:4]):
                self.assertEqual(x.shape,y.shape);self.assertEqual(x.tobytes(),y.tobytes())
            self.assertEqual(a[4]['score']['canonical_columns'],[0,1,2,3])
            ta,pa=self.fit(a);tb,pb=self.fit(b)
            self.assertEqual(ta.tobytes(),tb.tobytes());self.assertEqual(pa.tobytes(),pb.tobytes())

    def test_k1_global_and_full_are_bitwise_identical_with_missing_scores(self):
        data,rr,aa,cache=self.fixture()
        a=features.assemble(data,rr,aa,[1],True,'global_chemistry',cache=cache)
        b=features.assemble(data,rr,aa,[1],True,'full_profile',cache=cache)
        for x,y in zip(a[:4],b[:4]):
            self.assertEqual(x.shape,y.shape);self.assertEqual(x.tobytes(),y.tobytes())
        self.assertEqual(b[4]['score']['kept_profile_columns'],[])
        ta,pa=self.fit(a);tb,pb=self.fit(b)
        self.assertEqual(ta.tobytes(),tb.tobytes());self.assertEqual(pa.tobytes(),pb.tobytes())

    def test_only_exact_zero_training_profile_columns_are_dropped(self):
        data,rr,aa,cache=self.fixture()
        baseline=features.assemble(data,rr,aa,[0,1],True,'global_chemistry',cache=cache)
        # One held-out-only feature is unlearnable; the tiny training input stays.
        profile=np.zeros((12,2));profile[10,0]=1;profile[0,1]=1e-250
        original=features.score_features
        def scored(*args):
            base,_,native,residual,diag=original(*args)
            return base,profile,native,residual,diag
        with patch.object(features,'score_features',side_effect=scored):
            full=features.assemble(data,rr,aa,[1,0],True,'full_profile',cache=cache)
        self.assertEqual(full[4]['score']['kept_profile_columns'],[1])
        self.assertEqual(full[0].shape[1],baseline[0].shape[1]+1)
        np.testing.assert_array_equal(full[0][:,-1],profile[:,1])

if __name__=='__main__':unittest.main()
