import unittest
import numpy as np
import pandas as pd
from scipy.optimize import check_grad
from profile_information import model,features,evaluation,splits


class MathematicalContract(unittest.TestCase):
    def test_weighted_interval_gradient_and_missingness(self):
        rng=np.random.default_rng(4);x=rng.normal(size=(5,3));z=rng.normal(size=(4,2));q=rng.normal(size=(5,4,2))
        lo=rng.normal(size=(5,4));hi=lo.copy();lo[0,0]=-np.inf;hi[1,1]=np.inf
        w=np.ones((5,4));w[2,2]=0;lo[2,2]=-np.inf;hi[2,2]=np.inf
        size=1+3+2+3*2+2*(2+1);theta=rng.normal(size=size)
        fun=lambda t:model.objective(t,x,z,q,lo,hi,w,1.3)[0]
        grad=lambda t:model.objective(t,x,z,q,lo,hi,w,1.3)[1]
        self.assertLess(check_grad(fun,grad,theta),2e-5)
        lo2=lo.copy();hi2=hi.copy();lo2[2,2]=999;hi2[2,2]=1000
        a=model.objective(theta,x,z,q,lo,hi,w,1.3);b=model.objective(theta,x,z,q,lo2,hi2,w,1.3)
        self.assertEqual(a[0],b[0]);np.testing.assert_array_equal(a[1],b[1])

    def test_global_baseline_can_predict_nonzero_preference_change(self):
        x=np.array([[-2.],[-1.],[0.],[1.],[2.]])
        z=np.array([[-1.],[1.]])
        y=6+2*x@z.T+x
        q=np.zeros((5,2,0));params,diag=model.fit(x,z,q,y,y,np.ones_like(y),1e-5)
        pred=model.prediction(params,x,z,q)
        self.assertTrue(diag['success']);self.assertLess(np.max(np.abs(pred-y)),1e-4)
        d=(pred[4,0]-pred[0,0])-(pred[4,1]-pred[0,1])
        self.assertLess(d,-15.9)

    def test_row_mean_cancels_native_D(self):
        s=np.array([[3.,6.,9.],[4.,1.,7.]])
        r=s-s.mean(1,keepdims=True)
        fun=lambda a:(a[1,0]-a[0,0])-(a[1,1]-a[0,1])
        self.assertEqual(fun(s),fun(r))

    def test_strict_bound_boundary_order(self):
        self.assertEqual(int(evaluation.direction(0.,np.inf,True,False)),1)
        self.assertEqual(int(evaluation.direction(-np.inf,0.,False,True)),-1)
        self.assertEqual(int(evaluation.direction(0.,0.,False,False)),0)
        self.assertEqual(int(evaluation.direction(-1.,1.,False,False)),0)

    def test_source_decimal_contrast(self):
        self.assertEqual(evaluation.exact_delta_delta(['0.1','0.2','0.3','0.6']),0.)
        self.assertAlmostEqual(evaluation.exact_delta_delta([1,1,0.1,10]),2.)

    def test_standardization_uses_train_inputs_only(self):
        x=np.array([[1.,np.nan],[3.,2.],[100.,200.]])
        a,diag=features.standardize(x,[0,1]);x[2]=[-999,1e20];b,diag2=features.standardize(x,[0,1])
        np.testing.assert_array_equal(a[:2],b[:2]);self.assertEqual(diag,diag2)
        self.assertEqual(a[0,1],0)

    def test_kernel_map_reproduces_training_kernel(self):
        rng=np.random.default_rng(5);x=rng.normal(size=(8,3));k=x@x.T
        f,diag=features.kernel_map(k,k);kc=k-k.mean(0)-k.mean(1)[:,None]+k.mean()
        np.testing.assert_allclose(f@f.T,kc,atol=1e-12)

    def test_unknown_training_genes_are_not_known_target_tests(self):
        assays=pd.DataFrame({'gene':['a','b','c']})
        molecules=pd.DataFrame({'inner_fold_outer0':[0,1,1,0],'inner_fold_outer1':[0,1,1,0]})
        weights=np.array([[1,1,0],[1,1,0],[0,0,1],[0,0,1]],float)
        data=dict(assays=assays,molecules=molecules,weights=weights,eligible_rows=np.arange(4),outer_folds=np.array([0,0,1,1]),chemical_groups=np.array(['a','b','c','d']),row_molecule_indices=np.arange(4),column_genes=np.array(['a','b','c']))
        ss=list(splits.outer_splits(data,['new_molecules_known_proteins']))
        for s in ss:
            known=set(assays.iloc[s['train_assays']].gene)
            self.assertTrue(set(assays.iloc[s['test_assays']].gene)<=known)
            self.assertTrue((weights[np.ix_(s['train_rows'],s['train_assays'])]>0).any(1).all())


if __name__=='__main__':unittest.main()
