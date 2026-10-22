"""Synthetic-only technical checks; no saved scientific labels are fitted."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from profile_information import evaluation, features, model, splits
from profile_information import decoder_sensitivity as ds
from profile_information.run import write_json, mark_complete


def fixture():
    rng=np.random.default_rng(21);n=9;genes=np.array(['A','B','C'])
    molecules=pd.DataFrame({'standard_molecule_id':['m'+str(i) for i in range(n)],'chemical_group':['g'+str(i) for i in range(n)],'inner_fold_outer0':[0,0,1,1,1,1,1,1,1]})
    assay=pd.DataFrame({'gene':genes,'assay_id':['a','b','c'],'protein_index':np.arange(3),'mechanistic_class':['x']*3,'Format':['f']*3,'Readout':['r']*3,'Technology':['t']*3,'endpoint_interpretation':['b']*3,'assay_species':['human']*3})
    row_mol=np.r_[np.arange(n),0];nr=len(row_mol)
    scores=rng.normal(size=(nr,3));scores[4,1]=np.nan;scores[-1]=scores[0]+.2
    vw=np.ones(nr);vw[[0,nr-1]]=.5
    concentrations=np.array([[1+i,2+2*i+(i%2),3+(i%3)] for i in range(n)],float)
    y=6-np.log10(concentrations)
    lower=y[row_mol].copy();upper=lower.copy();lop=np.zeros_like(lower,bool);uop=lop.copy()
    lower[5,1]=-np.inf;uop[5,1]=True
    records=[]
    for i in range(n):
        for j in range(3):
            censored=(i==5 and j==1)
            records.append(dict(molecule_index=i,local_assay_index=j,lower=lower[i,j],upper=upper[i,j],lower_open=lop[i,j],upper_open=uop[i,j],source_relation='>' if censored else '=',ic50_uM=concentrations[i,j]))
    desc=rng.normal(size=(n,4));fp=rng.normal(size=(n,6));kt=features.rbf(fp,fp)
    return dict(dataset='synthetic',channel='synthetic',scope='native_genes',policy='all_records',molecules=molecules,assays=assay,records=pd.DataFrame(records),row_molecule_indices=row_mol,identity_rows={i:np.flatnonzero(row_mol==i) for i in range(n)},variant_weights=vw,column_genes=genes,column_protein_indices=np.arange(3),scores=scores,lower=lower,upper=upper,lower_open=lop,upper_open=uop,weights=np.broadcast_to(vw[:,None],(nr,3)).copy(),descriptors=desc,tanimoto=kt,protein_features=np.eye(3),chemical_groups=molecules.chemical_group.to_numpy()[row_mol],outer_folds=np.r_[np.ones(6,int),np.zeros(3,int),1],eligible_rows=np.arange(nr))


class DecoderSensitivityContract(unittest.TestCase):
    def test_weighted_map_reproduces_centered_psd_kernel(self):
        rng=np.random.default_rng(8);x=rng.normal(size=(6,4));k=features.rbf(x,x);w=np.array([.5,.5,1,1,1,1])
        f,diag=ds.weighted_kernel_map(k,k,w);v=w/w.sum();mu=k@v;kc=k-mu[:,None]-mu[None,:]+v@mu
        np.testing.assert_allclose(f@f.T,kc,atol=1e-12)
        np.testing.assert_allclose(v@f,0,atol=1e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(f@f.T).min(),-1e-12)

    def test_train_only_and_column_relabeling(self):
        data=fixture();rr=np.array([0,1,2,3,4,9]);cols=np.arange(3)
        a,da=ds.profile_blocks(data,rr,cols)
        data['scores'][[5,6,7,8]]*=10000
        b,db=ds.profile_blocks(data,rr,cols)
        self.assertEqual(da,db)
        for decoder in a:np.testing.assert_allclose(a[decoder][rr],b[decoder][rr],atol=1e-12)
        c,dc=ds.profile_blocks(data,rr,cols[::-1])
        np.testing.assert_allclose(b['linear']@b['linear'].T,c['linear']@c['linear'].T,atol=1e-12)
        np.testing.assert_allclose(b['rbf']@b['rbf'].T,c['rbf']@c['rbf'].T,atol=1e-10)
        self.assertEqual(dc['unique_source_rows'],len(rr))
        # Numerical construction check: unscaled standardized inputs get the
        # dimension divisor once, and this fixture has nontrivial distances.
        raw=data['scores'];g=np.nanmean(raw,axis=1);residual=raw-g[:,None]
        rz,_=features.standardize(residual,rr,data['variant_weights'][rr])
        av,_=features.standardize(np.isfinite(raw).astype(float),rr,data['variant_weights'][rr])
        expected=features.rbf(np.c_[rz,av][rr],np.c_[rz,av][rr])
        offdiag=expected[np.triu_indices(len(rr),1)]
        np.testing.assert_allclose(db['kernel']['training_offdiagonal_kernel']['quantiles'],np.quantile(offdiag,[0,.1,.25,.5,.75,.9,1]),atol=1e-12)
        self.assertLess(np.median(offdiag),.9)  # Synthetic numerical fixture only.
        with self.assertRaises(ValueError):ds.profile_blocks(data,[0,0,1],cols)

    def test_strength_zero_exact_baseline(self):
        data=fixture();split=dict(train_rows=np.array([0,1,2,3,4,5,9]),train_assays=np.arange(3),known=True)
        base=ds.assemble_base(data,split,np.arange(3),.5);algorithm={'optimizer':{'maxiter':5000}}
        a,pa,da=ds.fit_candidate(data,split,base,np.full((10,2),1e100),0.,1.,algorithm)
        x,z,q,_=base;rr=split['train_rows'];aa=split['train_assays']
        pb,db=model.fit(x[rr],z[aa],q[rr][:,aa],data['lower'][np.ix_(rr,aa)],data['upper'][np.ix_(rr,aa)],data['weights'][np.ix_(rr,aa)],1.,algorithm['optimizer'])
        np.testing.assert_array_equal(pa,pb);np.testing.assert_array_equal(a,model.prediction(pb,x,z,q))

    def test_objective_uses_exact_outer_direction_and_equal_assay_pairs(self):
        data=fixture();rows=np.arange(10);aa=np.arange(3);pred=np.nan_to_num(data['lower'],neginf=5)
        split=dict(test_rows=rows,test_assays=aa)
        trials,frames=ds.candidate_objectives(data,{'test':pred,'reverse':-pred},split)
        contexts=frames[2];c=contexts[contexts.model=='test'];v=c.D_sign_credit_sum/c.D_sign_credit_count
        self.assertAlmostEqual(trials[0]['validation_D_direction_concordance'],v.mean())
        self.assertAlmostEqual(trials[0]['validation_D_direction_concordance']+trials[1]['validation_D_direction_concordance'],1.)
        self.assertEqual(int(evaluation.direction(0,np.inf,True,False)),1)
        self.assertEqual(evaluation.exact_delta_delta([.1,.2,.3,.6]),0.)

    def test_ties_choose_baseline_and_undefined_is_not_relabelled(self):
        trials=[dict(candidate='baseline0' if s==0 else d,decoder=d,strength=s,validation_interval_mse=1.,validation_D_direction_concordance=np.nan) for d in ['linear','rbf'] for s in [0,.1,1,10]]
        choices=ds.choose_candidates(trials)
        self.assertEqual(choices['linear__pointwise']['strength'],0)
        self.assertIsNone(choices['linear__direction'])

    def test_synthetic_complete_split(self):
        data=fixture();outer=dict(id='mol_fold_0',regime='new_molecules_known_proteins',outer_fold=0,train_rows=np.array([0,1,2,3,4,5,9]),test_rows=np.array([6,7,8]),train_assays=np.arange(3),test_assays=np.arange(3))
        inner=splits.inner_split(data,outer);inner['known']=True
        config=dict(_analysis_signature='synthetic',resume=True,algorithm={'inner_validation_fold':0,'optimizer':{'maxiter':5000}},sensitivity={'strength_grid':[0.,.1,1.,10.]})
        with tempfile.TemporaryDirectory() as temp:
            frozen=Path(temp)/'frozen';out=Path(temp)/'out';frozen.mkdir()
            write_json(frozen/'split.json',outer);write_json(frozen/'tuning.json',dict(chosen_descriptor_weight=.5,chosen_alpha={'global_chemistry':1.},inner_split=inner))
            base=ds.assemble_base(data,outer,np.arange(3),.5)
            pred,par,diag=ds.fit_candidate(data,outer,base,None,0.,1.,config['algorithm'])
            np.savez_compressed(frozen/'predictions.npz',**{n:pred for n in ['global_chemistry','global_plus_selected_column','full_profile']})
            for n in ['global_chemistry','global_plus_selected_column','full_profile']:np.savez_compressed(frozen/(n+'_model.npz'),coefficients=par)
            mark_complete(frozen,{'status':'synthetic'},config)
            ds.run_split(data,outer,config,frozen,out)
            result=np.load(out/'predictions.npz');self.assertEqual(len(result.files),9)
            tuning=json.loads((out/'inner_tuning.json').read_text());self.assertEqual(len(tuning['trials']),8)
            for arm,choice in tuning['choices'].items():
                if choice['strength']==0:np.testing.assert_array_equal(result[arm],pred)
            self.assertTrue((out/'inner/candidate_coefficients.npz').exists())
            self.assertTrue((out/'complete.json').exists())
            def nonzero_choices(trials):
                return {d+'__'+o:next(t for t in trials if t['decoder']==d and t['strength']==.1) for d in ['linear','rbf'] for o in ['pointwise','direction']}
            with patch.object(ds,'choose_candidates',side_effect=nonzero_choices), patch.object(ds,'fit_candidate',wraps=ds.fit_candidate) as fits:
                ds.run_split(data,outer,config,frozen,Path(temp)/'nonzero')
                self.assertEqual(fits.call_count,9)  # Seven inner, two shared outer fits.
            nonzero=np.load(Path(temp)/'nonzero/predictions.npz')
            np.testing.assert_array_equal(nonzero['linear__pointwise'],nonzero['linear__direction'])
            np.testing.assert_array_equal(nonzero['rbf__pointwise'],nonzero['rbf__direction'])
            self.assertTrue((Path(temp)/'nonzero/linear_s0.1_model.npz').exists())


if __name__=='__main__':unittest.main()
