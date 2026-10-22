"""Relocated, storage-bounded public smoke; no scientific fits/full-study replay."""
from pathlib import Path
import hashlib,json,os,re,shutil,subprocess,sys,tempfile,time
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
LIMIT=500*1024*1024


def size(folder):return sum(p.stat().st_size for p in folder.rglob('*') if p.is_file())


def main():
    started=time.monotonic();qa=ROOT/'.qa';qa.mkdir(exist_ok=True)
    checks=[];peak=0;comparisons=[]
    with tempfile.TemporaryDirectory(prefix='public-relocated-',dir=qa) as temporary:
        fresh=Path(temporary)/'repo';fresh.mkdir()
        for name in ['src','growth_src','decoder_src','replay_src','data','geometry_inputs','reference','postprocessing','growth_analysis','growth_analysis_support','tests','growth_tests','decoder_tests','scripts','tools','configs','provenance']:
            shutil.copytree(ROOT/name,fresh/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc','*.log'))
        for name in ['pyproject.toml','requirements.txt','README.md','LICENSE','LICENSE-DATA.md','CITATION.cff']:shutil.copyfile(ROOT/name,fresh/name)
        (fresh/'release_assets').mkdir();shutil.copyfile(ROOT/'release_assets/fingerprint_run01.tar.gz',fresh/'release_assets/fingerprint_run01.tar.gz')
        environment=os.environ.copy();environment['TMPDIR']=str(fresh/'.tmp');Path(environment['TMPDIR']).mkdir();environment.pop('PYTHONPATH',None)
        def run(label,args):
            nonlocal peak
            now=time.monotonic();result=subprocess.run([sys.executable,*args],cwd=fresh,env=environment,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            (qa/(label+'.log')).write_text(result.stdout)
            peak=max(peak,size(fresh))
            if peak>LIMIT:raise RuntimeError('Temporary storage limit exceeded')
            checks.append(dict(check=label,command=['python',*args],exit_code=result.returncode,seconds=time.monotonic()-now,temporary_bytes_after=size(fresh)))
            if result.returncode:raise RuntimeError('Public smoke failed at '+label+'; see ignored QA log')
            print(json.dumps(dict(check=label,exit_code=0,temporary_bytes=peak)),flush=True)
        run('source_hashes',['scripts/verify_bundle.py'])
        run('unit_checks',['tools/run_checks.py'])
        run('geometry',['scripts/replay_geometry.py'])
        run('restore_fingerprint',['tools/package_assets.py','restore','release_assets/fingerprint_run01.tar.gz'])
        run('public_manifest',['tools/build_replay_manifests.py','--stage','postprocessing','--config','configs/postprocessing_smoke.json'])
        run('public_postprocessing',['scripts/replay_results.py','postprocessing','--config','configs/postprocessing_smoke.json'])
        run('public_measurement_strata',['scripts/replay_results.py','measurement-strata','--config','configs/measurement_strata_smoke.json'])
        dataset='fingerprint_run01/own44__own44__all_binding_genes__all_records'
        for actual_root,reference_root,names in [
            (fresh/'outputs/public_smoke_postprocessing',ROOT/'results/profile_postprocessing_run02',{'paired_estimates.csv','paired_uncertainty.csv','paired_association.csv','paired_association_uncertainty.csv'}),
            (fresh/'outputs/public_smoke_measurement_strata',ROOT/'results/profile_measurement_strata_run02',{'paired_estimates.csv','paired_uncertainty.csv','stratum_support.csv','stratum_status.csv'}),
        ]:
            files=[p for p in actual_root.rglob('*.csv') if p.name in names and dataset in str(p.relative_to(actual_root))]
            for actual in files:
                relative=actual.relative_to(actual_root);reference=reference_root/relative
                a=pd.read_csv(actual);b=pd.read_csv(reference)
                keys=[k for k in ['measurement_stratum','domain','weighting','model','metric','resampling'] if k in a]
                if not keys:keys=['measurement_stratum']
                a=a.sort_values(keys).reset_index(drop=True);b=b.sort_values(keys).reset_index(drop=True)
                if len(a)!=len(b):raise AssertionError('Public row membership differs: '+str(relative))
                pd.testing.assert_frame_equal(a[keys],b[keys],check_dtype=False)
                numeric=[c for c in a.select_dtypes(include=np.number) if c in b]
                np.testing.assert_allclose(a[numeric],b[numeric],rtol=2e-10,atol=2e-10,equal_nan=True)
                diff=np.abs(a[numeric].to_numpy(float)-b[numeric].to_numpy(float))
                comparisons.append(dict(table=str(relative),rows=len(a),numeric_columns=len(numeric),numeric_cells=int(a[numeric].notna().sum().sum()),max_abs_difference=float(np.nanmax(diff)) if diff.size and np.isfinite(diff).any() else 0.0,matched=True))
        peak=max(peak,size(fresh))
        if peak>LIMIT:raise RuntimeError('Temporary storage limit exceeded')
    receipt=dict(status='passed',scope='Fresh relocated code/data copy with installed pinned analysis dependencies; original private workspace is not accessed by the public commands. Only one representative public saved-prediction dataset is replayed.',checks=checks,peak_temporary_bytes=peak,temporary_limit_bytes=LIMIT,temporary_copy_removed=True,representative_dataset='fingerprint_run01/own44__own44__all_binding_genes__all_records',saved_splits=10,regimes=['new_molecules_known_proteins','new_proteins_known_molecules'],resamples=400,compact_comparisons=comparisons,numeric_tolerance=dict(rtol=2e-10,atol=2e-10,scientific_threshold=False),scientific_study_fits=False,synthetic_test_fits=True,full_public_replay=False,network_download_tested=False,archive_transport='Existing verified fingerprint release archive copied into the fresh public layout; archive checksum and restored scientific hashes checked.',seconds=time.monotonic()-started)
    (ROOT/'provenance/final_bounded_public_qa.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(status='passed',tables=len(comparisons),peak_bytes=peak,temporary_copy_removed=True)),flush=True)


if __name__=='__main__':main()
