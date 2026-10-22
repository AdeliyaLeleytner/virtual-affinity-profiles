"""Small direct fixtures for measurement-class support and resumable isolation."""
from pathlib import Path
import json
import tempfile
import unittest
import numpy as np
import pandas as pd
from test_postprocess import fixture
from measurement_strata import STRATA, METRICS, strata, support_counts, analyse, run_from_config, sha


def measurement_fixture():
    f=fixture()
    # Two exact rows (one exact no-change), one censor-involving identified row,
    # and one censor-involving unknown row. Unknown is not a negative label.
    f.loc[1,'observed_D_direction']=0
    f.loc[1,'pareto_identified']=False
    f.loc[2,'observed_D_direction']=1
    f.loc[2,'quantitative_context_matched']=False
    extra=f.iloc[[2]].copy();extra['observed_D_direction']=0
    f=pd.concat([f,extra],ignore_index=True)
    for model in ['global_chemistry','full_profile','native_score']:
        f.loc[[1,3],model+'__D_sign_credit']=np.nan
        f.loc[2,model+'__D_sign_credit']=1
        for key in ['pareto_correct','pareto_reverse','pareto_tradeoff','pareto_tie']:
            f.loc[1,model+'__'+key]=np.nan
    f['_row_id']=np.arange(len(f))
    return f


class MeasurementStrata(unittest.TestCase):
    def test_disjoint_strata_recover_all_identified_direction_support(self):
        f=measurement_fixture();parts=strata(f)
        self.assertEqual(set(parts[STRATA[0]]._row_id),{0,1})
        self.assertEqual(set(parts[STRATA[1]]._row_id),{2})
        combined=pd.concat(parts.values())
        self.assertEqual(set(combined.loc[combined.observed_D_direction.ne(0),'_row_id']),set(f.loc[f.observed_D_direction.ne(0),'_row_id']))
        self.assertNotIn(3,set(combined._row_id))

    def test_counts_keep_exact_ties_and_empty_metadata_domain(self):
        c=support_counts(measurement_fixture()).set_index(['measurement_stratum','domain'])
        exact=c.loc[STRATA[0],'all_source_operational']
        self.assertEqual((exact.observed_rectangle_rows,exact.observed_identified_D_rows,exact.observed_exact_nochange_D_rows),(2,1,1))
        censored=c.loc[STRATA[1],'all_source_operational']
        self.assertEqual((censored.observed_rectangle_rows,censored.observed_exact_rows,censored.observed_identified_D_rows),(1,0,1))
        self.assertEqual(c.loc[STRATA[1],'metadata_matched_operational'].observed_rectangle_rows,0)

    def test_only_direction_choice_outputs_and_paired_draws(self):
        point,draws,intervals,counts,status=analyse(measurement_fixture(),resamples=5,seed=13)
        for output in [point,draws,intervals]:
            self.assertTrue(set(output.metric).issubset(METRICS))
            self.assertEqual(set(output.measurement_stratum),set(STRATA))
        self.assertTrue(np.all(draws.arm_minus_baseline.dropna()==0))
        self.assertEqual(set(draws.resampling),{'chemical_group_fixed_panel','gene_context_fixed_chemistry','joint_chemical_gene_sensitivity'})
        self.assertEqual(set(point.weighting),{'actual_assay_pair_equal','raw_rectangle_equal'})

    def test_all_unknown_is_explicit_empty_support(self):
        point,draws,intervals,counts,status=analyse(measurement_fixture().iloc[[3]],resamples=5)
        self.assertTrue(point.empty and draws.empty and intervals.empty)
        self.assertEqual(counts.observed_rectangle_rows.sum(),0)
        self.assertTrue(status.status.eq("empty_observation_pool").all())

    def test_own_config_resume_and_input_immutability(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            root=Path(tmp)
            relative=Path('core_run02/toy__toy__native_genes__all_records/new_molecules_known_proteins/native_matched_common/continuous_rectangle_data.csv.gz')
            source=root/'input'/relative;source.parent.mkdir(parents=True)
            measurement_fixture().to_csv(source,index=False)
            config=root/'config.json';config.write_text(json.dumps(dict(input_dir='input',output_dir='output',resamples=2,seed=3)))
            before=sha(source)
            first=run_from_config(config,project_root=root)
            second=run_from_config(config,project_root=root)
            self.assertEqual((first['completed'],second['resumed']),(1,1))
            self.assertEqual(before,sha(source))
            estimates=pd.read_csv(root/'output'/relative.parent/'paired_estimates.csv')
            self.assertTrue(estimates.comparison_support.eq('native_matched_common').all())
            self.assertTrue(estimates.regime.eq('new_molecules_known_proteins').all())
            (root/'output'/relative.parent/'paired_estimates.csv').write_text('corrupted\n')
            repaired=run_from_config(config,project_root=root)
            self.assertEqual(repaired['completed'],1)


if __name__=='__main__':unittest.main()
