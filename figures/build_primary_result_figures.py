"""Display registered paired-result tables; no fits, resampling or selection.
python build_primary_result_figures.py --figure broad|spd --results-root PATH
Writes vector PDF/SVG, PNG previews, full source snapshots and precise plotted data.
Normal public cached rendering uses --from-cache; --results-root PATH can
refresh tables from a newly completed replay without fitting or resampling.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from orx_figstyle import use_style,save,PALETTE,WIDE

R=Path(__file__).resolve().parent
REGIMES=[('new_molecules_known_proteins','New molecules / known proteins'),
         ('new_proteins_known_molecules','New proteins / known molecules')]
METRICS=[('D_direction_concordance','Δ D-direction concordance (pp)'),
         ('pareto_correct_fraction','Δ Joint opposite-direction choice (pp)')]
ARMS={'full_profile':('Full profile',PALETTE['blue'],'o'),
      'target_residual':('Matching-target residual',PALETTE['orange'],'s'),
      'global_plus_selected_column':('Selected-column residual',PALETTE['purple'],'D'),
      'native_score':('Native score (unfitted)',PALETTE['green'],'^')}
SPD={'boltz_affinity':'Boltz-2 affinity','boltzina_vina_far_affinity':'Boltzina far',
     'boltzina_vina_zero_affinity':'Boltzina zero','gnina':'GNINA on Vina pose',
     'nesso_affinity':'Nesso-1 affinity','raw_ad4':'AutoDock4','raw_gnina_own':'GNINA own search','vina':'Vina'}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot(path,destination,receipt):
    before=sha(path);destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,destination)
    after=sha(path)
    if before!=after or sha(destination)!=after:raise ValueError(f'Source changed during snapshot: {path}')
    receipt.append(dict(source=str(path.resolve()),cached=str(destination.relative_to(R)),sha256=after))
    return pd.read_csv(destination,float_precision='round_trip')


def snapshot_json(path,destination,receipt):
    before=sha(path);destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,destination)
    after=sha(path)
    if before!=after or sha(destination)!=after:raise ValueError(f'Source changed during snapshot: {path}')
    receipt.append(dict(source=str(path.resolve()),cached=str(destination.relative_to(R)),sha256=after))
    return json.loads(destination.read_text())


def collect(args,specs,arms):
    receipts=[];rows=[];native=[];native_decisions=[];sensitivity=[]
    cache=R/'data'/('figure3_sources' if args.figure=='broad' else 'figure4_sources')
    for key,label in specs.items():
        name=(f'{key}__{key}__native_genes__all_records' if args.figure=='broad' else f'spd__{key}__native_genes__all_records')
        node=args.results_root/'core_run02'/name
        state=json.loads((node/'dataset_state.json').read_text())
        if not state.get('dataset_complete'):raise ValueError(f'Dataset incomplete: {name}')
        description_path=(args.prediction_root/name/'input_description.json' if args.prediction_root is not None else Path(state['source'])/'input_description.json')
        description=snapshot_json(description_path,cache/name/'input_description.json',receipts)
        absolute=snapshot(node/'within_assay_summary.csv',cache/name/'within_assay_summary.csv',receipts)
        absolute['source_label']=label;native.append(absolute)
        for regime,regime_label in REGIMES:
            native_path=node/regime/'native_matched_common'/'paired_estimates.csv'
            native_table=snapshot(native_path,cache/name/regime/'native_matched_paired_estimates.csv',receipts)
            native_table=native_table[native_table.model.str.startswith('native_')].copy()
            native_table['source_label']=label;native_table['source_key']=key;native_table['regime']=regime
            native_table['comparison_support']='native_matched_common';native_decisions.append(native_table)
            tables={}
            for comparison in (['learned_common'] if args.figure=='broad' else ['learned_common','native_matched_common']):
                file=node/regime/comparison/'paired_uncertainty.csv'
                full=snapshot(file,cache/name/regime/(comparison+'_paired_uncertainty.csv'),receipts)
                full['source_label']=label;full['source_key']=key;full['regime']=regime;full['comparison_support']=comparison
                full['input_identities']=description['input_identities'];full['experimental_eligible_identities']=description['evaluated_identities']
                sensitivity.append(full);tables[comparison]=full
            full=tables['learned_common' if args.figure=='broad' else 'native_matched_common']
            selected=full[full.domain.eq('all_source_operational') & full.weighting.eq('actual_assay_pair_equal') &
                          full.resampling.eq('chemical_group_fixed_panel') & full.model.isin(arms) & full.metric.isin([m for m,_ in METRICS])].copy()
            if len(selected)!=len(arms)*len(METRICS):raise ValueError(f'Missing or duplicated plotted rows: {file}')
            if not selected.n_resamples.eq(400).all():raise ValueError(f'Expected 400 registered draws: {file}')
            if not np.allclose(selected.estimate-selected.baseline_estimate,selected.arm_minus_baseline,atol=1e-12,rtol=0):raise ValueError('Paired delta mismatch')
            rows.append(selected)
    return pd.concat(rows,ignore_index=True),pd.concat(native,ignore_index=True),pd.concat(native_decisions,ignore_index=True),pd.concat(sensitivity,ignore_index=True),receipts


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--figure',choices=['broad','spd'],required=True)
    p.add_argument('--results-root',type=Path,default=Path('outputs/profile_postprocessing_run02'),help='Fresh replay output root, relative to the working directory')
    p.add_argument('--prediction-root',type=Path,help='Core prediction dataset root (e.g. outputs/core_run02); avoids dependence on dataset_state source paths')
    p.add_argument('--from-cache',action='store_true',help='Portable replay from accompanying consolidated CSVs; no upstream reads')
    args=p.parse_args()
    specs={'dockstring':'DOCKSTRING','own44':'Docking44'} if args.figure=='broad' else SPD
    arms=['full_profile','target_residual','global_plus_selected_column'] if args.figure=='broad' else ['full_profile','target_residual','native_score']
    name='figure3_broad_paired_value' if args.figure=='broad' else 'figure4_spd_paired_value'
    if args.from_cache:
        receipts=json.loads((R/(name+'_render_receipt.json')).read_text())['inputs']
        read=lambda suffix:pd.read_csv(R/'data'/(name+suffix),float_precision='round_trip')
        values,absolute,native_decisions,sensitivity=[read(suffix) for suffix in
            ['_values.csv','_within_assay_absolute.csv','_native_decision_absolute.csv','_all_weightings_domains.csv']]
    else:
        values,absolute,native_decisions,sensitivity,receipts=collect(args,specs,arms)
    values.to_csv(R/'data'/(name+'_values.csv'),index=False)
    absolute.to_csv(R/'data'/(name+'_within_assay_absolute.csv'),index=False)
    sensitivity.to_csv(R/'data'/(name+'_all_weightings_domains.csv'),index=False)
    # Explicit native rows retain their own common-support strata; never rank
    # native methods as though these cover the same molecules or assay pairs.
    absolute[absolute.model.str.startswith('native_')].to_csv(R/'data'/(name+'_native_absolute.csv'),index=False)
    native_decisions.to_csv(R/'data'/(name+'_native_decision_absolute.csv'),index=False)
    use_style();plt.rcParams.update({'font.size':7,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.0,'legend.fontsize':6.5})
    fig,axes=plt.subplots(2,2,figsize=(WIDE,6.2 if args.figure=='broad' else 6.6),sharex='col')
    fig.subplots_adjust(left=.29 if args.figure=='broad' else .25,right=.98,top=.94,bottom=.14,wspace=.90 if args.figure=='broad' else .72,hspace=.48)
    for i,(regime,regime_label) in enumerate(REGIMES):
        for j,(metric,xlabel) in enumerate(METRICS):
            ax=axes[i,j];sub=values[values.regime.eq(regime)&values.metric.eq(metric)]
            labels=[];ys=[]
            for source_index,(key,label) in enumerate(specs.items()):
                for arm_index,arm in enumerate(arms):
                    q=sub[sub.source_key.eq(key)&sub.model.eq(arm)].iloc[0]
                    if args.figure=='broad':
                        y=source_index*len(arms)+arm_index
                        labels.append(label+' · '+{'full_profile':'profile','target_residual':'target','global_plus_selected_column':'selected'}[arm]);ys.append(y)
                    else:y=source_index+(arm_index-1)*.22
                    color=ARMS[arm][1];marker=ARMS[arm][2]
                    ax.plot([100*q.paired_delta_q025,100*q.paired_delta_q975],[y,y],color=color,lw=1.15,zorder=2)
                    ax.scatter(100*q.arm_minus_baseline,y,color=color,marker=marker,s=23,edgecolors='white',linewidths=.35,zorder=3)
            if args.figure=='broad':ax.set_yticks(ys,labels);ax.set_ylim(5.5,-.5)
            else:ax.set_yticks(range(len(specs)),list(specs.values()));ax.set_ylim(len(specs)-.5,-.5)
            ax.axvline(0,color='#666666',lw=.8,ls='--',zorder=1);ax.set_xlabel(xlabel)
            ax.grid(axis='x',alpha=.25);ax.grid(axis='y',visible=False)
            ax.text(-.02,1.10,chr(97+2*i+j),transform=ax.transAxes,fontweight='bold',fontsize=9,ha='left',va='bottom')
            ax.text(.5,1.025,regime_label,transform=ax.transAxes,fontsize=6.5,ha='center',va='bottom')
            if args.figure=='broad':ax.axhline(2.5,color='#CCCCCC',lw=.55)
    fig.legend(handles=[Line2D([],[],color=ARMS[a][1],marker=ARMS[a][2],lw=1.15,label=ARMS[a][0]) for a in arms],
               loc='lower center',ncol=3,bbox_to_anchor=(.59,.025),frameon=False)
    save(fig,str(R/name));fig.savefig(R/(name+'.png'),dpi=250,facecolor='white');plt.close(fig)
    support=[]
    for (key,regime),group in values.groupby(['source_key','regime'],sort=False):
        row=group.iloc[0]
        # Arm and metric supports must be the same on the primary learned-common rows.
        for col in ['n_rectangles','n_actual_assay_pairs','n_molecules','n_identified_D','n_pareto','n_chemical_groups']:
            if group[col].nunique()!=1:raise ValueError(f'Primary supports differ: {key}/{regime}/{col}')
        support.append(dict(source=key,regime=regime,comparison_support='learned_common' if args.figure=='broad' else 'native_matched_common',**{k:int(row[k]) for k in ['n_rectangles','n_actual_assay_pairs','n_molecules','n_identified_D','n_pareto','n_chemical_groups','input_identities','experimental_eligible_identities']}))
    pd.DataFrame(support).to_csv(R/'data'/(name+'_support.csv'),index=False)
    caption=('**Paired changes in experimental decision accuracy after adding score-profile information.** '
      'Panels a,b show new molecules with known proteins; c,d show new proteins with known molecules. '
      'Left: agreement with the interval-identified direction of D, the difference between two molecular changes at distinct target assays; predicted ties receive half credit. '
      'Right: the joint choice is correct only when both predicted molecular changes have the observed opposite directions. '
      'Each point is the indicated arm minus the global-level-plus-chemistry baseline on common evaluated rectangles, in percentage points (pp); horizontal lines are the 2.5th–97.5th percentiles of 400 paired chemical-cluster resampling draws with target-assay contexts fixed. '
      'Actual assay pairs receive equal weight in the all-source operational IC50 contrast domain. Source variants retain the original identity weighting; source records and bounds are not pooled into a thermodynamic selectivity estimate. '
      'These intervals preserve molecule dependencies and are not rectangle-binomial intervals or a declared success threshold. '
      'The full-profile arm adds the full centered profile. The matching-target residual arm adds a source-matched column residual. ')
    if args.figure=='broad':
        caption+=('All displayed arms are fitted readouts on learned-common support. The selected-column residual arm uses training interval-order association for known proteins and the sequence-nearest available column for held-out proteins, without held-out labels. The latter can coincide with the matching-target residual when the nearest source column is the native match; equal values are not independent evidence. ')
    if args.figure=='spd':
        if values.input_identities.nunique()!=1 or values.experimental_eligible_identities.nunique()!=1:raise ValueError('SPD source input cohorts differ; caption requires explicit per-readout wording')
        input_n=int(values.input_identities.iloc[0]);eligible_n=int(values.experimental_eligible_identities.iloc[0])
        geometry_n=int(pd.read_csv(R/'data/spd_readouts.csv').n_molecules.iloc[0])
        caption+=('The three plotted series use native-matched-common support within each readout: all arms and the baseline are evaluated on the same outcomes. Native scores are unfitted; the profile and residual arms are fitted readouts. Native scores are not assigned calibrated RMSE. ')
        caption+=('All eight readouts are shown in their prespecified order. Boltzina shares Boltz-2 lineage, and the two GNINA readouts share a scoring architecture; these are not eight independent method replications. '
          f'The source score input contains {input_n} identities per readout; {eligible_n} identities have eligible experimental observations. These are distinct from the {geometry_n} identities complete across all twelve channels in the descriptive geometry analysis. Actual decision support is reported separately. '
          'Native availability and resulting decision support differ between readouts; this is not a method ranking on an identical complete cohort. The full learned-common comparison is preserved in the accompanying all-weightings/domains table, separately labelled. '
          'Comparison with the descriptive geometry participation ratios is qualitative: no correlation significance test or biological-dimensionality claim is made. ')
    caption+=('Exact plotted estimates, paired intervals, baseline accuracies, all weighting/domain sensitivities, support and within-assay absolute accuracy are supplied as CSV alongside this figure. '
      'The support table gives scored experimental identities, actual assay pairs and rectangles separately; rectangle counts are not independent sample sizes. No best or worst condition was selected.\n')
    (R/(name+'_caption.md')).write_text(caption)
    reconciled=R/'refreshed_final_figure_sources_verification.json';dispatch=R/'data/final_postprocessing_dispatch_receipt.json'
    receipt=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),registered_postprocessing_run='c133ae0c-9f6a-4da6-a321-7ed9f79e585a',
                 scope='Display of final completed registered postprocessing. Original source-cache provenance is retained and reconciled exactly against final run02; no new computation of scientific estimates.',
                 reconciliation={'file':str(reconciled.relative_to(R)),'sha256':sha(reconciled)},
                 final_dispatch={'file':str(dispatch.relative_to(R)),'sha256':sha(dispatch)},
                 inputs=receipts,selected_rows=len(values),support=support,
                 outputs={str(p.relative_to(R)):sha(p) for p in [R/(name+'.pdf'),R/(name+'.svg'),R/(name+'.png'),R/'data'/(name+'_values.csv')]})
    (R/(name+'_render_receipt.json')).write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'rendered':name,'rows':len(values),'support':support}))


if __name__=='__main__':main()
