"""Plot released frozen growth curves without selecting a best k or fitting.
python build_growth_figure.py --analysis-dir PATH --model full_profile
Requires completed growth_analysis output and an ORX provenance receipt supplied
by root. All six predeclared order curves are shown; none is an uncertainty CI.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.lines import Line2D
from orx_figstyle import use_style, save, PALETTE, WIDE

R=Path(__file__).resolve().parent
METRICS=[('D_direction_concordance','D-direction concordance'),
         ('pareto_correct_fraction','Joint opposite-direction choice')]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-dir',type=Path,required=True)
    p.add_argument('--run-log','--orx-receipt',dest='orx_receipt',type=Path,required=True,help='Completed growth-analysis stdout log; --orx-receipt is the archival alias')
    p.add_argument('--model',choices=['full_profile','global_chemistry'],default='full_profile')
    p.add_argument('--output-name')
    args=p.parse_args()
    receipt_path=args.analysis_dir/'growth_analysis_receipt.json'
    receipt=json.loads(receipt_path.read_text())
    if not args.orx_receipt.is_file() or '"growth_analysis_complete": true' not in args.orx_receipt.read_text():raise ValueError('Completed ORX execution receipt is required')
    if receipt['frozen_sources_before']!=receipt['frozen_sources_after']:raise ValueError('Frozen source verification differs before/after')
    csv=args.analysis_dir/'all_growth_preferences.csv'
    frame=pd.read_csv(csv,float_precision='round_trip')
    primary=frame[frame.domain.eq('all_source_operational') & frame.weighting.eq('actual_assay_pair_equal') & frame.model.eq(args.model) & frame.comparison_support.eq('learned_common')]
    sources=list(dict.fromkeys(primary.source_label))
    if len(sources)!=2:raise ValueError('Figure requires exactly two released source panels')
    seeds=receipt['config'].get('growth_analysis',receipt['config'])['seeds']
    policies={'random':PALETTE['blue'],'family_stratified':PALETTE['cyan']}
    styles=['-','--',':']
    use_style();plt.rcParams.update({'font.size':7,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'legend.fontsize':6.5})
    fig,axes=plt.subplots(2,2,figsize=(WIDE,6.25),sharey='col')
    fig.subplots_adjust(left=.11,right=.84,bottom=.21,top=.95,wspace=.98,hspace=.48)
    caption=[];plotted=[];paired_values=[];cost_values=[];source_paths=[csv,receipt_path,args.orx_receipt];query_axes=[]
    for i,source in enumerate(sources):
        sub=primary[primary.source_label.eq(source)]
        dataset=sub.dataset.iloc[0]
        def table(filename):
            path=args.analysis_dir/dataset/filename;source_paths.append(path)
            return pd.read_csv(path,float_precision='round_trip')
        paired=table('all_measurements_paired_uncertainty.csv')
        paired=paired[paired.domain.eq('all_source_operational') & paired.weighting.eq('actual_assay_pair_equal') & paired.resampling.eq('chemical_group_fixed_panel')]
        if not paired.n_resamples.eq(400).all():raise ValueError('Expected registered 400 paired draws')
        invariance=table('full_budget_invariance_checks.csv')
        if not invariance.predictions_bitwise_identical.all():raise ValueError('Canonical full-budget invariance failed')
        pair_cost=table('query_assay_pair_costs.csv');union_cost=table('query_profile_union_costs.csv')
        query_models=['chemistry_plus_target_score','chemistry_plus_selected_score']
        for model in query_models:
            q=pair_cost[pair_cost.model.eq(model)&pair_cost.n_saved_rectangles.gt(0)];u=union_cost[union_cost.model.eq(model)]
            cost_values.append(dict(source=source,model=model,training_candidate_min=int(u.training_score_candidate_columns.min()),training_candidate_max=int(u.training_score_candidate_columns.max()),pair_query_min=int(q.deployment_query_column_count.min()),pair_query_max=int(q.deployment_query_column_count.max()),whole_profile_union_min=int(u.deployment_whole_profile_union_count.min()),whole_profile_union_max=int(u.deployment_whole_profile_union_count.max()),evaluated_contrast_union_min=int(u.deployment_evaluated_contrast_union_count.min()),evaluated_contrast_union_max=int(u.deployment_evaluated_contrast_union_count.max())))
        if sub.order.nunique()!=6:raise ValueError('All six predeclared order curves are required')
        for j,(metric,label) in enumerate(METRICS):
            ax=axes[i,j];anchor=sub[sub.k.eq(0)][metric]
            if anchor.nunique(dropna=False)!=1:raise ValueError('k0 anchor differs across order copies')
            baseline=anchor.iloc[0]
            for policy,color in policies.items():
                for seed,style in zip(seeds,styles):
                    q=sub[sub.order.eq(f'{policy}_{seed}')].sort_values('k')
                    if set(q.k)!=set(range(int(sub.k.max())+1)):raise ValueError('Missing k values')
                    gain=100*(q[metric]-baseline)
                    ax.plot(q.k,gain,color=color,ls=style,lw=1.0,alpha=.85)
                    for (_,row),value in zip(q.iterrows(),gain):
                        plotted.append(dict(source=source,model=args.model,metric=metric,order=row.order,k=row.k,
                                            estimate=row[metric],zero_score_baseline=baseline,gain_percentage_points=value,
                                            n_assay_pairs=row.n_assay_pairs,n_rectangles=row.n_rectangles,
                                            n_identified_D=row.n_identified_D,n_pareto=row.n_pareto))
            ax.axhline(0,color='#555555',lw=.9,zorder=1)
            ax.scatter([0],[0],s=25,color='#555555',zorder=4)
            query=ax.inset_axes([1.05,0,.35,1],sharey=ax);query_axes.append(query)
            for position,(model,color,marker) in enumerate(zip(query_models,[PALETTE['orange'],PALETTE['purple']],['s','D'])):
                q=paired[paired.model.eq(model)&paired.metric.eq(metric)]
                if len(q)!=1:raise ValueError('Missing query comparator')
                q=q.iloc[0]
                if not np.isclose(q.baseline_estimate,baseline,atol=1e-12,rtol=0):raise ValueError('Query baseline differs from k0')
                low,point,high=100*np.array([q.paired_delta_q025,q.arm_minus_baseline,q.paired_delta_q975])
                query.vlines(position,low,high,color=color,lw=1.1);query.scatter(position,point,color=color,marker=marker,s=22,edgecolors='white',linewidths=.3,zorder=4)
                paired_values.append(dict(source=source,**q.to_dict()))
            endpoint='full_panel_profile' if args.model=='full_profile' else 'full_panel_level_chemistry'
            q=paired[paired.model.eq(endpoint)&paired.metric.eq(metric)].iloc[0]
            low,point,high=100*np.array([q.paired_delta_q025,q.arm_minus_baseline,q.paired_delta_q975])
            ax.vlines(sub.k.max(),low,high,color=PALETTE['blue'],lw=1.3,zorder=4,clip_on=False);ax.scatter(sub.k.max(),point,s=22,color=PALETTE['blue'],zorder=5,clip_on=False)
            paired_values.append(dict(source=source,**q.to_dict()))
            query.axhline(0,color='#555555',lw=.7);query.grid(axis='y',alpha=.2);query.grid(axis='x',visible=False)
            query.set_xlim(-.55,1.55);query.set_xticks([0,1],['M','S'],fontsize=6);query.tick_params(axis='y',labelleft=False,left=False);query.set_xlabel('Separate\nquery inputs',fontsize=6)
            ax.set_xlim(0,sub.k.max());ax.set_xlabel('Shared-panel columns\nk (count)')
            ax.set_ylabel(f'Δ {label}\n(percentage points)')
            ax.text(-.02,1.06,f'{chr(97+2*i+j)}  {source}',transform=ax.transAxes,fontweight='bold',fontsize=8,ha='left',va='bottom')
            ax.grid(axis='y',alpha=.25);ax.grid(axis='x',visible=False)
            support=sub.iloc[0]
            caption.append(f'{source}, {label}: k=0 estimate {baseline:.6f}; {int(support.n_assay_pairs)} actual assay pairs, {int(support.n_rectangles)} rectangles, {int(support.n_identified_D)} identified D directions and {int(support.n_pareto)} observed opposite-direction choices.')
    arm_label='Full profile' if args.model=='full_profile' else 'Level + chemistry'
    handles=[Line2D([],[],color=c,lw=1.2,label=arm_label+' · '+{'random':'random order','family_stratified':'class-stratified'}[policy]) for policy,c in policies.items()]
    handles += [Line2D([],[],color='#777777',ls=style,lw=1,label=f'Seed {seed}') for seed,style in zip(seeds,styles)]
    handles += [Line2D([],[],color='#555555',marker='o',lw=.9,label='Zero-score chemistry anchor')]
    handles += [Line2D([],[],color=PALETTE['orange'],marker='s',ls='',label='M: matching-score query'),Line2D([],[],color=PALETTE['purple'],marker='D',ls='',label='S: selected-score query')]
    fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.54,.02))
    name=args.output_name or ('figure5_panel_acquisition' if args.model=='full_profile' else 'figure_supplement_growth_level_chemistry')
    save(fig,str(R/name));fig.savefig(R/(name+'.png'),dpi=250,facecolor='white');plt.close(fig)
    values_path=R/'data'/(name+'_values.csv');pd.DataFrame(plotted).to_csv(values_path,index=False)
    pd.DataFrame(paired_values).to_csv(R/'data'/(name+'_paired_values.csv'),index=False);pd.DataFrame(cost_values).to_csv(R/'data'/(name+'_costs.csv'),index=False)
    text='**Acquired-panel decision curves relative to a score-free chemistry model.** '
    text+=('All six predeclared column orders are shown for the '+args.model.replace('_',' ')+' arm: three random and three round-robin broad-class orders, with the stated seeds. Every integer k from zero through the full source panel is retained, using all five held-out chemical folds. Curves are changes from the identical k=0 prediction anchor, in percentage points. They are composition trajectories, not sampling confidence intervals; no best k, acquisition budget or scientific success threshold is selected. Actual assay pairs receive equal weight in the all-source operational IC50 contrast domain. The D metric counts the correctly identified algebraic direction (predicted ties receive half credit); the joint-choice metric requires both opposing molecular changes to have the correct directions. Assay contexts, source variants and measurement bounds follow the frozen evaluation. Source-specific assay pair and rectangle support stays constant across k, orders and the two fitted arms. ')
    text+=' '.join(caption)
    text+=' Query-specific matching-score and label-selected-score comparators are separate analyses, never k=1 points on this shared-panel curve. Selected-score training searches the full candidate score panel; actual single-task, pair-task and whole-profile deployment unions are recorded separately. '
    text+='The full-panel endpoint bar and separate query bars are 2.5th–97.5th percentile intervals from 400 paired chemical-cluster resampling draws, conditional on fixed assay contexts. Query bars are on separate categorical axes, not k=1 points. Canonical full-budget predictions are bitwise identical across all six orders. '
    for cost in cost_values:text+=f"{cost['source']} {cost['model']}: training candidates {cost['training_candidate_min']}–{cost['training_candidate_max']}; per evaluated distinct-target pair {cost['pair_query_min']}–{cost['pair_query_max']} columns; whole queried-profile union {cost['whole_profile_union_min']}–{cost['whole_profile_union_max']}. "
    text+='No query-specific selected column is described as a free training acquisition. All unrounded values, supports and acquisition accounts are retained.\n'
    (R/(name+'_caption.md')).write_text(text)
    cache=R/'data'/('figure5_sources' if args.model=='full_profile' else 'growth_level_sources');cache.mkdir(parents=True,exist_ok=True)
    # A relocated package can reuse one cached source tree for both arms.
    # External analysis outputs are still snapshotted; cached inputs are not duplicated.
    if args.analysis_dir.resolve().is_relative_to((R/'data').resolve()):cache=args.analysis_dir.resolve()
    hashes={}
    for path in source_paths:
        before=hashlib.sha256(path.read_bytes()).hexdigest();relative=path.relative_to(args.analysis_dir) if path.is_relative_to(args.analysis_dir) else Path(path.name);dest=cache/relative;dest.parent.mkdir(parents=True,exist_ok=True)
        if path.resolve()!=dest.resolve():shutil.copyfile(path,dest)
        if hashlib.sha256(path.read_bytes()).hexdigest()!=before or hashlib.sha256(dest.read_bytes()).hexdigest()!=before:raise ValueError('Figure source changed during snapshot')
        hashes[str(path)]={'cached':str(dest.relative_to(R)),'sha256':before}
    (R/(name+'_render_receipt.json')).write_text(json.dumps(dict(inputs=hashes,model=args.model,scope='Display only; no fitting, cutoff or best-k selection',vector_outputs=['pdf','svg']),indent=2)+'\n')


if __name__=='__main__':main()
