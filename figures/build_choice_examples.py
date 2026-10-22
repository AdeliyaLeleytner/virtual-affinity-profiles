"""Post hoc deterministic measured-choice illustrations; no model fit.

Prepare from the complete saved census, then render only cached four-cell
records and every source-variant decision. The Nesso-only native change uses
its documented nominal pIC50 identity; heterogeneous methods are not converted.
"""
from pathlib import Path
from itertools import product
import argparse
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd

R=Path(__file__).resolve().parent
D=R/'data/figure6_choices'
A=R/'assets/figure6_choices'
REGIME_A='new_molecules_known_proteins'
REGIME_B='new_proteins_known_molecules'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(input_dir,prediction_root,postprocessing_root,data_loader):
    D.mkdir(parents=True,exist_ok=True);sources=[];selection=[];cases=[];cells=[];variants=[]
    def source(path):
        digest=sha(path);sources.append(dict(path=str(path),sha256=digest));return path
    base=postprocessing_root/'core_run02'
    specifications=[('A','dockstring','dockstring',REGIME_A,'learned_common'),
                    ('A','own44','own44',REGIME_A,'learned_common'),
                    ('B','spd','nesso_affinity',REGIME_B,'native_matched_common')]
    candidates={'A':[],'B':[]}
    for case,dataset,channel,regime,support in specifications:
        node=f'{dataset}__{channel}__native_genes__all_records'
        path=source(base/node/regime/support/'continuous_rectangle_data.csv.gz')
        cols=['split','assay_id_a','assay_id_b','gene_a','gene_b','molecule_a','molecule_b','molecule_id_a','molecule_id_b','chemical_group_a','chemical_group_b','chemical_tanimoto','pareto_identified','full_profile__pareto_correct','global_chemistry__pareto_correct','full_profile__D_sign_credit','global_chemistry__D_sign_credit']
        if case=='B':cols+=['native_score__pareto_correct','native_score__D_sign_credit']
        frame=pd.read_csv(path,usecols=cols,float_precision='round_trip');total=len(frame)
        if case=='A':
            eligible=frame.pareto_identified & frame.full_profile__pareto_correct.eq(1) & frame.global_chemistry__pareto_correct.eq(0)
            rule='Observed opposite directions; full profile correct on all source-variant mass, global baseline correct on none.'
        else:
            eligible=(frame.pareto_identified & frame.native_score__pareto_correct.eq(1) & frame.full_profile__pareto_correct.eq(0) & frame.native_score__D_sign_credit.eq(1) & frame.full_profile__D_sign_credit.eq(1))
            rule='Observed opposite directions; native correct and full profile incorrect, with both D directions correct on all source-variant mass.'
        chosen=frame[eligible].copy();chosen['same_chemical_group']=chosen.chemical_group_a.eq(chosen.chemical_group_b)
        chosen['dataset']=dataset;chosen['channel']=channel;chosen['regime']=regime;chosen['comparison_support']=support;chosen['source_node']=node;chosen['case']=case
        candidates[case].append(chosen)
        selection.append(dict(case=case,dataset=dataset,channel=channel,census_rows=total,logical_category_rows=len(chosen),same_chemical_group_rows=int(chosen.same_chemical_group.sum()),rule=rule,source_sha256=sha(path)))
    # This order is specified before examining numerical effect magnitudes.
    keys=['same_chemical_group','chemical_tanimoto','molecule_id_a','molecule_id_b','assay_id_a','assay_id_b','split','dataset']
    ascending=[False,False,True,True,True,True,True,True]
    spec=importlib.util.spec_from_file_location('choice_source_data',source(data_loader))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    for case in ['A','B']:
        pool=pd.concat(candidates[case],ignore_index=True).sort_values(keys,ascending=ascending,kind='stable')
        if pool.empty:
            selection.append(dict(case=case,empty_specified_category=True));continue
        row=pool.iloc[0].to_dict();selection.append(dict(case=case,pooled_eligible_rows=len(pool),ordering=keys,ascending=ascending,selected_identity_a=row['molecule_id_a'],selected_identity_b=row['molecule_id_b'],selected_tanimoto=row['chemical_tanimoto'],selected_same_chemical_group=bool(row['same_chemical_group']),selected_assays=[int(row['assay_id_a']),int(row['assay_id_b'])]))
        data=module.load(input_dir,row['dataset'],row['channel'] if row['dataset']=='spd' else None,'native_genes','all_records')
        for name in ['molecules.csv','assays.csv','binding_records.csv.gz','molecule_features.npz','proteins.csv','protein_features.npz',f'datasets/{row["dataset"]}.npz',f'datasets/{row["dataset"]}_rows.csv']:
            source(input_dir/name)
        rows=pd.read_csv(input_dir/'datasets'/f'{row["dataset"]}_rows.csv').set_index('row_variant_index')
        assay_map=dict(zip(data['assays'].assay_id,data['assays'].index))
        aa,ab=[int(assay_map[row[key]]) for key in ['assay_id_a','assay_id_b']]
        for end in ['a','b']:
            molecule=int(row['molecule_'+end]);m=data['molecules'].iloc[molecule]
            row['name_'+end]=m['name'];row['smiles_'+end]=m.reference_smiles
            for target,index in [('t',aa),('u',ab)]:
                records=data['records'][data['records'].molecule_index.eq(molecule)&data['records'].local_assay_index.eq(index)]
                if len(records)!=1:raise ValueError('Four source cells must each identify one retained actual-assay record')
                record=records.iloc[0]
                cells.append(dict(case=case,molecule_end=end,molecule_index=molecule,molecule_id=m.standard_molecule_id,molecule_name=m['name'],target=target,gene=data['assays'].iloc[index].gene,assay_id=int(data['assays'].iloc[index].assay_id),**{key:record[key] for key in ['ic50_uM','source_relation','source_row_id','lower','upper','lower_open','upper_open']}))
        raw_folder=prediction_root/row['source_node']/'core'/row['regime']/row['split']
        path=source(raw_folder/'predictions.npz');source(raw_folder/'split.json');complete=json.loads(source(raw_folder/'complete.json').read_text())
        if complete['artifact_sha256']['predictions.npz']!=sha(path):raise ValueError('Saved predictions differ from completion receipt')
        saved=np.load(path,allow_pickle=False);test=set(map(int,saved['test_rows']))
        names=[name for name in saved.files if name not in ['test_rows','test_assays'] and (row['comparison_support']=='native_matched_common' or not name.startswith('native_'))]
        mass=0.;pair_data=[]
        for va,vb in product(data['identity_rows'][int(row['molecule_a'])],data['identity_rows'][int(row['molecule_b'])]):
            if int(va) not in test or int(vb) not in test:continue
            if not all(np.isfinite(saved[name][np.ix_([va,vb],[aa,ab])]).all() for name in names):continue
            weight=float(data['variant_weights'][va]*data['variant_weights'][vb]);mass+=weight
            pair_data.append((int(va),int(vb),weight))
        for va,vb,weight in pair_data:
            for model in names:
                pa,pb=saved[model][vb,aa]-saved[model][va,aa],saved[model][vb,ab]-saved[model][va,ab]
                def representation(index):
                    rec=rows.loc[index]
                    return str(rec.get('source_smiles',rec.get('original_source_smiles',rec.get('reference_smiles',''))))
                variants.append(dict(case=case,model=model,source_variant_row_a=va,source_variant_row_b=vb,source_representation_a=representation(va),source_representation_b=representation(vb),source_variant_weight_a=float(data['variant_weights'][va]),source_variant_weight_b=float(data['variant_weights'][vb]),conditional_pair_weight=weight/mass,predicted_A_t=float(saved[model][va,aa]),predicted_A_u=float(saved[model][va,ab]),predicted_B_t=float(saved[model][vb,aa]),predicted_B_u=float(saved[model][vb,ab]),delta_t=float(pa),delta_u=float(pb),D=float(pa-pb),D_direction=int(np.sign(pa-pb)),units='native source-score units' if model.startswith('native_') else 'pIC50'))
        row['source_variant_pair_mass']=mass;row['n_source_variant_pairs']=len(pair_data);row['selection_status']='post hoc descriptive illustration, not representative effect estimate'
        cases.append(row)
    pd.DataFrame(cases).to_csv(D/'selected_cases.csv',index=False)
    pd.DataFrame(cells).to_csv(D/'actual_four_cells.csv',index=False)
    pd.DataFrame(variants).to_csv(D/'all_source_variant_decisions.csv',index=False)
    (D/'selection_receipt.json').write_text(json.dumps(dict(scope='Measured-outcome illustration; no new performance estimate or predictive fit.',post_hoc=True,selection_steps=selection,priority='Within the logical category: same inherited chemical group first, then highest cached chirality-aware radius-3 Morgan/2048-bit Tanimoto, then lexical identity/assay/split/source keys. No effect-magnitude ranking or minimum-effect gate.',source_files=sources),indent=2,default=str)+'\n')
    print(pd.DataFrame(cases)[['case','dataset','channel','name_a','name_b','gene_a','gene_b','assay_id_a','assay_id_b','same_chemical_group','chemical_tanimoto','n_source_variant_pairs']].to_string(index=False),flush=True)


def regenerate_molecule_assets(cases):
    """Optional, separately versioned RDKit regeneration of four cached icons."""
    from rdkit import Chem,rdBase
    from rdkit.Chem.Draw import rdMolDraw2D
    A.mkdir(parents=True,exist_ok=True)
    for _,row in cases.iterrows():
        for end in ['a','b']:
            mol=Chem.MolFromSmiles(row['smiles_'+end])
            if mol is None:raise ValueError('Invalid prepared molecular SMILES')
            for format_name in ['svg','png']:
                drawer=rdMolDraw2D.MolDraw2DSVG(600,300) if format_name=='svg' else rdMolDraw2D.MolDraw2DCairo(900,450)
                drawer.drawOptions().padding=.06;drawer.drawOptions().bondLineWidth=2;drawer.DrawMolecule(mol);drawer.FinishDrawing()
                content=drawer.GetDrawingText();path=A/(row.case+'_'+end+'.'+format_name)
                path.write_text(content) if isinstance(content,str) else path.write_bytes(content)
    receipt=dict(generator='RDKit MolDraw2DSVG and MolDraw2DCairo; cached prepared reference SMILES',rdkit=rdBase.rdkitVersion,source_selected_cases_sha256=sha(D/'selected_cases.csv'),reference_smiles={row.case+'_'+end:row['smiles_'+end] for _,row in cases.iterrows() for end in ['a','b']},files={p.name:sha(p) for p in A.iterdir() if p.suffix in ['.png','.svg']})
    (A/'molecular_depiction_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


def render(regenerate_molecules=False):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle,Patch
    from PIL import Image
    from orx_figstyle import use_style,save,PALETTE,WIDE
    cases=pd.read_csv(D/'selected_cases.csv',float_precision='round_trip')
    cells=pd.read_csv(D/'actual_four_cells.csv',float_precision='round_trip')
    decisions=pd.read_csv(D/'all_source_variant_decisions.csv',float_precision='round_trip')
    unit_identity=json.loads((D/'nesso_unit_identity.json').read_text())
    if regenerate_molecules:regenerate_molecule_assets(cases)
    depiction=json.loads((A/'molecular_depiction_receipt.json').read_text())
    if depiction['source_selected_cases_sha256']!=sha(D/'selected_cases.csv'):raise ValueError('Cached molecular icons do not match the selected-case table')
    for filename,digest in depiction['files'].items():
        if sha(A/filename)!=digest:raise ValueError('Cached molecular icon changed; use explicit optional regeneration')
    A.mkdir(parents=True,exist_ok=True);use_style();plt.rcParams.update({'font.size':7,'axes.labelsize':6.5,'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':6})
    fig=plt.figure(figsize=(WIDE,6.65));intervals=[];caption=[]
    def relation_value(row):return str(row.source_relation)+format(row.ic50_uM,'.12g')
    def axes_style(ax,xlim,ylim,gene_t,gene_u,native=False):
        ax.set_xlim(xlim);ax.set_ylim(ylim);ax.set_aspect('equal',adjustable='box');ax.axhline(0,color='#444444',lw=.65);ax.axvline(0,color='#444444',lw=.65)
        low=max(xlim[0],ylim[0]);high=min(xlim[1],ylim[1]);ax.plot([low,high],[low,high],color='#888888',ls='--',lw=.7)
        quantity='source score' if native else 'pIC50'
        ax.set_xlabel(f'Δ {quantity} ({gene_t})'+('\n(native units)' if native else ''))
        ax.set_ylabel(f'Δ {quantity} ({gene_u})'+('\n(native units)' if native else ''))
        ax.grid(False)
    for index,case in enumerate(['A','B']):
        x=.045 if index==0 else .555
        if case not in set(cases.case):
            fig.text(x,.94,case+'  No case in the specified category',fontsize=8);continue
        row=cases[cases.case.eq(case)].iloc[0];records=cells[cells.case.eq(case)]
        var=decisions[decisions.case.eq(case)]
        if not var.groupby('model').conditional_pair_weight.sum().apply(lambda v:np.isclose(v,1,atol=1e-12,rtol=0)).all():raise ValueError('Conditional source-variant weights do not sum to one')
        title='DOCKSTRING · new molecules' if case=='A' else 'Nesso-1 · new proteins'
        fig.text(x,.965,case.lower()+'  '+title,fontsize=8,fontweight='bold',va='top')
        for offset,end in [(0,'a'),(.235,'b')]:
            ax=fig.add_axes([x+offset,.79,.20,.14]);ax.imshow(Image.open(A/(case+'_'+end+'.png')));ax.axis('off')
            fig.text(x+offset+.10,.778,end.upper()+': '+str(row['name_'+end]),ha='center',fontsize=6.7)
        fig.text(x+.219,.858,'→',ha='center',va='center',fontsize=12,color='#666666')
        fig.text(x+.219,.745,f'Replacement A → B; Morgan Tanimoto {row.chemical_tanimoto:.6f}',ha='center',fontsize=6.3)
        table_ax=fig.add_axes([x,.589,.435,.132]);table_ax.axis('off')
        body=[]
        for end in ['a','b']:
            body.append([relation_value(records[(records.molecule_end.eq(end))&records.target.eq(t)].iloc[0]) for t in ['t','u']])
        labels=[f'{row.gene_a}\nassay {int(row.assay_id_a)}',f'{row.gene_b}\nassay {int(row.assay_id_b)}']
        table=table_ax.table(cellText=body,rowLabels=['A','B'],colLabels=labels,cellLoc='center',loc='center',bbox=[.10,0,.90,.91]);table.auto_set_font_size(False);table.set_fontsize(6.5)
        for (r,c),cell in table.get_celld().items():cell.set_edgecolor('#CCCCCC');cell.set_linewidth(.4);cell.set_facecolor('#F2F2F2' if r==0 else 'white')
        fig.text(x+.219,.727,'Actual source IC50 (μM; = or > retained)',ha='center',fontsize=6.3)
        by={(r.molecule_end,r.target):r for r in records.itertuples()}
        bounds=[]
        for target in ['t','u']:
            a,b=by['a',target],by['b',target]
            bounds.append((b.lower-a.upper,b.upper-a.lower,bool(b.lower_open or a.upper_open),bool(b.upper_open or a.lower_open)))
        intervals.append(dict(case=case,delta_t_lower=bounds[0][0],delta_t_upper=bounds[0][1],delta_t_lower_open=bounds[0][2],delta_t_upper_open=bounds[0][3],delta_u_lower=bounds[1][0],delta_u_upper=bounds[1][1],delta_u_lower_open=bounds[1][2],delta_u_upper_open=bounds[1][3]))
        if bounds[0][0]!=bounds[0][1]:raise ValueError('Selected illustration needs generic two-dimensional bound-region renderer')
        if case=='A':
            ax=fig.add_axes([.105,.19,.35,.35]);axes_style(ax,(-.8,.6),(-.25,1.15),row.gene_a,row.gene_b)
        else:
            if row.channel!='nesso_affinity':raise ValueError('Common nominal pIC50 change axis is documented only for this Nesso case')
            ax=fig.add_axes([.615,.19,.35,.35]);axes_style(ax,(-.36,.36),(-.36,.36),row.gene_a,row.gene_b)
        # Opposite-sign quadrant is determined by measured bounds, not a utility threshold.
        sign_t=1 if bounds[0][0]>0 else -1;sign_u=1 if bounds[1][0]>=0 else -1
        xmin,xmax=ax.get_xlim();ymin,ymax=ax.get_ylim();qx0,qx1=(0,xmax) if sign_t>0 else (xmin,0);qy0,qy1=(0,ymax) if sign_u>0 else (ymin,0)
        ax.add_patch(Rectangle((qx0,qy0),qx1-qx0,qy1-qy0,facecolor='#EEF3F6',edgecolor='none',zorder=0))
        observed_x=bounds[0][0];lo,hi,_,_=bounds[1]
        boundary=lo if np.isfinite(lo) else hi
        ax.scatter(observed_x,boundary,marker='o',s=28,facecolors='white',edgecolors='black',lw=.9,zorder=6)
        arrow_end=ymax-.04*(ymax-ymin) if np.isposinf(hi) else ymin+.04*(ymax-ymin)
        ax.annotate('',xy=(observed_x,arrow_end),xytext=(observed_x,boundary),arrowprops={'arrowstyle':'->','color':'black','lw':1.0})
        models=['global_chemistry','full_profile'] if case=='A' else ['full_profile','native_score']
        for model in models:
            selected=var[var.model.eq(model)];color={'full_profile':PALETTE['blue'],'native_score':PALETTE['green'],'global_chemistry':'#666666'}[model];marker={'full_profile':'o','native_score':'^','global_chemistry':'D'}[model]
            for record in selected.itertuples():ax.scatter(record.delta_t,record.delta_u,s=32,marker=marker,color=color,edgecolors='white',linewidths=.4,zorder=5)
        ax.text(.03,.98,'Calibrated predictions' if case=='A' else 'Predicted changes',transform=ax.transAxes,ha='left',va='top',fontsize=6)
        fig.text(x+.219,.112,'D = Δt − Δu; the diagonal has D = 0',ha='center',fontsize=6.3)
        caption.append(f"Case {case}: {row['name_a']} → {row['name_b']}, {row.gene_a}/{row.gene_b}, assays {int(row.assay_id_a)}/{int(row.assay_id_b)}, Morgan Tanimoto {row.chemical_tanimoto:.12f}, {int(row.n_source_variant_pairs)} source-variant pair with conditional weight 1. Observed Δt={bounds[0][0]:.12g}; observed Δu {'>' if np.isposinf(hi) else '<'}{boundary:.12g} pIC50.")
        for model in models:
            for record in var[var.model.eq(model)].itertuples():
                units='nominal predicted pIC50 changes' if model=='native_score' else record.units
                caption.append(f'{model}: Δt={record.delta_t:.12g}, Δu={record.delta_u:.12g}, D={record.D:.12g} ({units}); weight={record.conditional_pair_weight:.12g}.')
    handles=[Line2D([],[],marker='o',ls='',markerfacecolor='white',markeredgecolor='black',label='Observed open bound'),Line2D([],[],marker='o',ls='',color=PALETTE['blue'],label='Full profile'),Line2D([],[],marker='D',ls='',color='#666666',label='Global + chemistry'),Line2D([],[],marker='^',ls='',color=PALETTE['green'],label='Native Nesso (nominal pIC50)'),Line2D([],[],color='#888888',ls='--',lw=.7,label='D = 0'),Patch(facecolor='#EEF3F6',label='Observed opposing-sign quadrant')]
    fig.legend(handles=handles,loc='lower center',ncol=3,bbox_to_anchor=(.52,.03),frameon=False)
    name='figure6_measured_choice_examples';save(fig,str(R/name));fig.savefig(R/(name+'.png'),dpi=250,facecolor='white');plt.close(fig)
    pd.DataFrame(intervals).to_csv(D/'observed_delta_bounds.csv',index=False)
    text=('**The sign of a relative change does not identify both underlying molecular changes.** The source records compare replacements A→B at two actual assays. IC50 values retain their source relation and μM units. Higher pIC50 means lower IC50; Δ is B minus A. Open circles mark strict observed boundaries and arrows show the unbounded direction, not an invented point observation. Shading marks the opposing-sign quadrant identified by the measurements. The diagonal is D=Δt−Δu=0. All variant-specific model coordinates and conditional weights are retained; the selected cases each have one pair, so no prediction is averaged over source variants. Nesso documentation defines its raw affinity_pred_value as log10(IC50/μM). The cached oriented score is its negative, equal algebraically to nominal predicted pIC50 minus 6. The additive constant cancels in B−A, so Nesso native differences are exactly nominal predicted pIC50 differences and share the case-B axes. This uses the documented output definition and adds no empirical recalibration. No pIC50 unit assignment is made here for Vina, GNINA, AutoDock4 or other native readouts. ')
    text+=' '.join(caption)
    text+=(' These are post hoc descriptive illustrations. Case A was selected from both broad sources among observed opposing changes correctly predicted by the full profile and missed by the global baseline. Case B was selected from the Nesso-1 cold-protein census where native and full profile both predict the correct D direction but only native predicts both opposing changes. Within each specified category, selection prioritizes the inherited chemical group, then highest cached radius-3 chirality-aware Morgan/2048-bit Tanimoto, then lexical identity, assay, split and source keys; no effect-size ranking or minimum-effect gate was used. The replacement arrow compares recorded molecules and does not assert a medicinal-chemistry transformation. The assays define operational IC50 observations, not thermodynamic selectivity or clinical utility. The small full-profile Δt in A is shown exactly rather than assigned a practical-effect threshold. Numerical axes and labels remain vector; molecular pictures are original RDKit depictions of the cached prepared reference SMILES, with their SVG originals retained.\n')
    (R/(name+'_caption.md')).write_text(text)
    files=[D/'selected_cases.csv',D/'actual_four_cells.csv',D/'all_source_variant_decisions.csv',D/'observed_delta_bounds.csv',D/'selection_receipt.json',D/'nesso_unit_identity.json']
    (R/(name+'_render_receipt.json')).write_text(json.dumps(dict(scope='Measured-outcome illustration only; no new performance estimate or model fit.',molecule_rendering='Cached source icons; RDKit needed only for explicit optional regeneration',rdkit_asset_generator=depiction['rdkit'],inputs={str(p.relative_to(R)):sha(p) for p in files+[A/'molecular_depiction_receipt.json']},outputs={ext:sha(R/(name+'.'+ext)) for ext in ['pdf','svg','png']}),indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--mode',choices=['prepare','render','both'],default='render');parser.add_argument('--regenerate-molecules',action='store_true',help='Optional RDKit regeneration; normal rendering uses cached source icons')
    parser.add_argument('--input-dir',type=Path,help='Prepared public input directory, e.g. data')
    parser.add_argument('--prediction-root',type=Path,help='Core prediction dataset root, e.g. outputs/core_run02')
    parser.add_argument('--postprocessing-root',type=Path,help='Replay root containing core_run02, e.g. outputs/profile_postprocessing_run02')
    parser.add_argument('--data-loader',type=Path,help='Public data loader file, e.g. src/profile_information/data.py')
    args=parser.parse_args()
    if args.mode in ['prepare','both'] and any(value is None for value in [args.input_dir,args.prediction_root,args.postprocessing_root,args.data_loader]):parser.error('prepare/both requires explicit --input-dir, --prediction-root, --postprocessing-root and --data-loader; there is no private-workspace fallback')
    if args.mode in ['prepare','both']:prepare(args.input_dir,args.prediction_root,args.postprocessing_root,args.data_loader)
    if args.mode in ['render','both']:render(args.regenerate_molecules)


if __name__=='__main__':main()
