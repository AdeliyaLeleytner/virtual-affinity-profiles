"""Render portable, caption-led geometry figures from frozen CSV tables.
Run in the analysis environment: python build_geometry_figures.py
No model fitting, data retrieval, selection threshold or current results reading.
"""
from pathlib import Path
import hashlib,json,platform
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from orx_figstyle import use_style,save,PALETTE,WIDE
R=Path(__file__).resolve().parent;D=R/'data'
use_style();plt.rcParams.update({'font.size':7,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'legend.fontsize':6.5})
RAW=PALETTE['blue'];CENTER=PALETTE['orange'];DS=PALETTE['green'];OWN=PALETTE['purple']
def letter(ax,s):ax.text(-.02,1.04,s,transform=ax.transAxes,fontweight='bold',fontsize=9,ha='left',va='bottom')
def finish(fig,name):
 save(fig,str(R/name));fig.savefig(R/(name+'.png'),dpi=250,facecolor='white');plt.close(fig)
# 1. Eight correlated readouts; source composition; different support universes.
s=pd.read_csv(D/'spd_readouts.csv');names={'boltz_affinity':'Boltz-2 affinity','boltzina_vina_far_affinity':'Boltzina far','boltzina_vina_zero_affinity':'Boltzina zero','gnina':'GNINA on Vina pose','nesso_affinity':'Nesso-1 affinity','raw_ad4':'AutoDock4','raw_gnina_own':'GNINA own search','vina':'Vina'}
order=list(names);fig,axs=plt.subplots(2,2,figsize=(WIDE,5.65));fig.subplots_adjust(left=.22,right=.98,top=.95,bottom=.17,wspace=.35,hspace=.45)
for ax,metric,xlab in [(axs[0,0],'pr','Participation ratio (unitless)'),(axs[0,1],'pc1','Variance in PC1 (%)')]:
 for i,c in enumerate(order):
  q=s[s.channel.eq(c)].set_index('mode');v=[q.loc[k,metric]*(100 if metric=='pc1' else 1) for k in ['raw','molecule_centered']]
  ax.plot(v,[i,i],color='#888888',lw=.7,zorder=1);ax.scatter(v,[i,i],c=[RAW,CENTER],s=22,zorder=3,edgecolors='white',linewidths=.3)
 ax.set_yticks(range(len(order)),[names[c] for c in order] if metric=='pr' else []);ax.invert_yaxis();ax.set_xlabel(xlab);ax.grid(axis='x',alpha=.25);ax.grid(axis='y',visible=False);ax.set_xlim(left=0)
 if metric=='pc1':ax.set_xlim(0,100);ax.set_ylabel('Readout (same order as a)')
 letter(ax,'a' if metric=='pr' else 'b')
comp=pd.read_csv(D/'panel_composition.csv');cats=sorted(comp.broad_class.unique());ax=axs[1,0]
for delta,source,color in [(-.18,'DOCKSTRING',DS),(.18,'Docking44',OWN)]:
 v=comp[comp.source.eq(source)].set_index('broad_class').n_columns.reindex(cats,fill_value=0)
 bars=ax.barh(np.arange(len(cats))+delta,v,height=.34,color=color,edgecolor='white',lw=.4);ax.bar_label(bars,fmt='%d',padding=1,fontsize=5.5)
ax.set_yticks(range(len(cats)),cats);ax.invert_yaxis();ax.set_xlabel('Target columns (count)');ax.set_xlim(0,29);ax.grid(axis='y',visible=False);ax.grid(axis='x',alpha=.25);letter(ax,'c')
co=pd.read_csv(D/'cohort_support.csv');ax=axs[1,1]
for i,row in enumerate(co.itertuples()):
 for delta,k,color in [(-.17,'score_identities','#B5B5B5'),(.17,'binding_anchors','#333333')]:
  v=getattr(row,k);b=ax.barh(i+delta,v,height=.31,color=color);ax.bar_label(b,padding=2,fontsize=6)
ax.set_yticks(range(len(co)),co.cohort);ax.invert_yaxis();ax.set_xlabel('Normalized molecular identities (count)');ax.set_xlim(0,1050);ax.grid(axis='y',visible=False);ax.grid(axis='x',alpha=.25);letter(ax,'d')
handles=[Line2D([],[],marker='o',ls='',color=RAW,label='Raw score'),Line2D([],[],marker='o',ls='',color=CENTER,label='Row-centered score'),Patch(color=DS,label='DOCKSTRING'),Patch(color=OWN,label='Docking44'),Patch(color='#B5B5B5',label='Scored identities'),Patch(color='#333333',label='Binding anchors')]
fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.53,.015));finish(fig,'figure_geometry_readouts_support')
# 2. Mean agreement; full spectrum and PR; named same-gene agreement.
fig=plt.figure(figsize=(WIDE,5.8));gs=fig.add_gridspec(2,3,height_ratios=[1,1.1],left=.09,right=.96,bottom=.14,top=.95,wspace=.70,hspace=.66)
a=fig.add_subplot(gs[0,0]);q=pd.read_csv(D/'shared_means.csv');r=pd.read_csv(D/'shared_mean_correlations.csv',index_col=0).loc['DOCKSTRING','Docking44']
a.scatter(q.DOCKSTRING,q.Docking44,s=8,color=RAW,alpha=.45,edgecolors='none');a.set_xlabel('DOCKSTRING mean\n(oriented score units)');a.set_ylabel('Docking44 mean\n(oriented score units)');a.text(.04,.95,f'$r={r:.3f}$\n$n={len(q)}$',transform=a.transAxes,va='top',fontsize=7);a.grid(alpha=.2);letter(a,'a')
b=fig.add_subplot(gs[0,1]);sp=pd.read_csv(D/'full_spectra.csv');colors={'raw':RAW,'row_centered':CENTER,'rank_columns':'#777777','rank_centered':PALETTE['green']};labels={'raw':'Raw','row_centered':'Row-centered','rank_columns':'Column ranks','rank_centered':'Ranks, then centered'}
for view in colors:
 z=sp[sp.panel.eq('Joint_equal_columns')&sp.view.eq(view)];b.plot(z.component,z.cumulative*100,color=colors[view],ls='--' if view.startswith('rank') else '-',lw=1)
b.set_xlabel('Retained components (count)');b.set_ylabel('Cumulative variance (%)');b.set_xlim(1,102);b.set_ylim(0,102);b.grid(alpha=.2);letter(b,'b')
c=fig.add_subplot(gs[0,2]);su=pd.read_csv(D/'spectral_summary.csv');panels=['DOCKSTRING','Docking44','Joint_equal_columns']
for j,view in enumerate(colors):
 z=su[su.view.eq(view)].set_index('panel').loc[panels];c.scatter(z.pr,np.arange(3)+(j-1.5)*.14,s=18,c=colors[view],marker='s' if view.startswith('rank') else 'o',edgecolors='white',linewidth=.3)
 for i,v in enumerate(z.pr):
  if view in ['raw','row_centered']:c.annotate(f'{v:.2f}',(v,i+(j-1.5)*.14),xytext=(3,0),textcoords='offset points',va='center',fontsize=5.3,color=colors[view])
c.set_yticks(range(3),['DOCKSTRING','Docking44','Joint']);c.invert_yaxis();c.set_xlabel('Participation ratio\n(unitless)');c.set_xlim(0,29);c.grid(axis='x',alpha=.25);c.grid(axis='y',visible=False);letter(c,'c')
d=fig.add_subplot(gs[1,:]);shared=pd.read_csv(D/'shared_gene_main_figure.csv')
genes=sorted(shared.gene.unique());positions=np.arange(len(genes))
for i,gene in enumerate(genes):
 q=shared[shared.gene.eq(gene)].set_index('view')
 vals=[q.loc[v,'pearson'] for v in ['raw','row_centered']]
 d.plot([i,i],vals,color='#888888',lw=.8,zorder=1)
 d.scatter([i,i],vals,color=[RAW,CENTER],s=28,zorder=3,edgecolors='white',linewidths=.3)
d.set_xticks(positions,genes);d.set_xlim(-.6,len(genes)-.4);d.set_ylim(0,1)
d.set_xlabel('Shared nominal gene (source-specific receptor context)')
d.set_ylabel('Cross-source Pearson r\n(unitless)');d.grid(axis='y',alpha=.25);d.grid(axis='x',visible=False);letter(d,'d')
fig.legend(handles=[Line2D([],[],color=colors[v],lw=1,ls='--' if v.startswith('rank') else '-',label=labels[v]) for v in colors],loc='lower center',ncol=4,bbox_to_anchor=(.52,.005),frameon=False);finish(fig,'figure_geometry_shared_shape')
# Supplemental matrix: same complete identities, with explicit annotation strips.
meta=pd.read_csv(D/'column_metadata.csv').sort_values(['source','broad_class','nominal_gene'],kind='stable')
corr=pd.read_csv(D/'correlation_row_centered_complete.csv',index_col=0).loc[meta.column_id,meta.column_id]
classes=sorted(meta.broad_class.unique());class_colors=dict(zip(classes,list(PALETTE.values())+['#999999']))
fig=plt.figure(figsize=(WIDE,7.6));gs=fig.add_gridspec(2,2,width_ratios=[.023,1],height_ratios=[.023,1],left=.10,right=.84,bottom=.27,top=.91,wspace=.012,hspace=.012)
top=fig.add_subplot(gs[0,1]);left=fig.add_subplot(gs[1,0]);d=fig.add_subplot(gs[1,1]);rgba=np.array([class_colors[c] for c in meta.broad_class])
from matplotlib.colors import ListedColormap
strip=np.array([classes.index(c) for c in meta.broad_class]);cmap=ListedColormap([class_colors[c] for c in classes])
top.pcolormesh(strip[None,:],cmap=cmap,vmin=-.5,vmax=len(classes)-.5,rasterized=False);left.pcolormesh(strip[:,None],cmap=cmap,vmin=-.5,vmax=len(classes)-.5,rasterized=False);left.invert_yaxis()
for ax in [top,left]:ax.axis('off')
im=d.pcolormesh(corr.to_numpy(),cmap='RdBu_r',vmin=-1,vmax=1,shading='flat',rasterized=False);d.invert_yaxis();d.set_aspect('equal');d.grid(False)
at=0;middles=[];ticks=[]
for (src,cl),part in meta.groupby(['source','broad_class'],sort=False):
 middles.append(at+len(part)/2);ticks.append(f'{src}: {cl}');at+=len(part)
 d.axhline(at,color='#666666',lw=.35);d.axvline(at,color='#666666',lw=.35)
split=int(meta.source.eq('DOCKSTRING').sum())
for pos in [split]:d.axhline(pos,color='black',lw=1);d.axvline(pos,color='black',lw=1)
# Class strip and full class legend make numeric matrix indices interpretable.
d.set_xticks([0,25,50,75,102]);d.set_yticks([0,25,50,75,102]);d.tick_params(length=2,pad=2,labelsize=6)
d.set_xlabel('Target column (source, class, then nominal gene order)');d.set_ylabel('Target column (same order)')
top.text(split/2,2.5,'DOCKSTRING (58)',ha='center',va='center',fontsize=7)
top.text((split+len(meta))/2,2.5,'Docking44 (44)',ha='center',va='center',fontsize=7)
cb=fig.colorbar(im,ax=d,pad=.025,fraction=.025);cb.solids.set_rasterized(False);cb.set_label('Pearson r (unitless)',fontsize=7)
fig.legend(handles=[Patch(facecolor=class_colors[c],label=c) for c in classes],loc='lower center',ncol=3,bbox_to_anchor=(.48,.06),frameon=False,title='Broad protein class; same strip order on both axes')
finish(fig,'figure_supplement_residual_correlation')
receipt={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'style':'Vendored orx_figstyle.py; one use_style() initializer; portable DejaVu Sans selected explicitly; vector PDF/SVG, 250 dpi PNG previews.','scope':'Completed descriptive geometry only; deterministic estimates, no resampling intervals.','csv_inputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in D.glob('*.csv')}}
(R/'geometry_render_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
