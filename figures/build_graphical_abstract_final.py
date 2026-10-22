"""Final 920×300 graphical abstract using cached actual molecules/protein assets.
The archived draft and its generating script are never modified by this script.
"""
from pathlib import Path
import hashlib,io,json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
from orx_figstyle import use_style,PALETTE

R=Path(__file__).resolve().parent;D=R/'data';A=R/'assets'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    use_style();plt.rcParams.update({'savefig.bbox':None,'font.size':7,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6})
    molecules=pd.read_csv(D/'molecule_illustrations.csv');profiles=pd.read_csv(D/'illustration_profiles.csv')
    proteins=pd.read_csv(D/'protein_illustrations.csv');inputs=[D/'molecule_illustrations.csv',D/'illustration_profiles.csv',D/'protein_illustrations.csv',A/'structure_retrieval.json']
    colors=[PALETTE['blue'],PALETTE['orange']]
    fig=plt.figure(figsize=(9.2,3),dpi=100,facecolor='white')
    fig.text(.5,.965,'When do virtual affinity profiles improve biological decisions?',ha='center',va='top',fontsize=12,color='#17212B')
    for x,text in [(.105,'Molecules'),(.35,'Target proteins'),(.59,'Global level + shape'),(.85,'Measured preference tests')]:
        fig.text(x,.835,text,ha='center',fontsize=8.5,fontweight='bold',color='#17212B')
    for x in [.218,.483,.712]:
        ax=fig.add_axes([x-.012,.48,.024,.10]);ax.axis('off');ax.annotate('',(.95,.5),(.05,.5),arrowprops={'arrowstyle':'->','color':'#777777','lw':1.3})
    for index,row in enumerate(molecules.itertuples()):
        path=A/(row.name+'_molecule.png');inputs.append(path)
        ax=fig.add_axes([.013,.54-index*.24,.185,.20]);ax.imshow(Image.open(path));ax.axis('off')
        fig.text(.105,.515-index*.24,row.name,ha='center',fontsize=7.5,color=colors[index])
    for x,pdb,gene,label in [(.242,'6CM4','DRD2','GPCR'),(.322,'6PDJ','LCK','Kinase'),(.402,'3B66','AR','Nuclear receptor')]:
        match=proteins[proteins.source_target_key.eq(pdb.lower())]
        if len(match)!=1 or match.nominal_gene.iloc[0]!=gene:raise ValueError('Illustrative source structure mapping changed')
        path=A/(pdb+'_backbone.png');inputs.extend([path,A/(pdb+'.pdb')])
        ax=fig.add_axes([x,.40,.075,.34]);ax.imshow(Image.open(path));ax.axis('off')
        fig.text(x+.0375,.375,gene,ha='center',fontsize=7.5)
        fig.text(x+.0375,.325,pdb,ha='center',fontsize=6.5,color='#555555')
        fig.text(x+.0375,.275,label,ha='center',fontsize=6.0,color='#555555')
    for index,(field,ylabel) in enumerate([('score','Score'),('residual','Score − mean')]):
        ax=fig.add_axes([.51,.56-index*.24,.17,.16]);ax.grid(False)
        for i,row in enumerate(molecules.itertuples()):
            q=profiles[profiles.name.eq(row.name)];ax.plot(q.column,q[field],color=colors[i],lw=.9)
            if field=='score':ax.axhline(q.row_mean.iloc[0],color=colors[i],lw=.6,ls='--',alpha=.6)
        if field=='residual':ax.axhline(0,color='#777777',lw=.5)
        ax.set_ylabel(ylabel,fontsize=6.5,labelpad=1);ax.set_yticks([]);ax.set_xticks([]);ax.spines[['top','right']].set_visible(False)
        if index==1:ax.set_xlabel('Virtual target columns',fontsize=6.5,labelpad=2)
    fig.text(.85,.665,'Chemistry + global level\nvs profile vs query-target score',ha='center',va='center',fontsize=8,linespacing=1.35)
    fig.text(.85,.475,'Predict measured\npreference changes',ha='center',va='center',fontsize=8,linespacing=1.35)
    fig.text(.85,.30,'New molecules / known proteins\nNew proteins / known molecules',ha='center',va='center',fontsize=7,linespacing=1.35,color='#555555')
    fig.text(.5,.145,'Panel breadth does not guarantee biological specificity.',ha='center',va='center',fontsize=10,fontweight='bold',color='#17212B')
    fig.text(.5,.055,'Keep the query-target signal and test the decision.',ha='center',va='center',fontsize=9.5,color='#17212B')
    name='graphical_abstract_final'
    for extension in ['pdf','svg']:fig.savefig(R/(name+'.'+extension),facecolor='white')
    buffer=io.BytesIO();fig.savefig(buffer,format='png',dpi=100,facecolor='white');plt.close(fig);buffer.seek(0)
    picture=Image.open(buffer).convert('RGB').resize((920,300),Image.Resampling.LANCZOS)
    picture.quantize(colors=192,method=Image.Quantize.MEDIANCUT).save(R/(name+'.png'),optimize=True)
    if (R/(name+'.png')).stat().st_size>150000:raise ValueError('Graphical abstract exceeds requested 150KB')
    caption=('**When do virtual affinity profiles improve biological decisions?** The final graphical abstract shows real prepared molecules (butenafine and indomethacin, RDKit depictions of cached reference SMILES), source-deposited protein illustrations and the named molecules’ actual Docking44 score profiles and row means. The mean/residual plots illustrate the representation; they are not result-selected success cases. The experiment compares chemistry plus global score level, profile shape and query-specific scores on measured preference changes in distinct new-molecule/known-protein and new-protein/known-molecule regimes. The conclusion is scoped to the evaluated score panels and measured decision endpoints; it does not assert universal absence of useful profile information. '
      'The protein illustrations are cached smoothed C-alpha backbone traces of source-deposited DRD2/6CM4, LCK/6PDJ and AR/3B66 coordinates, using the first model’s longest chain. They are illustrative receptor examples, not new predictions or attestations that the deposited chains equal the experimental assay constructs. Deposition/retrieval provenance and chosen chains remain in assets/structure_retrieval.json and graphical_abstract_render_receipt.json. Numerical chart primitives and composition text remain vector in PDF/SVG; molecular/protein picture layers are intentionally raster illustrations. The archived graphical_abstract_draft files are preserved.\n')
    (R/(name+'_caption.md')).write_text(caption)
    receipt=dict(scope='Final display using existing assets and frozen illustrative profiles; no new data, fitting, structure prediction or image generation.',canvas_pixels=[920,300],png_bytes=(R/(name+'.png')).stat().st_size,inputs={str(p.relative_to(R)):sha(p) for p in inputs},outputs={ext:sha(R/(name+'.'+ext)) for ext in ['pdf','svg','png']},font='DejaVu Sans')
    (R/(name+'_render_receipt.json')).write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(rendered=name,pixels=picture.size,png_bytes=receipt['png_bytes'])))


if __name__=='__main__':main()
