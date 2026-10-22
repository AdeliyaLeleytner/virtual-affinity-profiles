"""Optional regeneration of cached molecular and source-coordinate pictures.
Uses retained reference SMILES and PDB files; final figure files are unchanged.
Requires the optional RDKit/SciPy environment, not the normal cached renderer.
No new molecular/protein predictions are produced.
"""
from pathlib import Path
import hashlib,json,platform,io
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.interpolate import make_interp_spline
from rdkit import Chem,rdBase
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
from orx_figstyle import use_style,PALETTE
R=Path(__file__).resolve().parent;D=R/'data';A=R/'assets'
use_style();plt.rcParams.update({'savefig.bbox':None,'font.size':7,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6})
# All annotation names and exact molecular representations come from prepared CSV.
m=pd.read_csv(D/'molecule_illustrations.csv');pro=pd.read_csv(D/'protein_illustrations.csv');profiles=pd.read_csv(D/'illustration_profiles.csv')
colors=[PALETTE['blue'],PALETTE['orange']]
for row in m.itertuples():
 mol=Chem.MolFromSmiles(row.reference_smiles);assert mol is not None
 drawer=rdMolDraw2D.MolDraw2DSVG(360,190);o=drawer.drawOptions();o.bondLineWidth=2;o.padding=.08;o.clearBackground=True;drawer.DrawMolecule(mol);drawer.FinishDrawing();(A/(row.name+'_molecule.svg')).write_text(drawer.GetDrawingText())
 drawer=rdMolDraw2D.MolDraw2DCairo(720,380);drawer.drawOptions().bondLineWidth=3;drawer.drawOptions().padding=.08;drawer.DrawMolecule(mol);drawer.FinishDrawing();(A/(row.name+'_molecule.png')).write_bytes(drawer.GetDrawingText())
# Protein picture: smoothed primary-coordinate backbone, split at chain gaps.
# No secondary-structure classifier, docking, sequence substitution or folding.
structures=[]
for pdb,gene,color in [('6CM4','DRD2',PALETTE['blue']),('6PDJ','LCK',PALETTE['green']),('3B66','AR',PALETTE['purple'])]:
 row=pro[pro.source_target_key.eq(pdb.lower())].iloc[0];assert row.nominal_gene==gene
 lines=(A/(pdb+'.pdb')).read_text().splitlines();chains={}
 for line in lines:
  if line.startswith('ENDMDL'):break
  if line.startswith('ATOM') and line[12:16].strip()=='CA' and line[16:17] in [' ','A']:
   chains.setdefault(line[21],[]).append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
 chain=max(chains,key=lambda c:len(chains[c]));xyz=np.array(chains[chain]);xyz-=xyz.mean(0);_,_,v=np.linalg.svd(xyz,full_matrices=False);xyz=xyz@v.T
 # Principal orientation is a deterministic display choice, not a fitted model.
 xyz=xyz[:,[1,0,2]];xyz[:,0]*=-1
 gaps=np.where(np.linalg.norm(np.diff(xyz,axis=0),axis=1)>5)[0]+1;chunks=np.split(xyz,gaps)
 fig,ax=plt.subplots(figsize=(1.2,1.8));fig.subplots_adjust(0,0,1,1);ax.set_aspect('equal');ax.axis('off')
 for chunk in chunks:
  if len(chunk)<4:continue
  t=np.arange(len(chunk));u=np.linspace(0,len(chunk)-1,len(chunk)*6);smooth=make_interp_spline(t,chunk,k=3)(u);segs=np.stack([smooth[:-1,:2],smooth[1:,:2]],axis=1);depth=(smooth[:-1,2]-xyz[:,2].min())/(np.ptp(xyz[:,2])+1e-12)
  from matplotlib.colors import to_rgb
  rgb=np.array(to_rgb(color));rgba=np.c_[np.clip(rgb[None,:]*(.6+.4*depth[:,None])+.15*(1-depth[:,None]),0,1),np.ones(len(depth))]
  # White contour separates crossing backbone strands; colored strokes are C-alpha traces.
  ax.add_collection(LineCollection(segs,colors='white',linewidths=2.1,capstyle='round',joinstyle='round'))
  ax.add_collection(LineCollection(segs,colors=rgba,linewidths=1.3,capstyle='round',joinstyle='round'))
 ax.autoscale();ax.margins(.08);fig.savefig(A/(pdb+'_backbone.png'),dpi=250,transparent=False,facecolor='white');fig.savefig(A/(pdb+'_backbone.svg'),transparent=False,facecolor='white');plt.close(fig)
 structures.append({'pdb_id':pdb,'gene':gene,'display_chain':chain,'display_CA_atoms':len(xyz),'rendering':'First model; longest deposited ATOM C-alpha chain; PCA orientation; cubic interpolation separately across <=5 Angstrom neighbors; white outlined backbone. Full deposited chain may contain engineered/fusion residues; not certified assay construct.'})

receipt={"scope":"Optional source-picture regeneration only; no new molecular/protein prediction or final graphical-abstract rendering.","rdkit":rdBase.rdkitVersion,"python":platform.python_version(),"structures":structures}
(A/"illustration_regeneration_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
