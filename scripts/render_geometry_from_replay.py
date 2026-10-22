"""Render geometry from freshly regenerated tables using unchanged figure code."""
from pathlib import Path
import argparse,hashlib,json,re,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,default=ROOT/'outputs/geometry/data');p.add_argument('--output-dir',type=Path,default=ROOT/'outputs/geometry_figures');args=p.parse_args()
    source=args.input_dir.resolve();output=args.output_dir.resolve();(output/'data').mkdir(parents=True,exist_ok=True)
    script=ROOT/'figures/build_geometry_figures.py'
    names=sorted(set(re.findall(r"D/'([^']+\.csv)'",script.read_text())))
    files=[]
    for name in names:
        path=source/name
        if not path.is_file():raise FileNotFoundError('Fresh geometry table required: '+name)
        target=output/'data'/name;shutil.copyfile(path,target)
        if sha(path)!=sha(target):raise AssertionError('Geometry copy mismatch: '+name)
        files.append(dict(table=name,sha256=sha(path)))
    for name in ['build_geometry_figures.py','orx_figstyle.py']:
        shutil.copyfile(ROOT/'figures'/name,output/name)
    subprocess.run([sys.executable,'build_geometry_figures.py'],cwd=output,check=True)
    receipt=dict(scope='Fresh scientific-array geometry tables rendered with unchanged published figure/style sources; no cached CSV fallback, fitting or network.',input_tables=files,renderer_sha256=sha(script),style_sha256=sha(ROOT/'figures/orx_figstyle.py'),outputs={p.name:sha(p) for p in output.glob('*.png')})
    (output/'fresh_geometry_render_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(rendered_geometry_figures=len(receipt['outputs']),input_tables=len(names))))


if __name__=='__main__':main()
