"""Replay local complementarity and a two-model coverage bound on one chart."""
from fractions import Fraction as F
from itertools import combinations, product
import hashlib
import json
import math
from pathlib import Path
import time

from verify_pointed_height_bounds import verify, require, homogeneous
from verify_pointed_minimal_neighbours import verify as verify_neighbours, matrix_mul


def run(folder,root):
    folder,root=Path(folder).resolve(),Path(root).resolve()
    files=[folder/'preconditioned_full-bounds.json']+list(sorted((folder/'neighbours').glob('neighbour-*-bounds.json')))
    require(len(files)==3,'frozen three-model set incomplete')
    packets=[json.loads(p.read_text()) for p in files]
    for p in packets:verify(p,root)
    verify_neighbours(folder,root)
    base=packets[0];M=list(map(F,base['input']['mapping']['matrix']))
    a,b,c,d=M;det=a*d-b*c;inv=[d/det,-b/det,-c/det,a/det]
    common=int(base['finite_gcd_divisor']);partitions=[]
    for data in packets[1:]:
        witness=data['input']['neighbour_witness']
        p,r=int(witness['prime']),int(witness['abscissa_residue'])
        red=list(map(F,witness['reduction_matrix']));a,b,c,d=red
        require(all(v.denominator==1 for v in red) and abs(a*d-b*c)==1,'reduction is not GL2(Z)')
        T=matrix_mul(inv,list(map(F,data['input']['mapping']['matrix'])));a,b,c,d=T
        require(abs(a*d-b*c)==p,'wrong neighbour index')
        adj=[d,-b,-c,a]
        for name in ('numerator','denominator'):
            transformed=homogeneous(list(map(F,data[name])),adj)
            original=list(map(F,base[name]))
            require(transformed==[p*p*v for v in original],'local covariance factor differs')
        row=next(v for v in base['local_trees'] if int(v['prime'])==p)
        node=row['affine']
        require(node['kind']=='split' and node['content_valuation']==0,'unavailable first residue partition')
        inside=next(c['node']['upper'] for c in node['children'] if c['residue']==r)
        outside=max([0,row['infinity']['upper']]+[c['node']['upper'] for c in node['children'] if c['residue']!=r])
        require(inside>=2 and outside>=0,'unsupported empty phase')
        partitions.append({'prime':p,'residue':r,'base_cancellation_caps':[outside,inside],
            'new_cancellation':'k_base + 2 - 4*inside'})
        common//=p**row['upper']
    cells=[]
    for phase in product((0,1),repeat=2):
        C=common
        for i,r in enumerate(partitions):C*=r['prime']**r['base_cancellation_caps'][phase[i]]
        bounds=[]
        for i,data in enumerate(packets):
            correction=F(1) if i==0 else F(partitions[i-1]['prime'])**(2-4*phase[i-1])
            bounds.append(F(C)*correction/F(data['real']['lower']))
        cells.append({'inside_residue_balls':list(phase),'D_by_model':list(map(str,bounds))})
    choices=[]
    for size in range(1,4):
        for indices in combinations(range(3),size):
            D=max(min(F(row['D_by_model'][i]) for i in indices) for row in cells)
            choices.append({'models':list(indices),'D':str(D),'squared_address_work_proxy':str(size*size*D)})
    best=min(choices,key=lambda r:(F(r['squared_address_work_proxy']),len(r['models']),r['models']))
    single=min((r for r in choices if len(r['models'])==1),key=lambda r:F(r['D']))
    factor=(float(F(single['D'])/F(best['D'])))**0.25
    return {'status':'PASS_INDEPENDENT_LOCAL_PORTFOLIO_REPLAY',
        'packet_bindings':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        'model_order':[p['input']['policy'] for p in packets], 'partitions':partitions,
        'joint_cells':cells,'subsets':choices,'selected_by_address_work_proxy':best,
        'best_single':single,'display':{'portfolio_B2':math.log(float(F(best['D'])))/4,
            'parameter_height_factor_over_best_single':factor,
            'address_work_ratio_to_best_single':len(best['models'])/factor**2},
        'theorem':'For each rational point P on the marked chart, at least one selected model has '
            'H(parameter)^4 <= D * H_x(2P-Q). A box of height H on each selected model therefore '
            'covers all such points with H_x(2P-Q)<=H^4/D, including projective parameters.',
        'boundary':'The same model must be selected across all places for a point. The max-over-cells '
            'of min-over-models respects that coupling. The model-count times H^2 proxy is not CPU '
            'or a rank-gain prediction. No new point search ran, no target point selected a model, '
            'and neither exhaustive model enumeration nor cross-curve transfer is claimed.'}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();require(not args.output.exists(),'preserve previous portfolio receipt')
    start=time.process_time();result=run(args.folder,Path(__file__).resolve().parents[2])
    result.update(cpu_seconds=time.process_time()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(result['status'],result['display'])
