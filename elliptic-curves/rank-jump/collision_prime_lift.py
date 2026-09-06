#!/usr/bin/env python3
"""Fixed collision support and exact real sign geometry for three quartets."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'COLLISION_PRIME_LIFT_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_soluble_quartet_compression_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_collision_prime_lift_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-collision-prime-lift-v1'


def worker(index):
    from sage.all import QQ,ZZ,AA,PolynomialRing,prod
    case=r.read(INPUT)['cases'][index];R=PolynomialRing(QQ,'t');qs=[R(c['form']) for c in case['covers']]
    assert all(q.degree()==2 and q.is_squarefree() for q in qs)
    wd=WORK/str(index);wd.mkdir(parents=True,exist_ok=True);resultants=[];support=set()
    for i in range(4):
        for j in range(i):
            value=ZZ(qs[i].resultant(qs[j]));assert value
            path=wd/f'resultant_{j}_{i}.json'
            if path.exists():row=r.read(path)
            else:
                factor=list(abs(value).factor());assert prod(p**e for p,e in factor)==abs(value)
                assert all(p.is_prime(proof=True) for p,e in factor)
                row={'indices':[j,i],'resultant':str(value),'prime_factorization':[[str(p),int(e)] for p,e in factor]}
                r.write_new(path,row)
            support.update(int(p) for p,e in row['prime_factorization']);resultants.append(row)
    # Order exact algebraic roots and sweep signs. No rational t is sampled.
    roots=[]
    for i,q in enumerate(qs):
        rr=sorted(q.roots(AA,multiplicities=False))
        roots.extend((a,i,j) for j,a in enumerate(rr))
    roots.sort(key=lambda row:row[0]);assert all(roots[i][0]!=roots[i-1][0] for i in range(1,len(roots)))
    signs=[1 if q[2]>0 else -1 for q in qs];cells=[];branch=[]
    def cell():return {'signs':list(signs),'product_positive':prod(signs)>0,'all_native_positive':all(s>0 for s in signs)}
    cells.append(cell())
    for a,i,j in roots:
        branch.append({'cover_index':i,'root_index':j,'root_approximation':str(a),
            'other_native_values_positive':all(signs[k]>0 for k in range(4) if k!=i)})
        signs[i]*=-1;cells.append(cell())
    real_surjective=all(not c['product_positive'] or c['all_native_positive'] for c in cells)
    real_surjective=real_surjective and all(b['other_native_values_positive'] for b in branch)
    infprod=prod(q[2] for q in qs)>0;inf_native=all(q[2]>0 for q in qs)
    real_surjective=real_surjective and (not infprod or inf_native)
    row={'id':case['id'],'pair_resultants':resultants,'collision_primes':[str(p) for p in sorted(support)],
        'real_branch_order':branch,'real_open_cells':cells,'real_product_points_at_infinity':bool(infprod),
        'real_native_lifts_at_infinity':bool(inf_native),'real_lift_map_surjective':bool(real_surjective)}
    if case['observed_parameter'] is not None:
        t=QQ(case['observed_parameter']);vals=[q(t) for q in qs]
        assert all(v>0 and v.is_square() for v in vals)
        parities=[[int(v.valuation(p)%2) for p in sorted(support)] for v in vals]
        assert all(not any(v) for v in parities)
        row['retained_parity_check']={'parameter':str(t),'all_values_positive':True,'all_collision_valuations_even':True}
    r.write_new(wd/'result.json',row)
    print(case['id'],'collision primes',len(support),'real map onto',real_surjective,flush=True)


def run():
    rows=[]
    for i in range(3):
        wd=WORK/str(i);wd.mkdir(parents=True,exist_ok=True);log=wd/'worker.log';ex=wd/'execution.json'
        if not log.exists():
            with log.open('x') as out:
                try:
                    p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker','--case',str(i)],stdout=out,stderr=out,timeout=60)
                    status={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
                except subprocess.TimeoutExpired:status={'status':'TIMEOUT'}
            r.write_new(ex,status)
        path=wd/'result.json'
        rows.append({'index':i,'execution':r.read(ex),'log':log.read_text(),'result':r.read(path) if path.exists() else {'status':'UNKNOWN'}})
        print(rows[-1]['execution'],rows[-1]['log'],flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.collision-prime-lift.v1','layer':'solubility','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['boundary']})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['run','worker']);p.add_argument('--case',type=int);a=p.parse_args()
    if a.mode=='run':run()
    else:worker(a.case)
