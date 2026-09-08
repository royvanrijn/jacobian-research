#!/usr/bin/env sage-python
"""Bounded exact conductor completion for the new3/17 rank22 curve."""
import argparse
from pathlib import Path
import sys
from sage.all import QQ, ZZ, EllipticCurve, prod
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ART=ROOT/'artifacts/generated-results/elliptic-curves'
BOUNDS=ROOT/'artifacts/local/elliptic-curves/next12-new-curve-conductor-bounds-v1/result.json'
POINTS=ART/'prospective_mw16_next12_results_v1.json'

def build(output):
    if output.exists():raise FileExistsError('preserve conductor certificate')
    bound=next(r for r in cert.read(BOUNDS)['rows'] if r['parameter']=='3/17')
    row=next(r for r in cert.read(POINTS)['curves'] if r['parameter']=='3/17')
    if row['rank_lower_bound']!=22 or row['curve']!=bound['curve']:raise ArithmeticError('rank22 input differs')
    residual=ZZ(bound['unprocessed_cofactor'])
    if residual.nbits()!=220:raise ArithmeticError('fixed220-bit cofactor gate differs')
    result={'schema':'elliptic-curves.next12-rank22-exact-conductor.v1','status':'RUNNING',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),BOUNDS,POINTS)},
        'parameter':'3/17','short_model':row['curve'],'rank_lower_bound':22,'residual':str(residual)}
    checkpoint(output,result)
    factors=list(residual.factor(proof=True))
    if prod(p**e for p,e in factors)!=residual or any(not p.is_prime(proof=True) for p,e in factors):raise ArithmeticError('residual prime factorization failed')
    result['residual_factorization']=[[str(p),int(e)] for p,e in factors];checkpoint(output,result)
    print('PROVED RESIDUAL FACTORIZATION',result['residual_factorization'],flush=True)
    integral=[1,0,0,-182451976602578656424609725499140,710003150253794219215652666162794189038512805392]
    E=EllipticCurve(QQ,integral);D=ZZ(E.discriminant());short=EllipticCurve(QQ,[QQ(q) for q in row['curve']])
    if E.c_invariants()!=short.c_invariants() or str(D)!=bound['displayed_discriminant']:raise ArithmeticError('integral model transport differs')
    all_factors=[(ZZ(r['prime']),ZZ(r['displayed_discriminant_valuation'])) for r in bound['local_data'] if r['displayed_discriminant_valuation']]+factors
    if prod(p**e for p,e in all_factors)!=abs(D):raise ArithmeticError('full discriminant factorization differs')
    local=[];N=ZZ(1)
    for p,e in all_factors:
        data=E.local_data(p,proof=True);f=int(data.conductor_valuation());N*=p**f
        minimal_val=int(data.minimal_model().discriminant().valuation(p))
        local.append({'prime':str(p),'discriminant_valuation':int(e),'minimal_discriminant_valuation':minimal_val,
            'conductor_valuation':f,'kodaira_symbol':str(data.kodaira_symbol())})
    result.update(status='PASS_EXACT_CONDUCTOR',integral_model=list(map(str,integral)),
        point_transport={'x':'X-1/12','y':'Y-x/2'},
        discriminant=str(D),discriminant_factorization=[[str(p),int(e)] for p,e in all_factors],
        local_data=local,conductor=str(N),integral_model_is_global_minimal=all(r['discriminant_valuation']==r['minimal_discriminant_valuation'] for r in local),
        claim_boundary='Exact conductor from a proven complete prime factorization and exact local Tate data; rank lower bound is supplied by the separately replayed rational point proof. No exact rank or universal novelty claim.')
    checkpoint(output,result);print('EXACT CONDUCTOR',N,'global minimal',result['integral_model_is_global_minimal'],flush=True)

def check(path):
    data=cert.read(path)
    if data['status']!='PASS_EXACT_CONDUCTOR':raise ArithmeticError('conductor unresolved')
    for p,h in data['sources'].items():
        if cert.hashed(ROOT/p)!=h:raise ArithmeticError('conductor source changed')
    E=EllipticCurve(QQ,[QQ(q) for q in data['integral_model']]);D=ZZ(E.discriminant())
    factors=[(ZZ(p),int(e)) for p,e in data['discriminant_factorization']]
    if any(not p.is_prime(proof=True) for p,e in factors) or prod(p**e for p,e in factors)!=abs(D):raise ArithmeticError('complete prime factorization failed')
    local=[];N=ZZ(1)
    for p,e in factors:
        row=E.local_data(p,proof=True);f=int(row.conductor_valuation());N*=p**f
        local.append({'prime':str(p),'discriminant_valuation':e,'minimal_discriminant_valuation':int(row.minimal_model().discriminant().valuation(p)),
            'conductor_valuation':f,'kodaira_symbol':str(row.kodaira_symbol())})
    if local!=data['local_data'] or str(N)!=data['conductor'] or all(r['discriminant_valuation']==r['minimal_discriminant_valuation'] for r in local)!=data['integral_model_is_global_minimal']:raise ArithmeticError('exact conductor replay differs')
    print('REPLAYED EXACT CONDUCTOR',N,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--check',type=Path);a=p.parse_args()
    check(a.check.resolve()) if a.check else build(a.output.resolve())
