#!/usr/bin/env sage-python
"""One exact halving layer in the explicit rank3 auxiliary subgroup."""
import argparse,json,sys
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve,prime_range
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/native-carrier-halving-v1'
OUT=ART/'native_carrier_halving_v1.json';INPUT=ART/'native_rank3_carrier_marked_point_v1.json';REPLAY=ART/'native_rank3_carrier_marked_point_replay_v1.json'
def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),INPUT,REPLAY,CAS/'verify_native_carrier_halving.py',CAS/'alternate_quartic_covers.py',CAS/'research_runtime/store.py']}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve auxiliary halving layer')
    data=cert.read(INPUT);replay=cert.read(REPLAY)
    if data['status']!='PASS' or replay['status']!='PASS' or data['rank_lower_bound']!=3:raise ArithmeticError('explicit auxiliary rank3 proof required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.native-carrier-halving.v1','sources':sources(),'masks':list(range(1,16)),'maximum_halving_layers':1,'no_root_prime_bound':997,'seconds_per_stage':120,'rss_bytes':1610612736,'gate':'The marked carrier has exact auxiliary rank3, but its three explicit independent points have not been proved to generate a2-saturated subgroup. A coefficient cube in an unsaturated subgroup may omit smaller carrier points and hence new parameter images. Test this finite algebraic gap before another parameter-image cube.','scope':'On the short auxiliary model, determine the single rational point of order2 by exact cubic factorization. For all15 nonzero parity combinations of the three recorded independent points and that torsion point, factor the exact duplication quartic and classify every rational half. Nonlinear factors require a finite prime with no root; rational factors require exact square tests and doubling. One layer only, no auxiliary point search, adaptive saturation loop, original-fibre search, new rank or full Mordell-Weil basis claim. No halves proves2-saturation of this subgroup including rational2-torsion; any half records a missing coset and leaves further saturation open.'})
def factor_proof(poly):
    fac=poly.factor();rows=[]
    for f,m in fac:
        record={'coefficients':list(map(str,f.list())),'multiplicity':int(m)}
        if f.degree()>1:
            found=None
            for p in prime_range(3,998):
                if any(c.denominator()%p==0 for c in f.list()) or f.leading_coefficient()%p==0:continue
                coeffs=[int(c.numerator()%p)*pow(int(c.denominator()%p),-1,int(p))%int(p) for c in f.list()]
                if all(sum(c*pow(x,i,int(p)) for i,c in enumerate(coeffs))%p for x in range(int(p))):found=int(p);break
            if found is None:raise ArithmeticError('no bounded no-rational-root witness; keep UNKNOWN')
            record['no_root_prime']=found
        rows.append(record)
    return {'coefficients':list(map(str,poly.list())),'unit':str(fac.unit()),'factors':rows}
def expected():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen halving inputs differ')
    src=cert.read(INPUT);E=EllipticCurve(QQ,[QQ(v) for v in src['short_curve']]);A,B=E.a4(),E.a6();R=PolynomialRing(QQ,'x');x=R.gen();cubic=x**3+A*x+B;torsion_proof=factor_proof(cubic);roots=cubic.roots(multiplicities=False)
    if len(roots)!=1:raise ArithmeticError('one rational2-torsion generator required')
    T=E(roots[0],0);basis=[E([QQ(v) for v in P]) for P in src['short_points'][:3]];points=basis+[T];rows=[]
    for mask in p['masks']:
        S=sum((P for i,P in enumerate(points) if mask>>i&1),E(0))
        if not S:raise ArithmeticError('independent free basis/parity combination collapsed')
        sx,sy=S.xy();quartic=x**4-4*sx*x**3-2*A*x**2+(-8*B-4*A*sx)*x+A*A-4*B*sx
        proof=factor_proof(quartic);halves=[];lifts=[]
        for q in quartic.roots(multiplicities=False):
            value=q**3+A*q+B;square=value.is_square();ys=sorted(set([value.sqrt(),-value.sqrt()])) if square else []
            candidates=[]
            for y in ys:
                Q=E(q,y)
                if 2*Q==S:halves.append(Q);candidates.append(list(map(str,Q.xy())))
            lifts.append({'x':str(q),'y_squared':str(value),'rational_square':bool(square),'halves':candidates})
        rows.append({'mask':mask,'sum_point':list(map(str,S.xy())),'duplication':proof,'rational_x_lifts':lifts,'halves':[list(map(str,Q.xy())) for Q in sorted(halves,key=lambda Q:Q.xy())]})
        print('AUXILIARY HALVING',mask,'RATIONAL HALVES',len(halves),flush=True)
    return {'schema':'elliptic-curves.native-carrier-halving-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'curve':src['short_curve'],'basis_points':src['short_points'][:3],'torsion_point':list(map(str,T.xy())),'torsion_cubic':torsion_proof,'rows':rows,'missing_coset_masks':[r['mask'] for r in rows if r['halves']],'two_saturated':all(not r['halves'] for r in rows),'claim_boundary':p['scope']}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','build','check']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    else:
        if a.stage=='build' and OUT.exists():raise FileExistsError('preserve halving certificate')
        d=expected()
        if a.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('halving replay differs')
        else:checkpoint(OUT,d)
