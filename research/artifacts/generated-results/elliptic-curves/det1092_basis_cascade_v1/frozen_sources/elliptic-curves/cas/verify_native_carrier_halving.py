#!/usr/bin/env python3
"""Sage-free complete rational-half classification from exact polynomial factors."""
import json
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import certify_compact_r17_candidates as cert
from alternate_quartic_covers import short_add,point_on_short_curve
from mod2_reduction_independence import _is_prime
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'native_carrier_halving_v1.json';OUT=ART/'native_carrier_halving_replay_v1.json'
def multiply(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def rational_roots(proof,expected):
    if list(map(F,proof['coefficients']))!=expected:raise ArithmeticError('duplication polynomial differs')
    product=[F(proof['unit'])];roots=[]
    for row in proof['factors']:
        f=list(map(F,row['coefficients']));m=row['multiplicity']
        if len(f)<2 or not f[-1] or type(m)!=int or m<1:raise ArithmeticError('invalid polynomial factor')
        for _ in range(m):product=multiply(product,f)
        if len(f)==2:roots.append(-f[0]/f[1])
        else:
            p=row['no_root_prime']
            if not _is_prime(p) or p>997 or any(c.denominator%p==0 for c in f) or f[-1].numerator%p==0:raise ArithmeticError('invalid no-root reduction')
            g=[c.numerator*pow(c.denominator,-1,p)%p for c in f]
            if any(sum(c*pow(x,i,p) for i,c in enumerate(g))%p==0 for x in range(p)):raise ArithmeticError('factor has a residue root')
    if product!=expected:raise ArithmeticError('complete exact factorization differs')
    return sorted(set(roots))

def main():
    d=cert.read(INPUT)
    if d['status']!='PASS' or any(cert.hashed(ROOT/n)!=h for n,h in d['sources'].items()):raise ArithmeticError('halving certificate binding differs')
    source=cert.read(ART/'native_rank3_carrier_marked_point_v1.json');model=list(map(F,d['curve']));A,B=model[3:]
    if d['curve']!=source['short_curve'] or d['basis_points']!=source['short_points'][:3] or source['rank_lower_bound']!=3:raise ArithmeticError('original independent auxiliary basis differs')
    basis=[tuple(map(F,P)) for P in d['basis_points']];T=tuple(map(F,d['torsion_point']));roots=rational_roots(d['torsion_cubic'],[B,A,F(0),F(1)])
    if roots!=[T[0]] or T[1]!=0 or not point_on_short_curve(model,T) or short_add(model,T,T) is not None:raise ArithmeticError('rational2-torsion classification differs')
    if [r['mask'] for r in d['rows']]!=list(range(1,16)):raise ArithmeticError('all15 nonzero parity classes required')
    for row in d['rows']:
        S=None
        for i,P in enumerate(basis+[T]):
            if row['mask']>>i&1:S=short_add(model,S,P)
        if S is None or list(map(str,S))!=row['sum_point']:raise ArithmeticError('exact parity sum differs')
        sx,sy=S;poly=[A*A-4*B*sx,-8*B-4*A*sx,-2*A,-4*sx,F(1)];xs=rational_roots(row['duplication'],poly);halves=[];lifts=[]
        for x in xs:
            v=x**3+A*x+B;n=isqrt(v.numerator) if v>=0 else -1;e=isqrt(v.denominator);square=v>=0 and n*n==v.numerator and e*e==v.denominator;ys=sorted({F(n,e),-F(n,e)}) if square else [];got=[]
            for y in ys:
                Q=(x,y)
                if short_add(model,Q,Q)==S:halves.append(Q);got.append(list(map(str,Q)))
            lifts.append({'x':str(x),'y_squared':str(v),'rational_square':square,'halves':got})
        if sorted(row['rational_x_lifts'],key=lambda r:F(r['x']))!=lifts or row['halves']!=[list(map(str,P)) for P in sorted(halves)]:raise ArithmeticError('complete rational lifts/doubling differ')
    masks=[r['mask'] for r in d['rows'] if r['halves']]
    if masks!=d['missing_coset_masks'] or d['two_saturated']!=(not masks):raise ArithmeticError('halving conclusion differs')
    result={'status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),INPUT,ROOT/'elliptic-curves/cas/alternate_quartic_covers.py']},'missing_coset_masks':masks,'two_saturated':not masks,'claim_boundary':d['claim_boundary']}
    if OUT.exists():
        if cert.read(OUT)!=result:raise ArithmeticError('retained independent replay differs')
    else:checkpoint(OUT,result)
    print('INDEPENDENT15 AUXILIARY HALVING CLASSES PASS; MISSING COSETS',masks,flush=True)
if __name__=='__main__':main()
