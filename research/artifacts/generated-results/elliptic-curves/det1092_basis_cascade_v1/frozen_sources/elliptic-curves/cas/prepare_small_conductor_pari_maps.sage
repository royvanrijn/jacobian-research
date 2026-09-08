#!/usr/bin/env sage-python
"""Freeze PARI-reduced horizontal maps from only the new rank26 curve generic points."""
import sys
from pathlib import Path
from math import lcm,isqrt
from sage.all import pari
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from half_lattice_pointed_sieve import linear_combination
from search_observability import transform,multiply
D=ROOT/'artifacts/local/elliptic-curves/small-conductor-adaptive-pari-maps-v1'
PARENT=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-small-conductor-followup-v2'

def run():
    protocol=cert.read(D/'protocol.json')
    if cert.hashed(Path(__file__).resolve())!=protocol['source_sha256']:raise ArithmeticError('map builder changed')
    source=cert.read(PARENT/'candidate-00/result.json');previous=cert.read(PARENT/'candidate-00/result.json')
    if cert.hashed(PARENT/'candidate-00/result.json')!=protocol['input_sha256'] or cert.hashed(PARENT/'candidate-00/result.json')!=protocol['centre_source_sha256']:raise ArithmeticError('generic map input changed')
    model=tuple(map(cert.F,source['curve']));points=tuple(tuple(map(cert.F,p)) for p in source['initial_state']['state']['reductions']['points']);rows=[]
    if (D/'maps.json').exists():raise FileExistsError('preserve map attempt')
    output={'status':'RUNNING','curve':source['curve'],'generic_points':source['generic_points'],'rows':rows}
    pari.allocatemem(256000000,silent=True)
    for c in previous['centres']:
        x,y=linear_combination(model,points,c['representative']);raw=(-3*x*x-4*model[3],-8*y,-6*x,cert.F(0),cert.F(1));den=lcm(*(v.denominator for v in raw))
        coefficients=[int(v*den*den) for v in raw];polynomial='+'.join(f'({v})*x^{i}' for i,v in enumerate(coefficients))
        result=pari('my(m1,m2,C1,C2);C1=hyperellminimalmodel(['+polynomial+',0],&m1);C2=hyperellred(C1,&m2);[C2,m1,m2]')
        first,second=result[1][1],result[2][1]
        matrices=[tuple(cert.F(str(m[i,j])) for i in range(2) for j in range(2)) for m in (first,second)];M=multiply(*matrices)
        P=[cert.F(str(result[0][0].polcoef(i))) for i in range(5)];Q=[cert.F(str(result[0][1].polcoef(i))) for i in range(3)]
        disc=[4*P[i]+sum(Q[j]*Q[i-j] for j in range(3) if 0<=i-j<3) for i in range(5)]
        transformed=transform(raw,M);k=next(i for i,v in enumerate(disc) if v);ratio=transformed[k]/disc[k]
        if ratio<=0 or any(a!=ratio*b for a,b in zip(transformed,disc)) or isqrt(ratio.numerator)**2!=ratio.numerator or isqrt(ratio.denominator)**2!=ratio.denominator:raise ArithmeticError('PARI horizontal quartic identity failed')
        rows.append({'centre':c,'raw_coefficients':list(map(str,raw)),'denominator_clearing':str(den),'first_matrix':list(map(str,matrices[0])),'second_matrix':list(map(str,matrices[1])),'matrix':list(map(str,M)),'reduced_P':list(map(str,P)),'reduced_Q':list(map(str,Q)),'discriminant_quartic':list(map(str,disc)),'square_ratio':str(ratio),'coordinate_policy':{'kind':'raw','matrix':list(map(str,M))}})
        checkpoint(D/'maps.json',output);print('FROZEN SMALL CONDUCTOR22 PARI MAP',len(rows),flush=True)
    output['status']='COMPLETE_DECLARED_MAPS';checkpoint(D/'maps.json',output)
if __name__=='__main__':run()
