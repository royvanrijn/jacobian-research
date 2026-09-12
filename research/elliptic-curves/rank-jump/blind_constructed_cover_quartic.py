#!/usr/bin/env python3
"""Verify conic/quartic transports and search a reduced binary quartic.

Uses only this blind arm's generated reductions and the frozen input indirectly
through their exactly verified degree-4 model. PARI hyperellratpoints searches
the abscissa height bound, stopping at one point. Infinity is checked exactly.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import re
import resource
import signal
import time

from blind_constructed_cover_check import OUT, ROOT, seq, dump, evaluate, primitive, det


def add(a,b):
    assert len(a)==len(b)
    return [x+y for x,y in zip(a,b)]


def scale(a,k):return [k*x for x in a]


def conv(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def power(a,n):
    out=[F(1)]
    for _ in range(n):out=conv(out,a)
    return out


def ternary(q,p):
    pairs=[(i,j) for i in range(3) for j in range(i,3)]
    out=[F(0)]*5
    for c,(i,j) in zip(q,pairs):out=add(out,scale(conv(p[i],p[j]),c))
    return out


def binary(q,T):
    aa=[T[0][0],T[1][0]];bb=[T[0][1],T[1][1]]
    out=[F(0)]*5
    for i,c in enumerate(q):out=add(out,scale(conv(power(aa,4-i),power(bb,i)),c))
    return out


def evalb(q,a,b):return sum(c*a**(len(q)-1-i)*b**i for i,c in enumerate(q))


def square_root(q):
    if q<0:return None
    a=isqrt(q.numerator);b=isqrt(q.denominator)
    return F(a,b) if a*a==q.numerator and b*b==q.denominator else None


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--column',type=int,choices=[6,7],required=True)
    ap.add_argument('--bound',type=int,choices=[1000,100000]);ns=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(60 if not ns.bound else 240)
    rows=[json.loads(x) for x in (OUT/'arithmetic-ledger.jsonl').read_text().splitlines()]
    used=sum(float(r['arithmetic_wall_seconds']) for r in rows)
    assert used+(240 if ns.bound else 60)<=1200
    start=time.monotonic();step=f'quartic-c{ns.column}-'+(f'b{ns.bound}' if ns.bound else 'check')
    try:
        c=json.loads((OUT/f'cover-worker-c{ns.column}-reduction.json').read_text())
        qs=[[F(x) for x in q] for q in c['reduced_quadrics']]
        assert qs[0][:4]==[0]*4 and qs[1][0]==1
        result=OUT/f'conic-quartic-minred-v3-c{ns.column}.result.txt'
        txt=result.read_text();assert 'VERIFIED_TRANSFORMATION' in txt and 'DONE' in txt
        pm=seq(txt,'CONIC_PARAMETRIZATION_MATRIX');assert len(pm)==9
        pm=[pm[3*i:3*i+3] for i in range(3)];assert det(pm)
        par=list(map(list,zip(*pm)))
        assert ternary(qs[0][4:],par)==[0]*5
        raw=seq(txt,'RAW_QUARTIC');reduced=seq(txt,'REDUCED_QUARTIC')
        assert len(raw)==len(reduced)==5
        lin=[sum(qs[1][i+1]*par[i][j] for i in range(3)) for j in range(3)]
        assert add(conv(lin,lin),scale(ternary(qs[1][4:],par),-4))==raw
        scalar=F(re.search(r'COMPOSED_SCALAR\s*([^\s]+)',txt).group(1))
        T=seq(txt,'COMPOSED_BINARY_MATRIX');assert len(T)==4
        T=[T[:2],T[2:]];assert det(T)
        assert scale(binary(raw,T),scalar*scalar)==reduced
        cert={'schema':'blind-cover-quartic.v1','column':ns.column,
              'degree4_reduction_certificate':f'cover-worker-c{ns.column}-reduction.json',
              'source_result':result.name,'source_sha256':hashlib.sha256(result.read_bytes()).hexdigest(),
              'conic_parametrization_matrix':pm,
              'conic_parametrization_convention':'[q1,q2,q3]=[a^2,a*b,b^2]*P',
              'raw_quartic_descending':raw,'reduced_quartic_descending':reduced,
              'degree2_scalar':scalar,'degree2_binary_matrix':T,
              'degree2_transformation':'[a,b]=[r,t]*T; reduced_F(r,t)=k^2*raw_F(a,b); raw_Y=reduced_Y/k',
              'degree4_recovery':'q0=(raw_Y-L(q1,q2,q3))/2; L is coefficients of q0*qj in the second reduced degree4 quadric',
              'linear_term_on_conic_param_descending':lin,
              'degree2_minimisation':'Magma Minimise with CrossTerms=false, then Reduce',
              'degree2_positive_level_primes':[2],
              'degree2_minimisation_boundary':'No-cross-term model may retain positive level at p=2; no local insolubility is inferred.',
              'exact_conic_parametrization_verified':True,'exact_quartic_discriminant_verified':True,
              'exact_binary_transformation_verified':True,'status':'REDUCED_QUARTIC_READY'}
        if ns.bound:
            searchstart=time.monotonic()
            yy=square_root(reduced[0])
            if yy is not None:
                rpt=[F(1),yy,F(0)];searchmethod='exact infinity square test'
            else:
                from cypari2 import Pari
                pari=Pari();pari.allocatemem(64000000,1073741824,silent=True)
                pp=pari('Polrev(['+','.join(str(q) for q in reversed(reduced))+'])')
                pts=pp.hyperellratpoints(ns.bound,1)
                searchmethod='PARI hyperellratpoints flag=1, after exact infinity test'
                rpt=None if len(pts)==0 else [F(str(pts[0][0])),F(str(pts[0][1])),F(1)]
                cert['pari_version']=str(pari('version()'))
            cert['search_wall_seconds']=time.monotonic()-searchstart
            cert['search_bound_abscissa_height']=ns.bound;cert['search_method']=searchmethod
            cert['status']='UNKNOWN_BOUNDED_MISS'
            if rpt is not None:
                r,Y,t=rpt;assert Y*Y==evalb(reduced,r,t)
                a=r*T[0][0]+t*T[1][0];b=r*T[0][1]+t*T[1][1]
                oldY=Y/scalar;assert oldY*oldY==evalb(raw,a,b)
                tail=[evalb(q,a,b) for q in par]
                L=sum(qs[1][i+1]*tail[i] for i in range(3))
                point=[(oldY-L)/2]+tail
                assert all(evaluate(q,point)==0 for q in qs)
                cert.update(status='RECOVERED_EXACT_ON_REDUCED_DEGREE4',
                            reduced_quartic_point=rpt,reduced_degree4_point=point,
                            reduced_degree4_primitive=[str(v) for v in primitive(point)],
                            reduced_binary_parameter=[r,t],raw_binary_parameter=[a,b])
                (OUT/f'quartic-c{ns.column}-b{ns.bound}.point.txt').write_text('RECOVERED_POINT\n[ '+', '.join(str(q) for q in point)+' ]\n')
        target=OUT/(step+'.json');target.write_text(json.dumps(dump(cert),indent=2)+'\n')
        print(json.dumps({'column':ns.column,'status':cert['status'],'certificate':str(target.relative_to(ROOT)),'point':dump(cert.get('reduced_degree4_point')),'search_wall_seconds':cert.get('search_wall_seconds')}),flush=True)
    finally:
        elapsed=time.monotonic()-start
        with (OUT/'arithmetic-ledger.jsonl').open('a') as h:
            h.write(json.dumps({'step':step,'kind':'search' if ns.bound else 'exact-verification',
                                'arithmetic_wall_seconds':elapsed,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})+'\n')
        print('Step wall seconds',elapsed,flush=True)


if __name__=='__main__':main()
