#!/usr/bin/env python3
"""Exact, standard-library-only checking of the blind cover reduction/recovery.

The only mathematical source input is the frozen covers.json. No point data is
read from elsewhere. The remaining inputs are this arm's logged Magma exports.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import gcd, lcm
from pathlib import Path
import re
import resource
import signal
import time

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'research/artifacts/local/elliptic-curves/constructed-class-blind-cover-v1'
INPUT = ROOT / 'research/artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/covers.json'
SHA = 'a8218511b1ac1b2def5dafe9c7803ebbdb0d33a7a0ffffca1c8e33a090efca1f'
PAIRS = [(i,j) for i in range(4) for j in range(i,4)]


def dump(obj):
    if isinstance(obj, F): return str(obj)
    if isinstance(obj, dict): return {k:dump(v) for k,v in obj.items()}
    if isinstance(obj, (list,tuple)): return [dump(v) for v in obj]
    return obj


def seq(text, marker):
    match = re.search(re.escape(marker) + r'\s*\[([^\]]*)\]', text)
    assert match, marker
    return [F(s.strip()) for s in match.group(1).split(',')]


def det(matrix):
    a = [[F(v) for v in row] for row in matrix]
    n = len(a); d=F(1)
    for i in range(n):
        if not a[i][i]:
            k=next((k for k in range(i+1,n) if a[k][i]),None)
            if k is None: return F(0)
            a[i],a[k]=a[k],a[i];d=-d
        d*=a[i][i]
        for j in range(i+1,n):
            r=a[j][i]/a[i][i]
            for k in range(i+1,n): a[j][k]-=r*a[i][k]
    return d


def mul(a,b,f):
    t=[F(0)]*5
    for i in range(3):
        for j in range(3): t[i+j]+=a[i]*b[j]
    for k in (4,3):
        for i in range(3): t[k-3+i]-=t[k]*f[i]
    return t[:3]


def norm(a,f):
    cols=[mul(a,[F(int(i==j)) for i in range(3)],f) for j in range(3)]
    return det(list(zip(*cols)))


def parse_quad(s):
    out=[F(0)]*10
    for term in re.findall(r'[+-]?[^+-]+',s.replace(' ','')):
        sign=-1 if term.startswith('-') else 1
        term=term.lstrip('+-')
        parts=term.split('*');c=F(sign)
        if parts[0][0].isdigit(): c*=F(parts.pop(0))
        inds=[]
        for p in parts:
            if '^' in p:
                v,e=p.split('^');inds += ['uvws'.index(v)]*int(e)
            else: inds.append('uvws'.index(p))
        assert len(inds)==2
        out[PAIRS.index(tuple(sorted(inds)))]+=c
    return out


def raw_quadrics(case,f):
    beta=[F(x) for x in case['beta_ascending']]
    assert norm(beta,f)==F(case['norm'])
    assert F(case['positive_norm_square_root'])**2==F(case['norm'])
    qs=[[F(0)]*10 for _ in range(3)]
    for n,(i,j) in enumerate(PAIRS):
        if j==3: continue
        a=[F(int(k==i)) for k in range(3)]
        b=[F(int(k==j)) for k in range(3)]
        z=mul(beta,mul(a,b,f),f)
        for k in range(3): qs[k][n]=z[k]*(1 if i==j else 2)
    source=[qs[2][:],qs[1][:]];source[1][-1]+=1
    assert source==[parse_quad(q) for q in case['cover_quadrics']]
    assert qs[0]==parse_quad(case['x_numerator'])
    return beta,qs,source


def substitute(q,S):
    # Original row coordinates = reduced row coordinates times S.
    return [sum(c*(S[k][i]*S[l][j] + (S[l][i]*S[k][j] if k!=l else 0))
                for c,(i,j) in zip(q,PAIRS)) for k,l in PAIRS]


def evaluate(q,p):
    return sum(c*p[i]*p[j] for c,(i,j) in zip(q,PAIRS))


def primitive(p):
    d=lcm(*(v.denominator for v in p))
    pp=[int(v*d) for v in p]
    g=gcd(*pp);pp=[v//g for v in pp]
    if next(v for v in pp if v)<0: pp=[-v for v in pp]
    return pp


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--column',type=int,choices=[6,7],required=True)
    ap.add_argument('--point-result');ns=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    signal.alarm(60)
    start=time.monotonic()
    step=('check-recovery-' if ns.point_result else 'check-reduction-')+str(ns.column)
    try:
        assert hashlib.sha256(INPUT.read_bytes()).hexdigest()==SHA
        data=json.loads(INPUT.read_text())
        case=next(c for c in data['cases'] if c['column']==ns.column)
        f=[F(x) for x in data['cubic_ascending']]
        beta,qs,source=raw_quadrics(case,f)
        result=OUT/f'compact-minred-c{ns.column}.result.txt'
        text=result.read_text()
        assert 'VERIFIED_TRANSFORMATION' in text and 'DONE' in text
        model=seq(text,'REDUCED_MODEL');assert len(model)==20
        A=seq(text,'COMPOSED_EQUATION_TRANS');assert len(A)==4
        A=[A[:2],A[2:]]
        S=seq(text,'COMPOSED_COORDINATE_TRANS');assert len(S)==16
        S=[S[4*i:4*i+4] for i in range(4)]
        assert det(A) and det(S)
        model=[model[:10],model[10:]]
        sub=[substitute(q,S) for q in source]
        transformed=[[sum(A[i][k]*sub[k][j] for k in range(2)) for j in range(10)] for i in range(2)]
        assert transformed==model, 'Magma tuple orientation/model coefficients do not verify'
        cert={'schema':'blind-cover-certificate.v1','column':ns.column,
              'input_path':str(INPUT.relative_to(ROOT)),'input_sha256':SHA,
              'backend':'Magma V2.29-10 standard public calculator',
              'algorithm':'Cremona-Fisher-Stoll genus-one degree-4 Minimise then Reduce',
              'minimisation_positive_level_primes':[],
              'source_monomial_order':PAIRS,'reduced_quadrics':model,
              'equation_matrix_A':A,'coordinate_matrix_S':S,
              'transformation_convention':'old_row = reduced_row * S; reduced_quadric_column = A * original_quadric_column(old_row)',
              'equation_matrix_determinant':det(A),'coordinate_matrix_determinant':det(S),
              'exact_polynomial_transport_verified':True,
              'source_beta_expansion_verified':True,'norm_invariant_verified':True,
              'raw_max_coefficient_bits':case['max_coefficient_bits'],
              'reduced_max_coefficient_bits':max(abs(c.numerator).bit_length() for q in model for c in q),
              'status':'REDUCED_SEARCH_NOT_YET_RUN','inputs':[str(INPUT.relative_to(ROOT)),str(result.relative_to(ROOT))]}
        if ns.point_result:
            resultp=Path(ns.point_result).resolve();assert resultp.parent==OUT
            pt=seq(resultp.read_text(),'RECOVERED_POINT');assert len(pt)==4
            assert all(evaluate(q,pt)==0 for q in model)
            old=[sum(pt[i]*S[i][j] for i in range(4)) for j in range(4)]
            op=primitive(old)
            assert all(evaluate(q,op)==0 for q in source)
            u,v,w,s=map(F,op);assert s
            bval=mul(beta,mul([u,v,w],[u,v,w],f),f)
            assert bval[2]==0 and bval[1]==-s*s
            X=bval[0]/s**2
            Y=F(case['positive_norm_square_root'])*norm([u,v,w],f)/s**3
            x=X/4;y=(Y-X)/8
            a1,a2,a3,a4,a6=map(F,data['original_curve'])
            assert y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6
            cert.update(status='RECOVERED_EXACT',reduced_point=pt,
                        reduced_point_primitive=[str(x) for x in primitive(pt)],
                        original_point_primitive=[str(x) for x in op],
                        original_point_scale_from_raw_transport=str(F(op[3],old[3])),
                        curve_X=X,curve_Y=Y,curve_x=x,curve_y=y,
                        normalised_xi_coefficients=[u/s,v/s,w/s],
                        xi_identity='beta*xi^2 = 4*x - theta',
                        exact_original_cover_point_verified=True,exact_curve_map_verified=True,
                        point_result_sha256=hashlib.sha256(resultp.read_bytes()).hexdigest())
            cert['inputs'].append(str(resultp.relative_to(ROOT)))
        suffix='recovery' if ns.point_result else 'reduction'
        target=OUT/f'cover-worker-c{ns.column}-{suffix}.json'
        target.write_text(json.dumps(dump(cert),indent=2)+'\n')
        print(json.dumps({'column':ns.column,'status':cert['status'],'reduced_max_coefficient_bits':cert['reduced_max_coefficient_bits'],'certificate':str(target.relative_to(ROOT))}))
    finally:
        elapsed=time.monotonic()-start
        with (OUT/'arithmetic-ledger.jsonl').open('a') as h:
            h.write(json.dumps({'step':step,'kind':'exact-verification','arithmetic_wall_seconds':elapsed,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})+'\n')
        print('Exact checking wall seconds',elapsed,flush=True)


if __name__=='__main__':main()
