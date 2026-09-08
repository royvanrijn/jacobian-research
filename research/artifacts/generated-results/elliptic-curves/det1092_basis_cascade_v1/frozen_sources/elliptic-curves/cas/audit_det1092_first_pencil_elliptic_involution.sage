#!/usr/bin/env sage-python
"""One exact degree15 invariant on the retained first-witness RR pencil.

Retrospective obstruction only. No point discovery or member selection.
The universal formula also supplies an equation for possible bielliptic loci.
"""
import hashlib
import json
from math import comb, factorial
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, gcd

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
INPUT = ART/'det1092_first_witness_pencil_genus_gate_v1.json'
OUT = ART/'det1092_first_pencil_elliptic_involution_v1.json'
LOCAL = ROOT/'artifacts/local/elliptic-curves/det1092-first-pencil-elliptic-involution-v1'
PRIMES = [17,47,53,61,67,71,79,83,89,101,107,113,127,137,149,179,191,197]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def emit_new(path, data):
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')


def falling(n,k):
    return 0 if n<k else factorial(n)//factorial(n-k)


def transvectant(f,g,k):
    """Unnormalized binary transvectant; lists ascend in the X exponent."""
    m,n=len(f)-1,len(g)-1
    out=[f[0]*0 for unused in range(m+n-2*k+1)]
    for i,fi in enumerate(f):
        for j,gj in enumerate(g):
            degree=i+j-k
            if degree<0 or degree>=len(out):
                continue
            weight=sum((-1)**s*comb(k,s)*falling(i,k-s)*falling(m-i,s)*falling(j,s)*falling(n-j,k-s) for s in range(k+1))
            if weight:
                out[degree]+=weight*fi*gj
    return out


def odd_invariant(f):
    i=transvectant(f,f,4)
    y1=transvectant(f,i,4)
    y2=transvectant(i,y1,2)
    y3=transvectant(i,y2,2)
    a,b,c=y1,y2,y3
    value=a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
    return value,[y1,y2,y3]


def run():
    if OUT.exists() or LOCAL.exists():
        raise FileExistsError('Preserve the completed or interrupted experiment')
    LOCAL.mkdir(parents=True)
    protocol={'classification':'retrospective exact diagnostic, not blind discovery',
              'input':str(INPUT.relative_to(ROOT)), 'input_sha256':sha(INPUT),
              'source_sha256':sha(Path(__file__)),
              'pencil':'Exactly the immutable pencil u=u0 through the historical first witness; no substitutions.',
              'invariant':'i=(F,F)_4; Y1=(F,i)_4; Y2=(i,Y1)_2; Y3=(i,Y2)_2; I15=det(coefficients(Y1,Y2,Y3)). Unnormalized transvectants.',
              'primes':PRIMES,
              'stop':'First same-degree no-root reduction of the exact primitive invariant, or all fixed primes exhausted. Keep all trials.',
              'limits':{'diagnostic_pencils':1,'exact_invariants':1,'prime_cap':18,'prime_bound':197,
                        'point_searches':0,'Selmer_runs':0,'class_group_runs':0}}
    emit_new(LOCAL/'protocol.json',protocol)
    d=json.loads(INPUT.read_text())
    V=PolynomialRing(ZZ,'v')
    f=[V(row) for row in d['q_coefficients_t_then_v']]
    assert len(f)==7
    invariant,covariants=odd_invariant(f)
    assert invariant.degree()<=60
    print('Exact invariant computed; degree',invariant.degree(),flush=True)
    if invariant:
        content=gcd(invariant.list())
        if invariant.leading_coefficient()<0:
            content=-content
        primitive=V(invariant/content)
    else:
        content=ZZ(0)
        primitive=V(0)
    saved={'primitive_coefficients':list(map(str,primitive.list())),
           'content':str(content),'degree':int(primitive.degree()),
           'covariant_coefficients':[[list(map(str,p.list())) for p in row] for row in covariants]}
    emit_new(LOCAL/'exact-invariant.json',saved)
    trials=[]
    found=False
    if invariant:
        for p in PRIMES:
            S=PolynomialRing(GF(p),'v')
            fp=S(primitive.list())
            roots=[int(a) for a in GF(p) if fp(a)==0]
            row={'prime':p,'degree':int(fp.degree()),'coefficients':list(map(int,fp.list())),
                 'roots':roots,'same_degree':bool(fp.degree()==primitive.degree())}
            trials.append(row)
            emit_new(LOCAL/('prime-%03d.json'%p),row)
            print('prime',p,'degree',fp.degree(),'roots',roots,flush=True)
            if row['same_degree'] and not roots:
                found=True
                break
    result={'classification':'retrospective verified application and new deduction',
            'status':'PASS_NO_RATIONAL_BIELLIPTIC_MEMBER_IN_FIRST_WITNESS_PENCIL' if found else 'INCONCLUSIVE_FIXED_PRIME_GATE',
            'protocol':protocol,'protocol_sha256':sha(LOCAL/'protocol.json'),
            'witness_coordinate':d['witness_coordinate'],'u0':d['u0'],
            'exact_invariant':saved,'modular_trials':trials,
            'no_rational_invariant_zero':found,
            'universal_necessary_condition':'I15(q(T;u,v))=0 is necessary for a smooth RR member to admit a degree2 map to an elliptic curve over the algebraic closure.',
            'boundary':'Only the fixed first-witness pencil is tested for rational invariant zeros. Other RR pencils and elliptic quotients of degree at least3 remain open. No Selmer or rational-solubility conclusion follows, and no execution policy uses the witness.'}
    emit_new(OUT,result)
    print(result['status'],flush=True)


if __name__=='__main__':
    run()
