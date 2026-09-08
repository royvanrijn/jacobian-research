#!/usr/bin/env sage-python
"""Replay with Sage's documented Mestre covariants and exact root gcd.

No constructor import. One immutable diagnostic pencil, no point searches.
"""
import hashlib
import inspect
import json
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, matrix, gcd, power_mod
from sage.schemes.hyperelliptic_curves.invariants import ubs

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
CERT=ART/'det1092_first_pencil_elliptic_involution_v1.json'
OUT=ART/'det1092_first_pencil_elliptic_involution_replay_v1.json'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def det3(rows):
    a,b,c=rows
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])


def covariant_rows(f):
    data=ubs(f)
    X,Z=f.parent().gens()
    monomials=[Z*Z,X*Z,X*X]
    return [[data[key].monomial_coefficient(m) for m in monomials] for key in ['y1','y2','y3']]


def verify():
    d=read(CERT)
    inp=ROOT/d['protocol']['input']
    assert sha(inp)==d['protocol']['input_sha256']
    constructor=ROOT/'elliptic-curves/cas/audit_det1092_first_pencil_elliptic_involution.sage'
    assert sha(constructor)==d['protocol']['source_sha256']
    pencil=read(inp)
    up=ART/'det1092_universal_rr_descent_preflight_v1.json'
    bp=ART/'det1092_first_centre_rr_net_replay_v1.json'
    universal=read(up)['universal_genus2']
    bridge=read(bp)
    V=PolynomialRing(QQ,'v')
    v=V.gen()
    S=PolynomialRing(V,['X','Z'])
    X,Z=S.gens()
    coeff=list(map(V,pencil['q_coefficients_t_then_v']))
    F=sum(coeff[i]*X**i*Z**(6-i) for i in range(7))
    source_R=PolynomialRing(QQ,'z')
    ur=bridge['chart_to_net_u']
    z0=QQ(bridge['witness_coordinate'])
    u0=source_R(ur['numerator'])(z0)/source_R(ur['denominator'])(z0)
    assert u0==QQ(pencil['u0'])==QQ(d['u0'])
    assert str(z0)==d['witness_coordinate']
    expected=[V(0) for i in range(7)]
    for row in universal['sparse_q']:
        expected[row['T']]+=QQ(row['coefficient'])*u0**row['u']*v**row['v']
    assert [QQ(universal['scale'])*p for p in expected]==[QQ(pencil['branch_squareclass_scale'])*p for p in coeff]
    # Documented normalized transvectants, independently implemented by Sage.
    normalized=covariant_rows(F)
    factors=[ZZ(129600)*8640]
    factors.append(ZZ(129600)*factors[0]*24)
    factors.append(ZZ(129600)*factors[1]*24)
    exact=d['exact_invariant']
    stored=[[V(c) for c in row] for row in exact['covariant_coefficients']]
    for i in range(3):
        assert [factors[i]*x for x in normalized[i]]==stored[i]
    raw=V(det3(stored))
    primitive=V(exact['primitive_coefficients'])
    content=ZZ(exact['content'])
    assert raw==content*primitive
    assert all(c.denominator()==1 for c in primitive.list())
    assert gcd([ZZ(c) for c in primitive.list()])==1
    assert primitive.degree()==exact['degree']==45 and primitive.leading_coefficient()>0
    assert len(d['modular_trials'])<=d['protocol']['limits']['prime_cap']
    for index,row in enumerate(d['modular_trials']):
        p=ZZ(row['prime'])
        assert p==d['protocol']['primes'][index] and p.is_prime()
        P=PolynomialRing(GF(p),'v')
        x=P.gen()
        f=P(primitive.list())
        assert list(map(int,f.list()))==row['coefficients']
        assert f.degree()==row['degree']
        assert row['same_degree']==(f.degree()==primitive.degree())
        rootpoly=f.gcd(power_mod(x,p,f)-x)
        assert rootpoly.degree()==len(row['roots'])
        assert all(f(a)==0 for a in row['roots']) and len(set(row['roots']))==len(row['roots'])
    last=d['modular_trials'][-1]
    assert last['prime']==127 and last['degree']==45 and not last['roots']
    assert d['no_rational_invariant_zero']
    # Formula regressions: even sextics vanish symbolically; arbitrary small
    # test forms transform with determinant weight45, including det=-1.
    K=PolynomialRing(QQ,['a0','a2','a4','a6'])
    R=PolynomialRing(K,['X','Z'])
    X,Z=R.gens()
    a0,a2,a4,a6=K.gens()
    even=a0*Z**6+a2*X**2*Z**4+a4*X**4*Z**2+a6*X**6
    assert det3(covariant_rows(even))==0
    R=PolynomialRing(QQ,['X','Z'])
    X,Z=R.gens()
    test=X**6+2*X**5*Z+3*X**4*Z**2+5*X**3*Z**3+7*X**2*Z**4+11*X*Z**5+13*Z**6
    invariant=det3(covariant_rows(test))
    assert invariant
    matrices=[[[1,1],[0,1]],[[0,1],[1,0]],[[2,1],[1,1]],[[2,0],[0,1]]]
    for rows in matrices:
        M=matrix(QQ,rows)
        transformed=test(M[0,0]*X+M[0,1]*Z,M[1,0]*X+M[1,1]*Z)
        assert det3(covariant_rows(transformed))==M.det()**45*invariant
    sage_source=Path(inspect.getsourcefile(ubs))
    return {'classification':'retrospective verified application and new deduction',
            'status':'PASS_INDEPENDENT_NO_RATIONAL_BIELLIPTIC_MEMBER_IN_FIRST_PENCIL',
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [CERT,inp,up,bp]},
            'checker_sha256':sha(Path(__file__)),
            'Sage_Mestre_covariants_source':str(sage_source),
            'Sage_Mestre_covariants_sha256':sha(sage_source),
            'exact_covariants_reproduced':True,'primitive_invariant_degree':45,
            'prime_trials_replayed':len(d['modular_trials']),
            'no_root_prime':127,'no_root_polynomial_mod127':last['coefficients'],
            'normalization_factors':list(map(str,factors)),
            'even_sextic_and_GL2_regressions':True,
            'conclusion':'Every rational affine member of the first-witness RR pencil has nonzero degree15 invariant and hence no degree2 elliptic quotient over the algebraic closure.',
            'point_searches':0,'Selmer_runs':0,'class_group_runs':0,
            'boundary':'Higher-degree elliptic quotients, other RR pencils and a predictive class construction remain open. This retrospective invariant-zero obstruction is not a solubility or Selmer certificate.'}


if __name__=='__main__':
    d=verify()
    text=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if OUT.exists():
        assert OUT.read_text()==text
    else:
        OUT.write_text(text)
    print(d['status'],flush=True)
