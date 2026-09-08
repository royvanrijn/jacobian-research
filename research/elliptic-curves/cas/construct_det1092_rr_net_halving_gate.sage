#!/usr/bin/env sage-python
"""Generic-input halving curve and singular-member incidence of the RR net.

One explicit quartic, one discriminant, <=18 small finite-field root tests.
No function-field genus algorithm, class groups, point searches, or sweeps.
"""
import hashlib,json,runpy
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
OUT=ART/'det1092_rr_net_halving_gate_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    pp=ART/'curve302_recovered_mw17_parent_v1.json';np=ART/'det1092_first_centre_rr_net_v1.json';loader=CAS/'load_curve302_recovered_parent.sage'
    parent=json.loads(pp.read_text());net=json.loads(np.read_text())
    E,base,_=runpy.run_path(str(loader))['load_curve302_recovered_parent'](pp)
    K=E.base_ring();R=K.ring();t=R.gen();w=-vector(ZZ,net['trace_word'])
    C=sum((n*P for n,P in zip(w,base)),E(0))
    a1,a2,a3,a4,a6=E.a_invariants();b2=E.b2()
    A=-E.c4()/48;B=-E.c6()/864;cx=C[0]+b2/12;cy=C[1]+(a1*C[0]+a3)/2
    assert cy*cy==cx**3+A*cx+B
    M=PolynomialRing(K,'m');m=M.gen()
    D=m**4-6*cx*m*m-8*cy*m-3*cx*cx-4*A
    assert D.discriminant()==256*E.discriminant()
    delta=R(E.discriminant());assert delta.degree()==24 and delta.gcd(delta.derivative()).degree()==0
    assert delta.gcd(R(cx.denominator()*cy.denominator())).degree()==0
    assert cx.numerator().degree()-cx.denominator().degree()<=4
    assert cy.numerator().degree()-cy.denominator().degree()<=6
    assert any(n%2 for n in w) and matrix(QQ,parent['generic_height_gram']).det()==1092
    AA=list(map(R,net['A']));BB=list(map(R,net['B']))
    N=-(BB[1]+(m-a1/2)*BB[2]);L=AA[1]+(m-a1/2)*AA[2]
    assert N.degree()==L.degree()==1
    assert N[1]*L[0]-N[0]*L[1]
    assert D.gcd(D.derivative()).degree()==0 and D.gcd(L).degree()==0
    primes=[17,47,53,61,67,71,79,83,89,101,107,113,127,137,149,179,191,197]
    witness=None
    for p in primes:
        vals=[QQ(v(0)) for v in D.list()]
        if any(v.denominator()%p==0 for v in vals):continue
        F=PolynomialRing(GF(p),'m');fp=F(vals)
        if fp.discriminant() and not fp.roots():
            witness={'prime':p,'quartic_coefficients':list(map(int,fp.list()))};break
    assert witness is not None
    def rec(f):return {'numerator':list(map(str,f.numerator().list())),'denominator':list(map(str,f.denominator().list()))}
    return {'classification':'verified application and new deduction','status':'PASS_EXPLICIT_HALVING_CURVE_AND_ZERO_LOCAL_OBSTRUCTION',
        'centre_word':list(map(int,w)),'short_A':rec(A),'short_B':rec(B),'centre_short':[rec(cx),rec(cy)],
        'halving_quartic_coefficients':[rec(v) for v in D.list()],
        'half_point_map':{'x':'(m^2-cx)/2','y':'m*(x-cx)-cy','assertion':'2Q=C on D(t,m)=0'},
        'discriminant_identity':'disc_m(D)=256*Delta_E',
        'net_r_numerator':[rec(v) for v in N.list()],'net_r_denominator':[rec(v) for v in L.list()],
        'singular_member_map':{'r':'N/L','v':'partial_t(r)-partial_m(r)*partial_t(D)/partial_m(D)','u':'r-t*v'},
        'open_conditions':['smooth parent fibre','finite centre and net coefficients','L != 0','partial_m(D) != 0','partial_m(r) != 0'],
        'geometric_halving_cover':{'degree':4,'connected':'Full geometric MW basis, centre word not even, and no geometric rational2torsion exclude components of degrees1 or2.',
            'branch_points':24,'ramification_profile_at_each':[2,1,1],'genus':9,
            'no_other_ramification':'Multiplication by2 is etale on smooth fibres; infinity is smooth.'},
        'zero_no_half_witness':witness,
        'limits':{'halving_curves':1,'quartic_discriminants':1,'finite_prime_tests_max':18,'finite_prime_bound':197,'point_searches':0,'parameter_sweeps':0},
        'boundary':'Dense-open singular-member incidence only. Exceptional net members, reducible curves and singularities over excluded fibres are not classified. This does not explain the independent302 point.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [pp,np,loader,Path(__file__)]}}
if __name__=='__main__':
    if OUT.exists():raise FileExistsError(OUT)
    d=build();OUT.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],d['zero_no_half_witness'],flush=True)
