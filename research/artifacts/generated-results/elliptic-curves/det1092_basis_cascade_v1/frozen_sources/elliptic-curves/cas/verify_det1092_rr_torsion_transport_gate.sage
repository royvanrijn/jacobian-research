#!/usr/bin/env sage-python
"""Independent bit-subset and exhaustive small-linear-map replay."""
import hashlib,json,signal
from itertools import product
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,matrix,identity_matrix
from sage.version import version as sage_version
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
CERT=ART/'det1092_rr_torsion_transport_gate_v1.json'
OUT=ART/'det1092_rr_torsion_transport_gate_replay_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    d=json.loads(CERT.read_text())
    for p,h in d['inputs'].items():assert sha(ROOT/p)==h
    src=ART/'det1092_marked_kummer_transport_v1.json';source=json.loads(src.read_text())
    assert d['status']=='PASS_TEN_RR_JACOBIANS_Q_SIMPLE_AND_NO_ELLIPTIC_TORSION_TRANSPORT'
    # Represent 2-torsion by even subsets modulo complementation, independently
    # of the constructor's four-coordinate quotient representation.
    even=[n for n in range(64) if n.bit_count()%2==0]
    def canon(n):return min(n,n^63)
    classes=sorted({canon(n) for n in even});assert len(classes)==16
    perm=d['permutation'];assert perm==[1,2,3,4,0,5]
    def action(n):return canon(sum(((n>>i)&1)<<perm[i] for i in range(6)))
    unseen=set(classes);orbits=[]
    while unseen:
        start=min(unseen);orbit=[];n=start
        while n not in orbit:orbit.append(n);n=action(n)
        assert n==start;unseen.difference_update(orbit);orbits.append(orbit)
    assert sorted(map(len,orbits))==[1,5,5,5]
    # A proper nonzero invariant subspace has 1,3,or7 nonzero vectors. None
    # can be a union of nonzero orbits, all of which have length5.
    assert all((2**r-1)%5 for r in [1,2,3])
    F=GF(2);R=PolynomialRing(F,'x');x=R.gen();J=matrix(F,d['J2_action'])
    basis=[(1<<i)|(1<<5) for i in range(4)]
    def encode(bits):
        n=0
        for b,e in zip(bits,basis):
            if b:n^=e
        return canon(n)
    for bits in product([0,1],repeat=4):
        image=J*matrix(F,4,1,bits)
        assert encode(image.list())==action(encode(bits))
    polynomial=R(d['characteristic_polynomial'])
    assert polynomial==x**4+x**3+x**2+x+1==R(J.charpoly().list())
    assert polynomial.gcd(x**2-x)==1 and polynomial.gcd(x**4-x)==1
    assert J**5==identity_matrix(F,4) and J!=identity_matrix(F,4)
    records=[];expected_actions=[]
    for bits in product([0,1],repeat=4):
        E=matrix(F,2,2,bits)
        if E.det():expected_actions.append(E)
    assert len(expected_actions)==len(d['elliptic_intertwining_tests'])==6
    # Only 256 maps in each direction per action, not a subgroup-state census.
    for E,saved in zip(expected_actions,d['elliptic_intertwining_tests']):
        assert E==matrix(F,saved['elliptic_action'])
        forward=0;backward=0
        for bits in product([0,1],repeat=8):
            T=matrix(F,2,4,bits);U=matrix(F,4,2,bits)
            if T*J==E*T:assert not T;forward+=1
            if J*U==U*E:assert not U;backward+=1
        assert (forward,backward)==(1,1)
        assert matrix(F,saved['J_to_E_constraints']).rank()==saved['J_to_E_rank']==8
        assert matrix(F,saved['E_to_J_constraints']).rank()==saved['E_to_J_rank']==8
        assert saved['Hom_dimensions']==[0,0]
        records.append({'elliptic_action':saved['elliptic_action'],'maps_tested_per_direction':256,
                        'forward_maps':forward,'backward_maps':backward,'only_zero_maps':True})
    cases=[];assert len(d['cases'])==len(source['cases'])==10
    for saved,original in zip(d['cases'],source['cases']):
        assert saved['q']==original['q'] and saved['label']==original['label']
        p=saved['p'];assert p==original['theta_primes']['one_plus_five']
        R=PolynomialRing(GF(p),'T');T=R.gen();q=R([QQ(a) for a in saved['q']])
        assert q.degree()==6 and q.gcd(q.derivative())==1
        factors=[]
        for item in saved['mod_p_factors']:
            assert item['multiplicity']==1;g=R(item['coefficients']);assert g.is_monic()
            if g.degree()==1:assert g.gcd(q)==g
            else:
                assert g.degree()==5
                # Rabin's degree-five irreducibility test, no factor() call.
                assert pow(T,p**5,g)==T%g and g.gcd(pow(T,p,g)-T)==1
            factors.append(g)
        assert sorted(g.degree() for g in factors)==[1,5]
        assert q==q.leading_coefficient()*factors[0]*factors[1]
        assert saved['J2_GQ_irreducible'] and saved['Q_simple']
        cases.append({'label':saved['label'],'p':p,'degree_pattern':[1,5],
                      'J2_irreducible':True,'no_elliptic_factor_over_Q':True})
    return {'classification':'verified application and new deduction',
            'status':'PASS_INDEPENDENT_TEN_RR_TORSION_TRANSPORT_OBSTRUCTIONS',
            'sage_version':sage_version,'subset_orbits':orbits,'intertwining_tests':records,'cases':cases,
            'new_degree_extension_deduction':{
                'statement':'Irreducibility and the no-elliptic-factor obstruction persist over every number field L whose degree over Q is not divisible by5.',
                'proof':'The image of G_L in the finite J[2] Galois image has index dividing [L:Q]. If5 does not divide that degree, its order is still divisible by5. Cauchy gives an element of order5. In dimension4 over F2 any such element has irreducible minimal polynomial Phi5, so J[2] remains irreducible.',
                'necessary_degree_for_an_elliptic_factor':'A multiple of5, not a sufficient condition.'},
            'scope':'No nonzero coefficient-induced H1/Selmer transport to elliptic2-torsion, and no elliptic quotient over Q or the stated extensions. This does not compute Selmer groups, prove absolute simplicity, or rule out nonlinear incidence predicates or other correspondences.',
            'limits':{**d['limits'],'finite_linear_maps_tested':3072,'point_searches':0},
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [CERT,src]},
            'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    signal.alarm(25);result=verify();payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==payload
    else:OUT.write_text(payload)
    print(result['status'],flush=True)
