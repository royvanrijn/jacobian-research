#!/usr/bin/env sage-python
"""Independent exact replay of the degree24 polynomial lift obstruction.

Rebuilds the cubic; certifies discriminant irreducibility from finite-field
Frobenius/gcd tests; verifies supplied square-root truncations and Bezout
identities. No constructor import, rational factorization, point search,
number field, class group, or Selmer dimension calculation.25-second cap.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,prod,power_mod
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_polynomial_lift_gate_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def provenance(d):
    for path,digest in d['inputs'].items():assert sha(ROOT/path)==digest,(path,'HASH_MISMATCH')
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True,default=int)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def finite_irreducibility(g,p):
    t=g.parent().gen();n=g.degree()
    assert g.is_monic() and n>0
    assert power_mod(t,p**n,g)==t%g
    for k in range(1,n):
        if n%k==0:assert (power_mod(t,p**k,g)-t).gcd(g)==1
def verify():
    paths=[DIR/'protocol.json',DIR/'construction.json']
    protocol,data=map(read,paths)
    provenance(protocol);provenance(data)
    assert protocol['limits']['polynomial_degree_bound']==24
    assert protocol['limits']['reduction_primes']==[149,151]
    parent_path=ART/'curve302_recovered_mw17_parent_v1.json'
    geometry_path=ART/'curve302_parent_geometric_picard19_v1.json'
    roster_path=ART/'det1092_rr_generic_point_controls_v2/protocol.json'
    parent,geometry,roster=map(read,[parent_path,geometry_path,roster_path])
    assert geometry['status']=='PASS'
    assert geometry['geometric_generic_MW_rank']==geometry['arithmetic_generic_MW_rank']==17
    assert geometry['full_geometric_basis_is_displayed_rational_basis'] is True
    R=PolynomialRing(QQ,'t');K=R.fraction_field()
    def rat(d):return K(R(d['numerator']))/R(d['denominator'])
    ai=list(map(rat,parent['a_invariants']));assert ai[:3]==[1,1,1]
    a,b,c=R(5),R(16*ai[3]+8),R(64*ai[4]+16)
    S=PolynomialRing(R,'X');X=S.gen();f=X**3+a*X*X+b*X+c
    assert [a.degree(),b.degree(),c.degree()]==[0,8,12]
    Delta=f.discriminant().monic()
    assert Delta.degree()==24 and Delta.gcd(Delta.derivative())==1 and Delta(0)
    assert [list(map(str,q.list())) for q in [a,b,c,Delta]]==[data[k] for k in ['a','b','c','monic_discriminant']]
    f0=PolynomialRing(GF(31),'X')([q(0) for q in f.list()])
    finite_irreducibility(f0,31)
    fK=PolynomialRing(K,'X')(f.list())
    for point in parent['basis_weierstrass_coordinates']:
        xp,yp=map(rat,point)
        assert fK(4*xp)==(8*yp+4*xp+4)**2
    # Universal norm identity is independent of the specialized discriminant.
    U=PolynomialRing(QQ,['a','b','c','k']);au,bu,cu,ku=U.gens()
    V=PolynomialRing(U,'z');z=V.gen();fu=z**3+au*z*z+bu*z+cu
    du=ku**3+au*ku*ku+bu*ku+cu
    alpha=du*(ku-z)
    assert matrix(U,3,3,lambda i,j:((alpha*z**j)%fu)[i]).det()==du**4
    templates={r['label']:r for r in data['templates']}
    assert set(templates)=={'double','simple'}
    roots={key:R(list(map(QQ,r['root_remainder']))) for key,r in templates.items()}
    r,s=roots['double'],roots['simple']
    assert r.degree()==s.degree()==23
    assert ((2*(a*a-3*b))*r-(9*c-a*b))%Delta==0
    assert (s+a+2*r)%Delta==0
    assert (s-r).gcd(Delta)==1
    # Exact factorization of the cubic over Q[t]/Delta; no root enumeration.
    difference=f-(X-r)**2*(X-s)
    assert all(q%Delta==0 for q in difference.list())
    exact_templates={}
    for label,z0 in roots.items():
        row=templates[label];A,B,C=[R(list(map(QQ,row[k]))) for k in ['A','B','C']]
        assert A==3*z0+a and B==3*z0*z0+2*a*z0+b
        assert Delta*C==z0**3+a*z0*z0+b*z0+c
        assert [A.degree(),B.degree(),C.degree()]==[23,46,45]
        order=row['known_residual_w_order'];assert order==({'double':2,'simple':3}[label])
        W=PolynomialRing(R,'w');w=W.gen()
        Q=W(Delta*Delta)+w*A*Delta+w*w*B+w**3*C
        hhat=W(Delta)+w*A/2
        if label=='simple':
            h2,rem=(B-A*A/4).quo_rem(2*Delta)
            assert rem==0 and h2.degree()<=22
            hhat+=w*w*h2
        # Proves rational divisibility of every square-truncation residual
        # by w^order, rather than inferring it from modular experiments.
        assert all((Q-hhat*hhat)[j]==0 for j in range(order))
        exact_templates[label]=(A,B,C,order)
    possible=set(range(25));prime_reports=[]
    for row,p in zip(data['primes'],[149,151]):
        assert row['prime']==p and ZZ(p).is_prime(proof=True)
        assert row==read(DIR/f'prime-{p}.json')
        F=GF(p);T=PolynomialRing(F,'t');t=T.gen()
        for q in [a,b,c,Delta,r,s]:assert all(v.denominator()%p for v in q)
        dp=T(Delta);assert dp.degree()==24 and dp.is_monic() and dp.gcd(dp.derivative())==1
        factors=[T(v) for v in row['factors']]
        assert prod(factors)==dp
        for h in factors:finite_irreducibility(h,p)
        deg=[int(h.degree()) for h in factors]
        assert deg==row['factor_degrees']
        sums={0}
        for n in deg:sums|={k+n for k in sums}
        assert sorted(sums)==row['possible_rational_factor_degrees']
        possible&=sums
        W=PolynomialRing(F,'w');w=W.gen();V=PolynomialRing(W,'t');tv=V.gen()
        def cv(q):return V(T(q).list())
        gate_reports=[]
        assert [g['label'] for g in row['square_gates']]==['double','simple']
        for gate in row['square_gates']:
            A,B,C,order=exact_templates[gate['label']]
            assert all(v.denominator()%p for q in [A,B,C] for v in q)
            Q=cv(Delta)**2+w*cv(A*Delta)+w*w*cv(B)+w**3*cv(C)
            h=V([W(v) for v in gate['monic_square_root_truncation']])
            assert h.degree()==24 and h.leading_coefficient()==1
            assert all(h[i].degree()<=24-i for i in range(25))
            residual=Q-h*h
            assert residual.degree()==23
            normalized=[W(v) for v in gate['normalized_residuals']]
            multipliers=[W(v) for v in gate['bezout_multipliers']]
            assert len(normalized)==len(multipliers)==24
            assert residual==V([w**order*q for q in normalized])
            assert sum(a*b for a,b in zip(normalized,multipliers))==1
            assert normalized[23].degree()==25-order
            assert gate['attained_w_degree_bound']==25-order
            assert gate['degree_preserving_residual_t_power']==23
            assert gate['known_residual_w_order']==order
            # Total degree(Q)<=48 implies deg_w(residual[t^23])<=25.
            # Attainment after division by w^order proves that this polynomial
            # loses no degree mod p. Bezout=1 then excludes a common factor
            # over Q and hence any nonzero characteristic-zero square parameter.
            gate_reports.append({'label':gate['label'],'residual_divisibility_order':order,
                'degree_preservation_bound_attained':25-order,
                'bezout_identity_verified':True,
                'nonzero_characteristic_zero_square_parameters':'NONE'})
        prime_reports.append({'prime':p,'factor_degrees':deg,'square_gates':gate_reports})
    assert possible=={0,24} and data['discriminant_factor_degree_intersection']==[0,24]
    cases=[]
    assert len(data['cases'])==len(roster['cases'])==9
    for i,(case,original) in enumerate(zip(data['cases'],roster['cases'])):
        assert case['index']==i and case['label']==original['label'] and case['parameter']==original['parameter']
        assert Delta(QQ(case['parameter'])) and case['smooth_fibre'] is True
        cases.append({'label':case['label'],'parameter':case['parameter'],
            'theorem_applies_to_specialization':True,'fibre_seed_existence':'NOT_EXCLUDED'})
    result={'status':'PASS_INDEPENDENT_POLYNOMIAL_LIFT_OBSTRUCTION',
        'classification':'verified application and new deduction',
        'polynomial_abscissa_degree_bound':24,
        'discriminant_irreducible_over_Q':True,
        'discriminant_root_remainder_degrees':[23,23],
        'prime_replays':prime_reports,'cases':cases,
        'theorem':'A polynomial k in Q[t] of degree<=24 with alpha=f(k)*(k-theta) unramified at every finite good parameter fibre has f(k) a square in Q[t] and alpha inherited from the full generic MW17 subgroup.',
        'proof_dependency':'The already proved full torsion-free geometric MW17 basis is rational; its proof is hash-pinned, not recomputed by this replay.',
        'logical_reduction':['Good-fibre unramifiedness forces f(k)=constant*square or constant*Delta*square.',
            'Discriminant irreducibility and the degree23 root remainders exclude the Delta case below degree23; odd degree excludes degree23.',
            'Degree24 forces k=r+u*Delta or k=s+u*Delta with nonzero constant u. Four finite-field Bezout identities, with exact divisibility and degree preservation, exclude both templates over characteristic zero.',
            'The remaining constant-square case defines a geometric section. All such sections are rational and in the displayed MW17, so the constant twist must be square.'],
        'boundaries':data['limitations'],
        'limits':protocol['limits'],
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[parent_path,geometry_path,roster_path,Path(__file__)]}}
    retain(DIR/'replay.json',result)
    print(result['status'],'degree_bound=24','four_Bezout_gates','nine_existing_addresses',flush=True)
if __name__=='__main__':
    signal.alarm(25);verify()
