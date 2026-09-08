#!/usr/bin/env sage-python
"""Independent chord elimination and cubic multiplication-table replay.

No constructor import, algebra square test, number field, local descent,
exceptional point, point search or later search artifact. 25-second cap.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve,matrix,vector

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_trace_twist_kummer_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def provenance(d):
    for path,digest in d['inputs'].items():assert sha(ROOT/path)==digest,(path,'HASH_MISMATCH')
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def product(P,Q,abc):
    a,b,c=abc
    coeff=[QQ(0)]*5
    for i in range(3):
        for j in range(3):coeff[i+j]+=P[i]*Q[j]
    for j in [4,3]:
        value=coeff[j]
        coeff[j-1]-=a*value;coeff[j-2]-=b*value;coeff[j-3]-=c*value
    return vector(QQ,coeff[:3])
def multiplication(P,abc):
    basis=[vector(QQ,[int(i==j) for i in range(3)]) for j in range(3)]
    return matrix(QQ,[product(P,b,abc) for b in basis]).transpose()
def verify():
    paths=[DIR/f for f in ['protocol.json','universal.json','construction.json','panel.json']]
    protocol,universal,data,panel=map(read,paths)
    for d in [protocol,data,panel]:provenance(d)
    assert protocol['limits']['old_addresses']==9 and protocol['limits']['new_addresses']==0
    assert protocol['limits']['point_searches']==protocol['limits']['exceptional_point_inputs']==0
    assert universal['variables']==['a','m','u','v','d','n']
    parent=read(ART/'curve302_recovered_mw17_parent_v1.json')
    cover=read(ART/'det1092_orbit8044_rank18_base_change_v2.json');provenance(cover)
    roster=read(ART/'det1092_rr_generic_point_controls_v2/protocol.json')

    # Verify the universal identities in a polynomial ring, not by sampling
    # parameters or using a square-class routine. No rational denominators.
    B=PolynomialRing(QQ,names=('a','m','u','v','d','n'))
    a,m,u,v,d,n=B.gens();W=PolynomialRing(B,'x');x=W.gen()
    g=x*x-2*u*x+u*u-v*v*d;c=m*m-a-2*u
    L=m*x+n;F=L*L+(x-c)*g
    A=(m*u+n)**2-v*v*d*(a+2*u)
    N=(m*u+n)*(x-u)+d*m*v*v
    assert F.degree()==3 and F[3]==1 and F[2]==a
    chord_residual=d*v*v*F-N*N-g*(d*v*v*x-A)
    assert not chord_residual
    # This is an exact polynomial identity, stronger than equality mod F.
    identity=(A-v*v*d*x)*L*L-(c-x)*N*N
    assert identity==(A-v*v*d*c)*F

    T=PolynomialRing(QQ,'t');K=T.fraction_field()
    def dec(d):return K(T(d['numerator']))/T(d['denominator'])
    ai=list(map(dec,parent['a_invariants']));E=EllipticCurve(K,ai)
    a,b,c0=[E.b2(),8*E.b4(),16*E.b6()]
    assert [a,b,c0]==list(map(dec,data['cubic_coefficients']))
    d=dec(data['d']);assert d==T(cover['curve_over_Q']['q_coefficients'])
    rr_c,rr_b,rr_a=map(dec,cover['lift']['residual_coefficients'])
    h=K(T(cover['splitting']['discriminant_square_factor']))
    assert rr_b*rr_b-4*rr_a*rr_c==h*h*d
    f0,f1,f2=[K(T(v)) for v in cover['lift']['line_coefficients']]
    u,v=map(dec,data['quadratic_point']['X']);w,z=map(dec,data['quadratic_point']['Y'])
    assert u==-2*rr_b/rr_a and v==2*h/rr_a
    y0=-(f0+f1*u/4)/f2;y1=-f1*v/(4*f2)
    assert w==8*y0+ai[0]*u+4*ai[2] and z==8*y1+ai[0]*v
    # Two coefficients of the elliptic equation in K[s]/(s^2-d).
    assert w*w+z*z*d==u**3+a*u*u+b*u+c0+(3*u+a)*v*v*d
    assert 2*w*z==(3*u*u+2*a*u+b)*v+v**3*d
    n,m=map(dec,data['line']);assert z==m*v and w==m*u+n and v
    cx,cy=map(dec,data['trace']);assert cx==m*m-a-2*u and cy==-(m*cx+n)
    basis=[E(list(map(dec,P))) for P in parent['basis_weierstrass_coordinates']]
    word=cover['lift']['trace_word'];assert data['trace_word']==word
    Z=sum((int(j)*P for j,P in zip(word,basis)),E(0))
    assert cx==4*Z[0] and cy==8*Z[1]+4*ai[0]*Z[0]+4*ai[2]
    r,eta=map(dec,data['anti_twist_point'])
    assert r==(w/v)**2/d-a-2*u and eta==(w/v)*(u-r)/d-z
    assert d*eta*eta==r**3+a*r*r+b*r+c0
    NN=list(map(dec,data['square_numerator']));LL=list(map(dec,data['square_denominator']))
    assert NN==[-w*u/v+d*z,w/v] and LL==[n,m]
    # All displayed maps are exact rational functions; record their finite
    # affine exclusions rather than silently dropping undefined parameters.
    expressions=[*ai,a,b,c0,d,u,v,w,z,cx,cy,r,eta,*NN,*LL]
    exclusions=[f.denominator() for f in expressions]
    exclusions += [f.numerator() for f in [d,v,cy,eta,K(E.discriminant())]]
    assert all(exclusions)

    assert len(panel['cases'])==len(roster['cases'])==9
    cases=[]
    for i,(case,row) in enumerate(zip(roster['cases'],panel['cases'])):
        assert row==read(DIR/f'case-{i:02d}.json')
        assert row['index']==i and row['label']==case['label'] and row['parameter']==case['parameter']
        tau=QQ(case['parameter']);assert all(p(tau) for p in exclusions)
        abc=[f(tau) for f in [a,b,c0]];dv=d(tau);rv,ev=r(tau),eta(tau)
        cxv,cyv=cx(tau),cy(tau)
        assert list(map(str,abc))==row['cubic_coefficients']
        assert row['d']==str(dv) and row['trace']==list(map(str,[cxv,cyv]))
        assert row['anti_twist_point']==list(map(str,[rv,ev]))
        gamma=list(map(QQ,row['square_witness']));gamma+= [QQ(0)]*(3-len(gamma))
        gamma=vector(QQ,gamma);assert len(gamma)==3
        Lv=vector(QQ,[f(tau) for f in LL]+[0]);Nv=vector(QQ,[f(tau) for f in NN]+[0])
        betaZ=vector(QQ,[cxv,-1,0]);betaR=vector(QQ,[dv*rv,-dv,0])
        ML,MN,MZ=[multiplication(f,abc) for f in [Lv,Nv,betaZ]]
        assert ML.det() and MN.det() and MZ.det()
        assert ML*gamma==Nv
        assert product(betaZ,product(gamma,gamma,abc),abc)==betaR
        assert MZ.det()==cyv*cyv
        assert multiplication(betaR,abc).det()==dv**4*ev*ev
        # Rational trace belongs to the displayed generic subgroup by the
        # exact function-field word identity above. No rank oracle is used.
        sn,sd=map(int,row['nonsquare_bounds']);dn,dd=dv.numerator(),dv.denominator()
        assert dv>0 and sn*sn<=dn<(sn+1)**2 and sd*sd<=dd<(sd+1)**2
        assert sn*sn!=dn or sd*sd!=dd
        cases.append({'index':i,'label':case['label'],
          'unit_norms':list(map(str,[ML.det(),MN.det(),MZ.det()])),
          'trace_word':word,'cover_split':False,
          'relative_twist_to_elliptic_2_kummer_class':'ZERO_INHERITED_TRACE',
          'elliptic_seed_existence':'UNKNOWN_NOT_EXCLUDED'})
    result={'status':'PASS_INDEPENDENT_TRACE_TWIST_KUMMER_OBSTRUCTION',
      'classification':'new explicit algebraic deduction and verified application',
      'universal_polynomial_identity':True,
      'statement':'For a quadratic point Q with rational trace Z, the anti-trace R=Q-sigma(Q) is rational on the d-twist, and its norm-square class d*(x_R-theta) equals (x_Z-theta) times the displayed square. For this generic conic it is inherited on all nine frozen fibres, not a non-generic elliptic seed class.',
      'cases':cases,'limits':protocol['limits'],
      'excluded_affine_parameters':[list(map(str,p.list())) for p in exclusions],
      'boundary':'Obstructs this direct 2-torsion-module transfer, not rational points on either twist, other Kummer classes, higher descent, or seeds from a split branch. No Selmer dimension or nonzero Sha class is inferred.',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[Path(__file__)]}}
    retain(DIR/'replay.json',result)
    print(result['status'],'nine inherited classes; no exceptional inputs',flush=True)

if __name__=='__main__':
    signal.alarm(25);verify()
