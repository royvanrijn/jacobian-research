#!/usr/bin/env sage-python
"""Independent common-quadric and integral-lattice replay.

Companion matrices only: no number-field or ideal backend, producer imports,
point search, factorization, unit group or class group.25s cap.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,lcm
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_seed_norm_quadric_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
generic=read(OUT/'generic-quadric.json');summary=read(OUT/'summary.json')
for data in [read(OUT/'protocol.json'),generic,summary]:
    for name,digest in data['inputs'].items():assert sha(ROOT/name)==digest
arithmetic=read(ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json')
R=PolynomialRing(QQ,'z');f=R(arithmetic['cubic_ascending'])
assert list(map(str,f.list()))==generic['cubic'] and f.degree()==3 and f.is_monic()
assert f.change_ring(GF(31)).is_irreducible()
T=matrix(QQ,[[0,0,-f[0]],[1,0,-f[1]],[0,1,-f[2]]])
one=vector(QQ,[1,0,0]);theta=vector(QQ,[0,1,0]);Id=matrix.identity(QQ,3)
def mult(v):return v[0]*Id+v[1]*T+v[2]*T*T
def tr(v):return mult(v).trace()
def norm(v):return mult(v).det()
def inverse(v):return mult(v).solve_right(one)
def product(v,w):return mult(v)*w
e1=theta-tr(theta)/3*one;e2=T*T*one-tr(T*T*one)/3*one
assert [list(map(str,e)) for e in [e1,e2]]==generic['trace_zero_basis']
Q=matrix(QQ,generic['quadric_Gram'])
assert Q[:2,:2]==matrix(QQ,2,2,lambda i,j:-tr(product([e1,e2][i],[e1,e2][j]))/2)
assert Q[2:,2:]==matrix(QQ,[[1,QQ(1)/2],[QQ(1)/2,1]])
assert all(Q[i,j]==0 for i in range(2) for j in range(2,4)) and Q==Q.transpose()
assert Q.det()/f.discriminant()==QQ(1)/16 and not Q.det().is_square()
assert f.discriminant()>0 and Q[:1,:1].det()<0 and Q[:2,:2].det()>0
u0=T**3*one/norm(theta)
assert list(map(str,u0))==generic['norm_one_projection_source'] and norm(u0)==1
zeta=inverse(u0-one);b=tr(zeta)/3;v=zeta-b*one
p=vector(QQ,generic['quadric_projection_point'])
assert p[0]*e1+p[1]*e2==v and p[2]==b+1 and p[3]==b and p*Q*p==0
assert product(v+p[2]*one,inverse(v+p[3]*one))==u0
def polynomial(ring,rows):return sum(QQ(row['coefficient'])*ring.monomial(*row['exponents']) for row in rows)
P=PolynomialRing(QQ,['v1','v2','a','b']);v1,v2,aa,bb=P.gens()
qvec=vector(P,P.gens());quadric=qvec*Q*qvec
vp=v1*vector(P,e1)+v2*vector(P,e2)
plus=vp+vector(P,[aa,0,0]);minus=vp+vector(P,[bb,0,0])
def polymult(v):return v[0]*Id+v[1]*T+v[2]*T*T
den=polynomial(P,generic['torus_denominator'])
num=vector(P,[polynomial(P,row) for row in generic['torus_numerator']])
assert den==polymult(minus).det()
assert polymult(minus)*num==den*plus
assert polymult(plus).det()-den==(aa-bb)*quadric
assert polymult(num).det()==den**2*(den+(aa-bb)*quadric)
S=PolynomialRing(QQ,['r','s','h']);r,s,h=S.gens()
j=generic['projection_hyperplane_zero_coordinate'];assert p[j]
assert j==next(i for i,c in enumerate(p) if c)
y=vector(S,[0]*4)
for i,var in zip([i for i in range(4) if i!=j],S.gens()):y[i]=var
Phi=vector(S,[polynomial(S,row) for row in generic['quadric_parameterization']])
assert Phi==(y*Q*y)*vector(S,p)-2*(vector(S,p)*Q*y)*y
assert Phi*Q*Phi==0
assert Phi-Phi[j]/p[j]*vector(S,p)==-2*(vector(S,p)*Q*y)*y
# Conversely the second-intersection formula recovers every Q-point away
# from the tangent plane; the exact residual is q(X)*p.
projected=qvec-qvec[j]/p[j]*vector(P,p)
recovered=(projected*Q*projected)*vector(P,p)-2*(vector(P,p)*Q*projected)*projected
assert recovered+2*(vector(P,p)*Q*qvec)*qvec==quadric*vector(P,p)
zk=[vector(QQ,R(text).list()+[0]*(3-len(R(text).list()))) for text in arithmetic['maximal_order_basis']]
Z=matrix(QQ,zk).transpose();cases=[]
for index in range(2):
    data=read(OUT/('integral-lattice-%02d.json'%index))
    for name,digest in data['inputs'].items():assert sha(ROOT/name)==digest
    relative=read(ART/('det1092_seed_artin_v2/relative-%02d.json'%index))
    ideal=matrix(QQ,relative['ideal']);B=Z*ideal;N=QQ(relative['norm'])
    beta=vector(QQ,relative['beta']);a0=N*inverse(beta)
    assert norm(beta)==N*N and norm(a0)==N and abs(ideal.det())==N
    assert B==matrix(QQ,data['power_basis_matrix'])
    assert list(map(str,a0))==data['rational_norm_witness']
    L=matrix(QQ,data['torus_to_norm_coordinates'])
    assert B*L==mult(a0)
    scale=ZZ(data['integrality_scale']);C=matrix(ZZ,data['integrality_integer_matrix'])
    assert scale==lcm(v.denominator() for v in L.list()) and C==scale*L
    for row,u in zip(data['exceptional_points'],[one,u0]):
        coordinates=L*u
        assert list(map(str,coordinates))==row['norm_form_coordinates']
        assert norm(B*coordinates)==N
        assert row['integral']==all(v.denominator()==1 for v in coordinates)
        assert not row['integral']
    # Exact composition into the old normalized norm polynomial, over the
    # polynomial ring rather than by testing a finite parameter sample.
    lifted=B*vector(P,L*num)
    assert polymult(lifted).det()==N*den**2*(den+(aa-bb)*quadric)
    cases.append({'index':index,'integral_transport_verified':True,
                  'identity_and_projection_exceptions_integral':False})
report={'status':'PASS_INDEPENDENT_COMMON_QUADRIC_AND_INTEGRAL_ADMISSIBILITY',
    'classification':'new explicit reduction and verified application, no principality solution',
    'generic_quadric_rank':4,'generic_quadric_signature':[2,2],
    'quadric_discriminant_ratio':'1/16','rational_quadric_point_from_cubic_only':True,
    'complete_norm_one_chart_with_explicit_exceptions':True,
    'cases':cases,
    'written_scope':'Nonsquare discriminant excludes a second rational point in the tangent plane, so projection misses only its rational center. All other nonidentity norm-one elements have the explicit inverse. The two omitted norm-form points are verified nonintegral for each ideal.',
    'boundary':'The common rational model forgets the ideal class unless the integral matrices are retained. Neither integral norm equation, unit-versus-ideal fork, prospective302 seed, or elliptic Selmer/Sha class is constructed.',
    'inputs':{str(path.relative_to(ROOT)):sha(path) for path in [OUT/'protocol.json',OUT/'generic-quadric.json',
        *[OUT/('integral-lattice-%02d.json'%i) for i in range(2)],Path(__file__)]}}
dest=OUT/'independent-replay.json'
if dest.exists():assert read(dest)==report
else:
    with dest.open('x') as stream:json.dump(report,stream,indent=2,sort_keys=True);stream.write('\n')
print(report['status'],'two lattices, four nonintegral chart exceptions',flush=True)
