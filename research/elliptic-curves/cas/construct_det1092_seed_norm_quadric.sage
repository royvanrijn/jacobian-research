#!/usr/bin/env sage-python
"""A generic cubic norm-one quadric; two retrospective integral lattices.

25s cap. No unit/class group, point search, factorization or norm search.
The generic quadric and projection point are fixed before relative ideals.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,NumberField,matrix,vector,lcm
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_seed_norm_quadric_v1'
ARITH=ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
REL=[ART/('det1092_seed_artin_v2/relative-%02d.json'%i) for i in range(2)]
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists():assert read(p)==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def rows(M):return [[str(x) for x in row] for row in M.rows()]
def coefficients(a):return [str(a[i]) for i in range(3)]
def polyrecord(f):return [{'exponents':list(map(int,e)),'coefficient':str(c)} for e,c in sorted(f.dict().items())]
OUT.mkdir(exist_ok=True)
save('protocol.json',{'classification':'generic rational quadric, then retrospective integral admissibility',
    'rule':'Build the trace-zero norm-one quadric from the cubic alone. Use u0=theta^3/Norm(theta) as the projection point source, never a missing point or unit search. Then transport only the two already frozen relative norm equations to that same quadric. Retain the unit and projection-point exceptions explicitly; do not replace a failed constructor.',
    'limits':{'seconds':25,'quadrics':1,'relative_ideals':2,'class_groups':0,'unit_groups':0,
              'point_searches':0,'norm_searches':0,'factorizations':0,'later_cascade_inputs':0},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [ARITH,*REL,Path(__file__)]}})
arithmetic=read(ARITH);R=PolynomialRing(QQ,'z');f=R(arithmetic['cubic_ascending'])
K=NumberField(f,'theta');theta=K.gen()
e1=theta-theta.trace()/3;e2=theta**2-(theta**2).trace()/3
zero_basis=[e1,e2];assert e1.trace()==e2.trace()==0
trace=matrix(QQ,2,2,lambda i,j:(zero_basis[i]*zero_basis[j]).trace())
Q=matrix(QQ,4,4)
Q[:2,:2]=-trace/2;Q[2,2]=Q[3,3]=1;Q[2,3]=Q[3,2]=QQ(1)/2
assert Q.det()!=0 and not Q.det().is_square()
assert Q.det()/f.discriminant()==QQ(1)/16
u0=theta**3/theta.norm();assert u0.norm()==1 and u0!=1
zeta=1/(u0-1);b=zeta.trace()/3;a=b+1;v=zeta-b
V=matrix(QQ,[e1.vector(),e2.vector()]).transpose()
coords=V.solve_right(v.vector());p=vector(QQ,[*coords,a,b])
assert p*Q*p==0 and (v+a)/(v+b)==u0
projection_index=next(i for i,c in enumerate(p) if c)
S=PolynomialRing(QQ,['r','s','h']);parameters=S.gens()
other=[i for i in range(4) if i!=projection_index]
y=vector(S,[0]*4)
for i,value in zip(other,parameters):y[i]=value
Qy=y*Q*y;Bpy=vector(S,p)*Q*y
Phi=Qy*vector(S,p)-2*Bpy*y
assert Phi*Q*Phi==0 and Phi[projection_index]==Qy*p[projection_index]
P=PolynomialRing(QQ,['v1','v2','a','b']);v1,v2,aa,bb=P.gens()
vp=vector(P,[v1*e1[i]+v2*e2[i] for i in range(3)])
def multiplication(coords):
    return matrix(P,3,3,lambda i,j:sum(coords[k]*(theta**(k+j))[i] for k in range(3)))
plus=vp+vector(P,[aa,0,0]);minus=vp+vector(P,[bb,0,0])
M=multiplication(minus);den=M.det();num=M.adjugate()*plus
quadratic=vector(P,P.gens())*Q*vector(P,P.gens())
assert multiplication(plus).det()-den==(aa-bb)*quadratic
generic={'status':'PASS_GENERIC_NORM_ONE_QUADRIC_AND_PARAMETRIZATION',
    'classification':'generic equation-only model; no relative ideal or seed used',
    'cubic':list(map(str,f.list())),'trace_zero_basis':[coefficients(e1),coefficients(e2)],
    'quadric_Gram':rows(Q),'determinant_discriminant_ratio':str(Q.det()/f.discriminant()),
    'norm_one_projection_source':coefficients(u0),'quadric_projection_point':list(map(str,p)),
    'projection_hyperplane_zero_coordinate':projection_index,
    'parameter_order':['r','s','h'],'quadric_parameterization':[polyrecord(x) for x in Phi],
    'quadric_coordinate_order':['v1','v2','a','b'],
    'torus_numerator':[polyrecord(x) for x in num],'torus_denominator':polyrecord(den),
    'inverse':'For u!=1, zeta=1/(u-1), b=Tr(zeta)/3, a=b+1, v=zeta-b. Express v in the displayed trace-zero basis.',
    'projection_exception':'The only rational point in the tangent plane at p is p, since det(Q) is nonsquare. All other rational quadric points are covered by the displayed projection chart.',
    'torus_exceptions':'a=b contracts to u=1 where defined; Norm(v+b)=0 is excluded. Treat u=1 and u=u0 separately in integral admissibility.',
    'inputs':{str(ARITH.relative_to(ROOT)):sha(ARITH)}}
save('generic-quadric.json',generic)
print(generic['status'],'signature(2,2), discriminant ratio1/16',flush=True)
# Only now use the two seed-relative lattices, never to choose the quadric.
zk=[K(R(x)) for x in arithmetic['maximal_order_basis']]
Z=matrix(QQ,[x.vector() for x in zk]).transpose()
cases=[]
for i,path in enumerate(REL):
    relative=read(path);ideal=matrix(QQ,relative['ideal']);power_basis=Z*ideal
    beta=K(R(relative['beta']));N=QQ(relative['norm']);a0=K(N)/beta
    assert abs(ideal.det())==N and a0.norm()==N
    mult_a0=matrix(QQ,3,3,lambda k,j:(a0*theta**j)[k])
    lattice=power_basis.inverse()*mult_a0
    scale=ZZ(lcm(c.denominator() for c in lattice.list()))
    integral_matrix=matrix(ZZ,scale*lattice)
    exceptional=[]
    for name,u in [('identity',K.one()),('projection_source',u0)]:
        point=lattice*u.vector()
        assert K(sum((point[j]*sum((power_basis[k,j]*theta**k for k in range(3)),K.zero()) for j in range(3)),K.zero())).norm()==N
        exceptional.append({'name':name,'norm_form_coordinates':list(map(str,point)),
                            'integral':all(x.denominator()==1 for x in point)})
    case={'index':i,'ideal_norm':str(N),'rational_norm_witness':coefficients(a0),
          'power_basis_matrix':rows(power_basis),'torus_to_norm_coordinates':rows(lattice),
          'integrality_scale':str(scale),'integrality_integer_matrix':rows(integral_matrix),
          'exceptional_points':exceptional,
          'criterion':'On a rational quadric point outside the listed denominators, multiply the three cubic torus numerators by the integer matrix. All three quotients by integrality_scale*torus_denominator must be integers. This is equivalent to the original integral norm equation; parametrization alone does not impose it.',
          'inputs':{str(path.relative_to(ROOT)):sha(path)}}
    save('integral-lattice-%02d.json'%i,case);cases.append(case)
save('summary.json',{'status':'PASS_COMMON_QUADRIC_TWO_DISTINCT_INTEGRAL_LATTICES',
    'classification':'exact reduction, neither principality question solved',
    'rational_model_depends_only_on_cubic':True,
    'relative_ideal_information_retained_in_integrality_matrices':True,
    'exceptional_integral_solutions':[any(p['integral'] for p in c['exceptional_points']) for c in cases],
    'boundary':'The common rational quadric supplies no seed, unit, ideal relation or rank gain. Dropping the integral lattice makes every two-torsion ideal case rationally soluble.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'generic-quadric.json',*[OUT/('integral-lattice-%02d.json'%i) for i in range(2)]]}})
print('PASS_COMMON_QUADRIC_TWO_DISTINCT_INTEGRAL_LATTICES',
      [[p['integral'] for p in c['exceptional_points']] for c in cases],flush=True)
