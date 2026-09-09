#!/usr/bin/env sage-python
"""Independent complete-group and homogeneous-congruence replay.

No constructor imports or saved Kummer character rows are used for rank.
The infinitely many rational terms are certified symbolically, not sampled.
"""
import hashlib, json, signal
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, gcd, prod
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_local_shadow_seeds_v2'
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
protocol=read(OUT/'protocol.json'); data=read(OUT/'construction.json')
for name,digest in protocol['inputs'].items(): assert sha(ROOT/name)==digest
assert protocol['script_sha256']==sha(ROOT/'elliptic-curves/cas/construct_det1092_local_shadow_seeds.sage')
assert data['protocol_sha256']==sha(OUT/'protocol.json')
source=read(ART/'det1092_conic_seed_progression_v1/input.json')
frame=read(ART/'det1092_split_descent_v1/dependent-conic/generic-frame.json')
R=PolynomialRing(QQ,'u'); F=R.fraction_field()
def dec(v): return F(R(v['numerator']))/R(v['denominator'])
ai=list(map(dec,source['weierstrass_functions'])); E=EllipticCurve(F,ai)
original=[[dec(v) for v in P] for P in source['point_functions']]
k=QQ(data['constant_short_model_scale'])
AB=[-k**4*E.c4()/48,-k**6*E.c6()/864]
expected=[[k*k*(x+E.b2()/12),k**3*(y+(ai[0]*x+ai[2])/2)] for x,y in original]
maps=data['map_coefficients']
pairs=[[R(v) for v in pair] for pair in maps['short_coefficient_pairs']]
triples=[[R(v) for v in P] for P in maps['projective_points']]
t_pair=[R(v) for v in maps['base_parameter_pair']]
assert len(pairs)==2 and len(triples)==18
for ps in pairs+triples+[t_pair]:
    assert all(v in ZZ for p in ps for v in p)
    assert gcd(ps).degree()==0
    assert gcd([ZZ(v) for p in ps for v in p])==1
for pair,v in zip(pairs,AB): assert pair[0]==v*pair[1]
for ps,(x,y) in zip(triples,expected):
    assert ps[0]==x*ps[2] and ps[1]==y*ps[2]
    assert y*y==x**3+AB[0]*x+AB[1]
T=dec(source['base_map']); assert t_pair[0]==T*t_pair[1]
assert max(T.numerator().degree(),T.denominator().degree())==2
r=QQ(data['dependent_conic_slope']); a,b=r.numerator(),r.denominator()
assert T(r)==QQ(data['dependent_original_parameter'])
Ed=EllipticCurve(QQ,[v(r) for v in AB])
Pd=[Ed([x(r),y(r)]) for x,y in expected]
assert list(map(str,Ed.a_invariants()))==frame['curve']==data['control_anchor_equation']
assert [list(map(str,P.xy())) for P in Pd[:17]]==frame['basis']
relation=data['dependence']
assert ZZ(relation['multiplier'])*Pd[-1]==sum((ZZ(w)*Q for w,Q in zip(relation['word'],Pd[:17])),Ed(0))
assert relation['multiplier']>0

L,P,q0,period=map(ZZ,[data['L'],data['P'],data['q0'],data['period']])
S=data['control_primes']; proof=data['proof_primes']
assert S==sorted(row['prime'] for row in frame['records']) and len(S)==25
assert set(S).isdisjoint(proof) and all(179<p<=1009 and ZZ(p).is_prime() for p in proof)
assert len(set(proof))==19 and P==prod(ZZ(p) for p in proof)
assert L==prod(ZZ(row['p'])**row['exponent'] for row in data['control_stability'])
assert period==L*P and 0<q0<period and gcd(L,P)==1
assert q0%L==1 and (a*q0+b*L)%P==0 and gcd(b*q0,P)==1
# Homogeneous u_n=[a(q0+period*n)+bL : b(q0+period*n)].
U0=a*q0+b*L; U1=a*period; V0=b*q0; V1=b*period
assert (U0-a)%L==(V0-b)%L==U1%L==V1%L==0
assert U0%P==U1%P==V1%P==0 and gcd(V0,P)==1
assert U1*V0-U0*V1==-b*b*L*period and U1*V0-U0*V1!=0

def homogeneous(ps,x,z):
    d=max(p.degree() for p in ps)
    return [ZZ(sum(c*x**j*z**(d-j) for j,c in enumerate(p))) for p in ps]
def residue(v,p):
    v=QQ(v); assert v.denominator()%p
    return int(v.numerator()*v.denominator().inverse_mod(p)%p)
def point_reduction(ps,x,z,p):
    values=homogeneous(ps,x,z)
    common=min(v.valuation(p) for v in values if v)
    values=[residue(QQ(v)/ZZ(p)**common,p) for v in values]
    assert any(values)
    if values[2]==0:
        assert values[0]==0 and values[1]
        return None
    inv=pow(values[2],-1,p)
    return (values[0]*inv%p,values[1]*inv%p)

def finite_group(p,A,B):
    roots={}
    for y in range(p): roots.setdefault(y*y%p,[]).append(y)
    points=[None]+[(x,y) for x in range(p) for y in roots.get((x**3+A*x+B)%p,[])]
    pool=set(points)
    def add(P,Q):
        if P is None: return Q
        if Q is None: return P
        x,y=P; xx,yy=Q
        if x==xx and (y+yy)%p==0: return None
        slope=((3*x*x+A)*pow(2*y,-1,p) if x==xx else (yy-y)*pow(xx-x,-1,p))%p
        xxx=(slope*slope-x-xx)%p
        result=(xxx,(slope*(x-xxx)-y)%p); assert result in pool
        return result
    doubles={add(Q,Q) for Q in points}
    labels={Q:0 for Q in doubles}; reps=[None]; dim=0
    for Q in points:
        if Q in labels: continue
        additions=[]
        for code,V in enumerate(reps):
            W=add(V,Q); additions.append(W)
            for D in doubles:
                Z=add(W,D); assert Z not in labels
                labels[Z]=code|(1<<dim)
        reps+=additions; dim+=1
    assert set(labels)==pool and len(points)==len(doubles)*2**dim
    return points,doubles,labels,reps,dim

control_cert=[]; control_rows=[]
for p,row in zip(S,data['control_stability']):
    assert row['p']==p and L.valuation(p)==row['exponent']
    values=[homogeneous(pair,a,b) for pair in pairs]
    den=[int(v[1].valuation(p)) for v in values]
    commons=[int(min(v.valuation(p) for v in homogeneous(Q,a,b) if v)) for Q in triples]
    par=min(v.valuation(p) for v in homogeneous(t_pair,a,b) if v)
    assert den==row['coefficient_denominator_orders']
    assert commons==row['point_common_orders'] and par==row['parameter_common_order']
    assert all(row['exponent']>v for v in den+commons+[par])
    assert Ed.discriminant().valuation(p)==0
    A0,B0=[residue(v(r),p) for v in AB]
    allpoints,doubles,labels,reps,dim=finite_group(p,A0,B0)
    reductions=[point_reduction(Q,a,b,p) for Q in triples]
    assert all(Q in labels for Q in reductions)
    codes=[labels[Q] for Q in reductions]
    block=[[(v>>j)&1 for v in codes] for j in range(dim)]; control_rows+=block
    control_cert.append({'p':p,'group_order':len(allpoints),'dimension':dim,
                         'marked_point_reductions':reductions,'matrix_rows':block,
                         'preserved_for_all_integer_n':True})
M=matrix(GF(2),len(control_rows),18,[x for row in control_rows for x in row])
assert M.rank()==17 and M[:,:17].rank()==17

proof_cert=[]; proof_rows=[]
E0=EllipticCurve(QQ,[v(0) for v in AB])
for p in proof:
    assert all(pair[1](0)%p for pair in pairs+[t_pair])
    assert all(any(v(0)%p for v in Q) for Q in triples)
    assert E0.discriminant().valuation(p)==0
    A0,B0=[residue(v(0),p) for v in AB]
    allpoints,doubles,labels,reps,dim=finite_group(p,A0,B0)
    reductions=[point_reduction(Q,ZZ(0),ZZ(1),p) for Q in triples]
    assert all(Q in labels for Q in reductions)
    codes=[labels[Q] for Q in reductions]
    block=[[(v>>j)&1 for v in codes] for j in range(dim)];proof_rows+=block
    if p==data['no_two_torsion_prime']: assert len(allpoints)%2==1
    proof_cert.append({'p':p,'curve_coefficients':[A0,B0],'all_points':allpoints,
                       'doubled_points':sorted(doubles,key=lambda Q:(Q is not None,Q)),
                       'coset_representatives':reps,'marked_point_reductions':reductions,
                       'matrix_rows':block,'preserved_for_all_integer_n':True})
M=matrix(GF(2),len(proof_rows),18,[x for row in proof_rows for x in row])
assert M.rank()==18 and M[:,:17].rank()==17
assert data['no_two_torsion_prime'] in proof

result={
    'status':'PASS_INDEPENDENT_UNIFORM_LOCAL_SHADOW_RANK18',
    'classification':'new constructive theorem, independently replayed',
    'sequence':data['sequence'],'r':str(r),'L':str(L),'P':str(P),'q0':str(q0),
    'all_integer_terms_defined_and_smooth':True,
    'all_integer_terms_certified_rank_lower_bound':18,
    'control_footprint_rank_on_18_columns':17,
    'disjoint_proof_footprint_rank_on_18_columns':18,
    'limiting_branch_is_inherited':True,
    'limiting_elliptic_fibre_rank':'NOT_BOUNDED_ABOVE',
    'control_certificates':control_cert,'proof_certificates':proof_cert,
    'uniformity':'Homogeneous parameter pairs are congruent to[a:b] at the required control powers and to[0:unit] at proof primes. Primitive map orders preserve every reduction. Nonzero proof reductions rule out poles, singular curves and undefined point triples for every integer n. The nonconstant Mobius parameter sequence converges in R to r; the degree2 original base map gives infinitely many distinct elliptic parameters.',
    'scope':'The fixed25-prime observations cannot distinguish these independent branches from the dependent branch. No assertion about all local places, finite data chosen adaptively, generic Selmer dimensions, or302 member selection.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',OUT/'construction.json',Path(__file__)]}}
result=json.loads(json.dumps(result))  # Normalize finite-point tuples to JSON arrays.
destination=OUT/'independent-replay-v2.json'
if destination.exists(): assert read(destination)==result
else:
    with destination.open('x') as stream: json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
print(result['status'],'local rank17, proof rank18,25+19 primes',flush=True)
