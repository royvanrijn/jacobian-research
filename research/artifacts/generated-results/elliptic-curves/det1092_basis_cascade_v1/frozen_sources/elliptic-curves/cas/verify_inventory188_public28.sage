#!/usr/bin/env python3
"""Independent Sage group enumeration of E(Fp)/2 for the public28 reproduction."""
import json,hashlib
from pathlib import Path
from sage.all import QQ,ZZ,GF,EllipticCurve,Matrix,lcm,gcd,PolynomialRing

ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'inventory188_public28_reproduction_v1.json'
OUT=ART/'inventory188_public28_sage_replay_v1.json'
if OUT.exists():raise FileExistsError('preserve independent group replay')
d=json.loads(INPUT.read_bytes())
for n,h in d['sources'].items():
    if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h:raise ArithmeticError('bound public evidence changed')
E=EllipticCurve(QQ,d['curve']);C=EllipticCurve(QQ,d['public_curve']);u=QQ(d['transport_scale_u'])
points=[E([QQ(x),QQ(y)]) for x,y in d['points']]
for raw,expected in zip(d['public_points'],d['transported_public_points']):
    P=C([QQ(x) for x in raw]);x=(P[0]+C.b2()/12)/u**2;y=(P[1]+(C.a1()*P[0]+C.a3())/2)/u**3
    if E([x,y])!=E([QQ(c) for c in expected]):raise ArithmeticError('public transport differs')
if len(points)!=55 or points[27:]!=[E([QQ(x),QQ(y)]) for x,y in d['transported_public_points']]:
    raise ArithmeticError('complete public/local point roster required')
projective=[]
for P in points:
    den=lcm([c.denominator() for c in P]);v=[ZZ(c*den) for c in P];g=gcd(v)
    projective.append([c//g for c in v])
rows=[];groups=[]
for prime in [s['prime'] for s in d['union_signatures']]:
    F=GF(prime);e=EllipticCurve(F,[F(c) for c in E.a_invariants()])
    if e.discriminant()==0:raise ArithmeticError('bad reduction in proof')
    elements=e.points();key=lambda P:tuple(int(c) for c in P)
    doubles={key(2*P):2*P for P in elements};reps=[e(0)]
    masks={key(P):0 for P in doubles.values()}
    for P in elements:
        if key(P) in masks:continue
        bit=len(reps);extra=[R+P for R in reps]
        for i,R in enumerate(extra):
            for T in doubles.values():masks[key(R+T)]=i+bit
        reps+=extra
    if len(masks)!=len(elements) or len(reps) not in (1,2,4):raise ArithmeticError('quotient enumeration incomplete')
    dimension=ZZ(len(reps)).valuation(2)
    reduced=[e([F(c) for c in P]) for P in projective]
    rows.extend([[(masks[key(P)]>>i)&1 for P in reduced] for i in range(dimension)])
    groups.append({'prime':prime,'order':len(elements),'double_subgroup_order':len(doubles),'quotient_dimension':int(dimension)})
M=Matrix(GF(2),rows);chosen=d['independent_column_indices']
union_rank=int(M.rank());chosen_rank=int(M.matrix_from_columns(chosen).rank());public_rank=int(M.matrix_from_columns(list(range(27,55))).rank())
if union_rank!=d['rank_lower_bound'] or chosen_rank!=len(chosen) or public_rank!=28:
    raise ArithmeticError('independent enumerated quotient ranks differ')
p=d['rank_certificate']['no_rational_2_torsion_prime'];F=GF(p);R=PolynomialRing(F,'x');x=R.gen()
if not (x**3+F(E.a4())*x+F(E.a6())).is_irreducible():raise ArithmeticError('rational2-torsion exclusion failed')
result={'schema':'elliptic-curves.inventory188-public28-sage.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve(),INPUT)},
        'groups':groups,'union_finite_rank':union_rank,'exported_basis_rank':chosen_rank,'public_basis_rank':public_rank,
        'no_two_torsion_prime':p,'claim_boundary':'Independent Sage point membership and exact source transport, complete finite group enumeration, quotient-by-doubles matrices and irreducible2-division cubic. This proves the public28 and union lower bounds, not an upper rank, new discovery or prospective point-search recovery.'}
with OUT.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print('INDEPENDENT SAGE PUBLIC28 AND UNION',union_rank,'PASS')
