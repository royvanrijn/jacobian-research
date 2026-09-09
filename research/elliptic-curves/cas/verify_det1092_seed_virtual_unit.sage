#!/usr/bin/env sage-python
"""Independent rational-lattice replay; no PARI number-field/ideal calls.

Uses the previously certified maximal order and prime ideals, but checks
valuations by integral lattice membership, products by cubic multiplication
and integer HNF, and squareclass independence by finite-field characters.
Run with timeout 25s.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,lcm
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_seed_half_ideal_v1'
ARITH=ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
INPUTS=[ARITH,OUT/'302-first-unlock.json',OUT/'generic-projection.json',OUT/'virtual-unit.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
d,s,projection,result=[json.loads(p.read_text()) for p in INPUTS]
R=PolynomialRing(QQ,'z');z=R.gen();f=R(d['cubic_ascending'])
assert f.change_ring(GF(31)).is_irreducible() and f.discriminant()>0
zk=[R(v) for v in d['maximal_order_basis']]
def vec(a):return vector(QQ,[a[i] for i in range(3)])
Z=matrix(QQ,[vec(a) for a in zk]).transpose();Zi=Z.inverse()
def mul(a,b):return (a*b)%f
def algebra(v):return sum((v[i]*zk[i] for i in range(3)),R.zero())
def canonical(A):
    scale=lcm([q.denominator() for q in A.list()])
    return matrix(ZZ,(scale*A).transpose()).hermite_form(include_zero_rows=False).transpose()/scale
def ideal_product(A,B):
    cols=[Zi*vec(mul(algebra(a),algebra(b))) for a in A.columns() for b in B.columns()]
    return canonical(matrix(QQ,cols).transpose())
def principal(a):return canonical(matrix(QQ,[Zi*vec(mul(a,b)) for b in zk]).transpose())
def contains(A,v):return all(q.denominator()==1 for q in A.solve_right(v))
def norm(a):return matrix(QQ,[vec(mul(a,z**i)) for i in range(3)]).transpose().det()
def parse_hnf(text):return matrix(QQ,[[QQ(x) for x in row.split(',')] for row in text[1:-1].split(';')])
O=matrix.identity(QQ,3)
assert ideal_product(O,O)==O
assert abs(Z.det())==1/QQ(d['defining_order_index'])
assert f.discriminant()*Z.det()**2==QQ(d['field_discriminant'])
polys=[R(x['beta_ascending']) for x in d['generic_classes']]+[R(s['alpha'])]
valuations=projection['generic_valuations']+[projection['seed_valuations']]
places=[]
for block in d['prime_decomposition']:
    p=block['p'];product=O
    for index,row in enumerate(block['primes']):
        P=parse_hnf(row['hnf']);e=row['ramification_index'];degree=row['residue_degree']
        assert abs(P.det())==QQ(p)**degree
        for _ in range(e):product=ideal_product(product,P)
        places.append((p,index,e,degree,P))
    assert product==p*O
assert [list(row[:4]) for row in places]==projection['places']
for j,(p,index,e,degree,P) in enumerate(places):
    maxpower=max(row[j] for row in valuations)+1;powers=[O]
    for _ in range(maxpower):powers.append(ideal_product(powers[-1],P))
    for a,row in zip(polys,valuations):
        v=row[j];coordinates=Zi*vec(a)
        assert contains(powers[v],coordinates) and not contains(powers[v+1],coordinates)
# Independently show that the factor-free half-ideal defect is supported at S.
for a in polys:
    N=ZZ(norm(a));assert N.is_square() and N>0;y=N.sqrt()
    I=canonical(principal(a).augment(y*O));Isq=ideal_product(I,I)
    ainv=a.inverse_mod(f)
    defect=matrix(QQ,[Zi*vec(mul(ainv,algebra(c))) for c in Isq.columns()]).transpose()
    assert defect.denominator()==1
    remainder=abs(I.det())**2/N
    for p in d['S_finite']:remainder/=QQ(p)**remainder.valuation(p)
    assert remainder==1
# Frozen good-prime pool only. An unusable prime is omitted, never replaced.
chars=[[] for _ in polys];used=[];omitted=[]
for block in d['independence_blocks']:
    p=block['prime'];roots=block['roots'];F=GF(p)
    assert len(roots)==3 and all(f.change_ring(F)(r)==0 for r in roots)
    values=[[a.change_ring(F)(r) for r in roots] for a in polys]
    if any(x==0 for row in values for x in row):omitted.append(p);continue
    for dest,row in zip(chars,values):dest.extend(int(not x.is_square()) for x in row)
    used.append(p)
characters=matrix(GF(2),chars);assert characters[:17,:].rank()==17 and characters.rank()==18
A=matrix(GF(2),valuations[:17]);full=matrix(GF(2),valuations)
assert A.rank()==full.rank()==9
word=vector(GF(2),result['generic_word']);assert word*A==vector(GF(2),valuations[17])
alpha=polys[-1]
for bit,a in zip(word,polys[:17]):
    if bit:alpha=mul(alpha,a)
assert alpha==R(result['alpha'])
beta=R(result['beta']);gamma=R(result['gamma'])
assert beta==mul(alpha,mul(gamma,gamma))
J=matrix(QQ,result['square_root_ideal']);Jred=matrix(QQ,result['reduced_square_root_ideal'])
assert ideal_product(J,J)==principal(alpha)
assert ideal_product(Jred,Jred)==principal(beta)
assert ideal_product(J,principal(gamma))==canonical(Jred)
assert Jred.denominator()==1 and norm(beta)==QQ(result['norm_square_root'])**2
assert abs(Jred.det())==QQ(result['reduced_ideal_norm'])
report={'status':'PASS_INDEPENDENT_CUBIC_LATTICE_AND_CHARACTER_REPLAY',
    'classification':'verified application and new lower-bound deduction',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [*INPUTS,Path(__file__)]},
    'generic_squareclass_rank':17,'with_seed_squareclass_rank':18,
    'ideal_parity_rank':9,'generic_even_valuation_dimension':8,
    'with_seed_even_valuation_dimension':9,
    'norm_positive_unit_squareclass_dimension':2,'class_group_2_rank_lower_bound':7,
    'good_character_primes_used':used,'good_character_primes_omitted':omitted,
    'prime_ideal_valuation_checks':len(places)*len(polys),
    'square_ideal_identities':2,'squareclass_identity':True,
    'boundary':'No class-group order or generators, no individual nonprincipal-ideal proof, no new strict class, no Selmer/Sha dimension, and no prospective seed selection.'}
path=OUT/'independent-replay.json'
if path.exists():assert json.loads(path.read_text())==report
else:
    with path.open('x') as stream:json.dump(report,stream,indent=2,sort_keys=True);stream.write('\n')
print(report['status'],'valuation checks',report['prime_ideal_valuation_checks'],
      'class2 rank >=',report['class_group_2_rank_lower_bound'],flush=True)
