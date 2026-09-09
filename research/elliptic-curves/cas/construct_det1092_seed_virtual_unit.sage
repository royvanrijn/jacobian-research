#!/usr/bin/env sage-python
"""Make the seed's even-valuation representative; exact cubic LLL only.

One normalized class, one three-dimensional LLL, first vector only. No
unit/class group, principal-ideal solver, factorization, or point search.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,NumberField,pari,matrix,lcm
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_seed_half_ideal_v1'
ARITH=ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
SEED=OUT/'302-first-unlock.json';PROJECTION=OUT/'generic-projection.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists():assert json.loads(p.read_text())==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def rows(A):return [[str(x) for x in row] for row in A.rows()]
def coeff(a):return [str(a[i]) for i in range(3)]
save('virtual-unit-protocol.json',{
    'classification':'retrospective construction using the fixed generic parity solution',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [ARITH,SEED,PROJECTION,Path(__file__)]},
    'limits':{'seconds':25,'normalized_classes':1,'LLL_dimension':3,'LLL_vectors_tested':1,
              'class_groups':0,'unit_groups':0,'point_searches':0},
    'selection':'First vector of exact trace-LLL of the inverse of the certified square-root ideal. Not a shortest-vector assertion.'})
d=json.loads(ARITH.read_text());s=json.loads(SEED.read_text());p=json.loads(PROJECTION.read_text())
word=p['cancelling_generic_word'];R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending'])
K=NumberField(f,'theta');nf=pari.nfinit([pari(f),d['S_finite']])
assert ZZ(nf.disc())==ZZ(d['field_discriminant']) and ZZ(nf[3])==ZZ(d['defining_order_index'])
zk=[K(R(str(v))) for v in nf.nf_get_zk()]
def pa(a):return pari.Mod(pari(R(a.list())),pari(f))
alpha=K.one();I=pari.idealhnf(nf,1);factors=[K(R(s['alpha']))]
factors += [K(R(x['beta_ascending'])) for bit,x in zip(word,d['generic_classes']) if bit]
for a in factors:
    norm=ZZ(a.norm());assert norm.is_square()
    alpha*=a;I=pari.idealmul(nf,I,pari.idealhnf(nf,pa(a),pari(norm.sqrt())))
J=I;corrections=[]
for prime in d['S_finite']:
    for index,P in enumerate(pari.idealprimedec(nf,prime)):
        v=int(pari.idealval(nf,pa(alpha),P));w=int(pari.idealval(nf,I,P));assert v%2==0
        e=v//2-w
        if e:
            J=pari.idealmul(nf,J,pari.idealpow(nf,P,e));corrections.append([prime,index,e])
assert pari.idealpow(nf,J,2)==pari.idealhnf(nf,pa(alpha))
inv=matrix(QQ,pari.idealinv(nf,J).sage())
basis=[sum((inv[i,j]*zk[i] for i in range(3)),K.zero()) for j in range(3)]
G=matrix(QQ,3,3,lambda i,j:(basis[i]*basis[j]).trace())
assert all(G[:i,:i].det()>0 for i in [1,2,3])
scale=lcm([x.denominator() for x in G.list()])
U=matrix(ZZ,pari.qflllgram(pari(matrix(ZZ,scale*G))).sage());assert abs(U.det())==1
gamma=sum((U[i,0]*basis[i] for i in range(3)),K.zero())
beta=alpha*gamma**2
Jred=pari.idealmul(nf,J,pa(gamma))
assert pari.idealpow(nf,Jred,2)==pari.idealhnf(nf,pa(beta))
assert matrix(QQ,Jred.sage()).denominator()==1
N=QQ(beta.norm());assert N.is_square()
answer={'status':'PASS_EXPLICIT_SEED_VIRTUAL_UNIT',
    'classification':'retrospective exact construction and new deduction',
    'checker_sha256':sha(Path(__file__)),
    'generic_word':word,'alpha':coeff(alpha),'gamma':coeff(gamma),'beta':coeff(beta),
    'square_root_ideal':rows(matrix(QQ,J.sage())),
    'reduced_square_root_ideal':rows(matrix(QQ,Jred.sage())),
    'inverse_ideal_basis':rows(inv),'trace_Gram':rows(G),'LLL_transform':rows(U),
    'bad_prime_corrections':corrections,
    'reduced_ideal_norm':str(pari.idealnorm(nf,Jred)),
    'norm_square_root':str(N.sqrt()),'unit_found':pari.idealnorm(nf,Jred)==1,
    'identity':'beta = alpha_seed * product(alpha_generic_i^word_i) * gamma^2; (beta)=Jred^2',
    'squareclass_conclusion':'beta lies outside the generic squareclass image, and all its finite prime-ideal valuations are even.',
    'class_group_lower_bound':7,
    'class_group_argument':'Generic17 and seed span18 independent norm-square classes. Their ideal-parity map has rank9. The9-dimensional even-valuation kernel maps to Cl(K)[2], with kernel at most2 from norm-positive units in a totally real cubic.',
    'boundary':'Virtual unit means even valuations, not an actual unit. The individual ideal class may still be trivial or generic. No prospective way to select this seed-derived beta has been established.'}
save('virtual-unit.json',answer)
print(answer['status'],'reduced ideal norm',answer['reduced_ideal_norm'],'unit',answer['unit_found'],flush=True)
print('beta coefficient numerator digits',[len(str(abs(QQ(c).numerator()))) for c in answer['beta']],flush=True)
