#!/usr/bin/env sage-python
"""Project the first seed's bad-prime parity against generic17, <=25s.

All places and generic classes were frozen before this diagnostic. No new
point, class group, prime search, or family parameter is introduced.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,pari,matrix,vector
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_seed_half_ideal_v1'
ARITH=ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
SEED=OUT/'302-first-unlock.json'
PARITY=OUT/'ideal-parity.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists():assert json.loads(p.read_text())==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
save('generic-projection-protocol.json',{
    'classification':'retrospective evaluation using frozen generic classes and places',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [ARITH,SEED,PARITY,Path(__file__)]},
    'limits':{'seconds':25,'generic_classes':17,'exceptional_classes':1,
              'rational_primes':20,'class_groups':0,'point_searches':0},
    'selection':'Existing generic17 order, fixed ascending bad primes, canonical F2 linear solve. No place tuning.'})
d=json.loads(ARITH.read_text());s=json.loads(SEED.read_text())
R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending']);F=GF(2)
nf=pari.nfinit([pari(f),d['S_finite']]);assert ZZ(nf[3])==ZZ(d['defining_order_index'])
places=[(p,i,P) for p in d['S_finite'] for i,P in enumerate(pari.idealprimedec(nf,p))]
polys=[R(x['beta_ascending']) for x in d['generic_classes']]+[R(s['alpha'])]
valuations=[[int(pari.idealval(nf,pari.Mod(pari(a),pari(f)),P)) for p,i,P in places] for a in polys]
A=matrix(F,valuations[:17]);v=vector(F,valuations[17]);aug=A.stack(matrix(F,[v]))
in_image=A.rank()==aug.rank()
answer={'status':'PASS_GENERIC_IDEAL_PARITY_PROJECTION',
    'classification':'verified application','checker_sha256':sha(Path(__file__)),
    'places':[[p,i,int(P[2]),int(P[3])] for p,i,P in places],
    'generic_valuations':valuations[:17],'seed_valuations':valuations[17],
    'generic_rank':int(A.rank()),'augmented_rank':int(aug.rank()),
    'generic_kernel_basis':[list(map(int,row)) for row in A.left_kernel().basis()],
    'seed_in_generic_ideal_parity_image':in_image}
if in_image:
    word=A.transpose().solve_right(v);assert word*A==v
    answer['cancelling_generic_word']=list(map(int,word))
    answer['interpretation']='After multiplication by this generic word, the seed defines a globally even-valuation squareclass. Its possible ideal-class or unit origin is not determined.'
else:
    dual=aug.solve_right(vector(F,[0]*17+[1]));assert A*dual==0 and v*dual==1
    answer['separating_functional']=list(map(int,dual))
    answer['separating_places']=[answer['places'][i] for i,b in enumerate(dual) if b]
    answer['interpretation']='No generic Kummer multiplier can turn the seed into an everywhere-even-valuation squareclass. It is already distinguished in bad-prime ideal parity, before an ordinary class-group projection.'
save('generic-projection.json',answer)
print('generic rank',A.rank(),'augmented',aug.rank(),'seed in generic parity',in_image,flush=True)
print('word or places',answer.get('cancelling_generic_word',answer.get('separating_places')),flush=True)
