#!/usr/bin/env sage-python
"""Numerical proposals followed by exact rational group equality, with fixed bounds."""
import json,sys,hashlib,math
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from pathlib import Path
from sage.all import QQ,EllipticCurve,matrix,vector,RealField
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
from research_runtime.store import checkpoint
D=ROOT/'artifacts/local/elliptic-curves/det1092-unit-section-relations-v2'
p=json.loads((D/'protocol.json').read_bytes());assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==p['script_sha256']
source=ROOT/p['seed'];assert hashlib.sha256(source.read_bytes()).hexdigest()==p['seed_sha256'];seed=json.loads(source.read_bytes())
E=EllipticCurve(QQ,seed['curve']);points=[E(P) for P in seed['points']];indices=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,16];basis=[points[i] for i in indices]
geometry=SourceFileLoader('section_relation_metric',str(CAS/'prospective_half_lattice_v3.sage')).load_module()
gram,asym=geometry.canonical_height_gram(tuple(map(Fraction,seed['curve'])),[tuple(map(Fraction,P)) for P in seed['points']]);real=RealField(384)
G=matrix(real,[[str(gram[i][j]) for j in indices] for i in indices]);relations=[]
for j in [14,15]:
 v=G.solve_right(vector(real,[str(gram[i][j]) for i in indices]));q=[Fraction(str(x)).limit_denominator(p['maximum_denominator']) for x in v];den=math.lcm(*(x.denominator for x in q));word=[int(x*den) for x in q]
 assert den<=p['maximum_denominator'] and all(abs(x)<=p['maximum_coefficient'] for x in word)
 assert den*points[j]==sum((word[i]*basis[i] for i in range(15)),E(0))
 relations.append(dict(section_index=j,denominator=den,word=word));print('EXACT RELATION',j,den,word,flush=True)
out=dict(status='PASS',curve=seed['curve'],points=seed['points'],relations=relations,independent_indices=indices,section_span_rank=15,maximum_gram_asymmetry=str(asym),boundary='The two omitted section values are in the rational span of the chosen fifteen, by exact Sage group identities. Rank15 below is certified separately by finite reductions. Numerical heights only propose identities; no whole-fibre rank upper bound or new point search.')
assert not (D/'result.json').exists();checkpoint(D/'result.json',out)
