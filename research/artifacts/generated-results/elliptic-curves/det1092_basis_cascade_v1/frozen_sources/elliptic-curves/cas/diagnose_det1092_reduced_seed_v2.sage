#!/usr/bin/env sage-python
"""Bounded exact finite-kernel and odd-prime check of the frozen unit seed."""
import json,sys
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,EllipticCurve,matrix
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
import audit_recorded_point_mod2_rank_v3 as m2
import audit_retained_cloud_modl as ml
from research_runtime.memory_store import MemoryFactStore
from research_runtime.finite_reduction import ReductionCache
from research_runtime.store import checkpoint
D=ROOT/'artifacts/local/elliptic-curves/det1092-reduced-seed-diagnosis-v2'
p=json.loads((D/'protocol.json').read_bytes());assert p['script_sha256']==cert.hashed(Path(__file__))
source=ROOT/p['parent'];assert cert.hashed(source)==p['parent_sha256'];data=cert.read(source)
R=PolynomialRing(QQ,'t');val=lambda v: R(v['numerator'])(1)/R(v['denominator'])(1)
a=[val(v) for v in data['a_invariants']];E=EllipticCurve(QQ,a);points=[E([val(v) for v in P]) for P in data['basis_weierstrass_coordinates']];model=tuple(cert.F(str(v)) for v in a);xy=[tuple(cert.F(str(v)) for v in P.xy()) for P in points]
rows=[];cache=ReductionCache(MemoryFactStore());sigs=[]
for q in m2._primes_up_to(997):
 if q==2:continue
 try:sg=m2.signature(cache,model,xy,q)
 except ValueError:continue
 rows.extend(sg.rows)
assert rows and all(len(r)==17 for r in rows)
M=matrix(GF(2),rows);ker=M.right_kernel().basis();out=dict(status='RUNNING',finite_rank=int(M.rank()),kernel_words=[list(map(int,v)) for v in ker],relations=[]);checkpoint(D/'result.json',out)
for v in ker:
 word=list(map(int,v));S=sum((word[i]*points[i] for i in range(17)),E(0));row=dict(word=word,is_zero=bool(S.is_zero()))
 if not S.is_zero():row['sum']=list(map(str,S.xy()))
 out['relations'].append(row);checkpoint(D/'result.json',out)
 print('KERNEL',word,'ZERO',S.is_zero(),flush=True)
inp=dict(status='COMPLETE_DECLARED_FINITE_AUDIT',family='det1092-reduced-chart',parameter='1',curve=list(map(str,a)),points=[list(map(str,P.xy())) for P in points],rank_lower_bound=int(M.rank()))
checkpoint(D/'input.json',inp);ml.build(D/'input.json',D/'odd.json');ml.check(D/'odd.json');out['status']='PASS';out['odd_ranks']=[r['finite_column_rank'] for r in cert.read(D/'odd.json')['audits']];checkpoint(D/'result.json',out)
