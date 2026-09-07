#!/usr/bin/env sage-python
"""Exact specialized subgroups; bounded relation proposals verified as group identities."""
import sys,math
from pathlib import Path
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from sage.all import QQ,PolynomialRing,EllipticCurve,RealField,matrix,vector
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
from research_runtime.search_state import raw_state
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.store import checkpoint,digest
D=ROOT/'artifacts/local/elliptic-curves/mixed-reduced-parent-exposure-v1'
def main():
 p=cert.read(D/'prepare-protocol.json');assert all(cert.hashed(ROOT/n)==h for n,h in p['sources'].items());data=cert.read(D/'intake-input.json');assert cert.hashed(D/'intake-input.json')==p['input_sha256']
 R=PolynomialRing(QQ,'s');K=R.fraction_field();parent=data['parent'];generic_a=[K(R(v['numerator'])/R(v['denominator'])) for v in parent['a_invariants']];generic_P=[[K(R(v['numerator'])/R(v['denominator'])) for v in P] for P in parent['basis_weierstrass_coordinates']]
 geometry=SourceFileLoader('mixed_seed_metric',str(CAS/'prospective_half_lattice_v3.sage')).load_module();rows=[]
 for row in data['rows']:
  new=row['parent']=='det1092'
  if new:
   t=QQ(row['parameter']);a=[v(t) for v in generic_a];E=EllipticCurve(QQ,a);original=[]
   for x,y in generic_P:
    if x.denominator()(t)==0 or y.denominator()(t)==0:
     assert x.denominator()(t)==y.denominator()(t)==0;original.append(E(0))
    else:original.append(E([x(t),y(t)]))
   candidates=[P for P in original if not P.is_zero()]
  else:
   old=data['retained'][row['id']];a=list(map(QQ,old['curve']));E=EllipticCurve(QQ,a);candidates=[E(P) for P in old['points']];original=candidates
  model=tuple(F(str(v)) for v in a);xy=[tuple(F(str(v)) for v in P.xy()) for P in candidates]
  state=raw_state(model,xy,cache=QuotientOnlyReductionCache(MemoryFactStore()),prime_bound=1000);basis=[E(P) for P in state.basis];r=len(basis);assert 12<=r<=len(candidates)
  if not new:assert r==26 and basis==candidates
  proof=checked_rank(model,[tuple(F(str(v)) for v in P.xy()) for P in basis],state.reductions.primes,state.no_two_torsion_prime)
  relations=[]
  if new:
   missing=[(j,P) for j,P in enumerate(original) if not P.is_zero() and P not in basis]
   if missing:
    pts=[*basis,*[P for j,P in missing]];gram,asym=geometry.canonical_height_gram(model,[tuple(F(str(v)) for v in P.xy()) for P in pts]);real=RealField(384);G=matrix(real,[[str(gram[i][j]) for j in range(r)] for i in range(r)])
   for j,P in enumerate(original):
    if P.is_zero():den=1;word=[0]*r
    elif P in basis:den=1;word=[int(Q==P) for Q in basis]
    else:
     k=next(i for i,(ix,Q) in enumerate(missing) if ix==j);coords=G.solve_right(vector(real,[str(gram[i][r+k]) for i in range(r)]));fractions=[F(str(v)).limit_denominator(p['maximum_denominator']) for v in coords];den=math.lcm(*(v.denominator for v in fractions));word=[int(v*den) for v in fractions]
    assert den<=p['maximum_denominator'] and max(map(abs,word))<=p['maximum_coefficient'];assert den*P==sum((word[i]*basis[i] for i in range(r)),E(0))
    relations.append(dict(section_index=j,denominator=den,word=word))
   generic=[list(map(str,P.xy())) for P in basis]
  else:generic=old['generic_points']
  points=[list(map(str,P.xy())) for P in basis];assert points[:len(generic)]==generic
  seed=dict(family=row['family'],parameter=row['parameter'],curve=list(map(str,a)),points=points,generic_points=generic,rank_certificate=proof,source_section_points=[None if P.is_zero() else list(map(str,P.xy())) for P in original],section_relations=relations,initial_rank=r,model_coefficient_bits=max(max(abs(v.numerator()).nbits(),v.denominator().nbits()) for v in a))
  out=D/row['id']/'seed.json';assert not out.exists();checkpoint(out,seed);rows.append(dict(row,initial_rank=r,generic_dimension=len(generic),mask_floor=0 if new else len(generic),seed_sha256=cert.hashed(out)));checkpoint(D/'intake-result.json',dict(status='RUNNING',rows=rows));print('EXACT SEED',row['id'],r,'from',len(original),flush=True)
 checkpoint(D/'intake-result.json',dict(status='PASS',rows=rows,prepare_protocol_hash=digest(p)))
if __name__=='__main__':main()
