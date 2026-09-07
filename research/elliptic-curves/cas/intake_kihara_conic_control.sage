#!/usr/bin/env sage-python
"""One balanced rational conic control; no score or point enumeration."""
import json,sys,hashlib
from dataclasses import asdict
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,gcd,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/kihara-conic-control-v1';OUT=ART/'kihara_conic_control_v1.json';sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from audit_recorded_point_mod2_rank_v3 import signature,insert,_primes_up_to
from memory_rank_certificate import checked_rank

def main():
 protocol=json.loads((D/'protocol.json').read_text())
 for n,h in protocol['sources'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h
 data=json.loads((ART/'kihara_split_infinity_v1.json').read_text());c=QQ(data['conic_constant']);d=c.squarefree_part();r=(c/d).sqrt();v=ZZ(d).isqrt()+1;assert d>0 and (v-1)**2<d<v*v;s0=r*v
 R=PolynomialRing(QQ,'s');s=R.gen();K=R.fraction_field();A=R(data['A']);B=R(data['B']);scale=K(data['integral_function_field_scale'])(s0)
 # First undo the harmless construction denominator scale, then clear only
 # the necessary denominator valuations. Remove common powers only at the
 # same known denominator primes and2,3; no discriminant factorization.
 av=A(s0)/scale**4;bv=B(s0)/scale**6;primes=set(ZZ(av.denominator()).prime_divisors()+ZZ(bv.denominator()).prime_divisors()+[2,3]);u=QQ(1)
 for p in primes:u*=QQ(p)**max((-av.valuation(p)+3)//4,(-bv.valuation(p)+5)//6)
 assert (av*u**4).denominator()==(bv*u**6).denominator()==1
 E=EllipticCurve(QQ,[av*u**4,bv*u**6]);transport=u/scale
 pts=[E([transport**2*K(data['sections'][i][0])(s0),transport**3*K(data['sections'][i][1])(s0)]) for i in data['basis_indices']]
 model=tuple(cert.F(str(c)) for c in E.a_invariants());cloud=[tuple(cert.F(str(c)) for c in p.xy()) for p in pts]
 cache=ReductionCache(MemoryFactStore());pivots={};signatures=[];torsion=None
 for prime in _primes_up_to(997):
  if prime==2:continue
  if torsion is None and cert.short_curve_has_no_rational_2_torsion_modular_certificate(model,prime):torsion=prime
  try:sig=signature(cache,model,cloud,prime)
  except ValueError:continue
  before=len(pivots)
  for row in sig.rows:insert(pivots,row)
  if len(pivots)>before:signatures.append(asdict(sig))
 proof=checked_rank(model,cloud,[s['prime'] for s in signatures],torsion) if len(pivots)==len(cloud) and torsion else None
 result={'schema':'elliptic-curves.kihara-conic-control.v1','status':'PASS' if proof else 'UNRESOLVED_FINITE_QUOTIENT','conic_squarefree_part':str(d),'conic_scale':str(r),'balanced_integer_parameter':int(v),'conic_parameter_s':str(s0),'old_parent_T':str(K(data['base_T'])(s0)),
  'curve':list(map(str,model)),'points':[[str(c) for c in p.xy()] for p in pts],'coefficient_bits':max(abs(ZZ(c)).nbits() for c in E.a_invariants()),'point_transport_from_polynomial_model':str(transport),'prime_bound':997,'finite_mod2_rank':len(pivots),'rank_certificate':proof,'rank_lower_bound':13 if proof else None,
  'scope':'One deterministic balanced conic parameter: ceil(sqrt(squarefree part(c))) after removing square scaling. No score, parameter population, point enumeration, catalogue-informed selection or rank gain from search. Exact13-section evaluation and finite rank proof when successful; the conic direction is built into the new seed.'}
 assert not OUT.exists();OUT.write_text(json.dumps(result,indent=2)+'\n');print(result['status'],'v',v,'bits',result['coefficient_bits'],'rank',len(pivots),flush=True)
if __name__=='__main__':main()
