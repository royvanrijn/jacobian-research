#!/usr/bin/env sage-python
"""Four fixed fresh Kihara fibres: exact seed models and finite rank intake."""
import json,hashlib,sys,argparse
from pathlib import Path
from fractions import Fraction
from dataclasses import asdict
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,gcd,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas';sys.path[:0]=[str(CAS),str(ROOT/'elliptic-curves')]
from ecsearch.kihara import kihara_rank14_replay
import certify_compact_r17_candidates as cert
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from audit_recorded_point_mod2_rank_v3 import signature,insert,_primes_up_to
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint
D=ROOT/'artifacts/local/elliptic-curves/kihara-fresh-fibres-v1';SOURCE=ART/'kihara_global_parent_geometry_v1.json'
def intake(index):
 protocol=json.loads((D/'protocol.json').read_text());t=QQ(protocol['parameters'][index]);out=D/('fibre'+str(index))/'seed.json';assert not out.exists()
 R=PolynomialRing(QQ,'t');g=json.loads(SOURCE.read_text())['geometry'];A=R(g['minimal_A']);B=R(g['minimal_B']);scale=R(g['removed_scale'])
 # Remove only the common constant fourth/sixth powers, exactly.
 ca=gcd(A.list());cb=gcd(B.list());constant=QQ(1)
 primes=set(ZZ(ca.numerator()).prime_divisors()+ZZ(cb.numerator()).prime_divisors()+ZZ(ca.denominator()).prime_divisors()+ZZ(cb.denominator()).prime_divisors())
 for p in primes:constant*=QQ(p)**min(ca.valuation(p)//4,cb.valuation(p)//6)
 A=A/constant**4;B=B/constant**6
 original=kihara_rank14_replay(Fraction(str(t)));c=QQ(original.weierstrass_coefficients[1]);u=3*t**5/(scale(t)*constant)
 E=EllipticCurve(QQ,[A(t),B(t)]);points=[E([u*u*(QQ(x)+c/3),u**3*QQ(y)]) for x,y in original.weierstrass_points]
 assert len(points)==14
 # Integral raw short model, without global factorization/minimalization.
 den=lcm([v.denominator() for v in E.a_invariants()]);Eint=EllipticCurve(QQ,[den**4*E.a4(),den**6*E.a6()]);points=[Eint([den**2*P[0],den**3*P[1]]) for P in points]
 model=tuple(cert.F(str(c)) for c in Eint.a_invariants());cloud=[tuple(cert.F(str(c)) for c in P.xy()) for P in points]
 cache=ReductionCache(MemoryFactStore());pivots={};signatures=[];torsion=None
 for prime in _primes_up_to(997):
  if prime==2:continue
  if torsion is None and cert.short_curve_has_no_rational_2_torsion_modular_certificate(model,prime):torsion=prime
  try:sig=signature(cache,model,cloud,prime)
  except ValueError:continue
  before=len(pivots)
  for bits in sig.rows:insert(pivots,bits)
  if len(pivots)>before:signatures.append(asdict(sig))
 rank=len(pivots);proof=checked_rank(model,cloud,[s['prime'] for s in signatures],torsion) if rank==14 and torsion is not None else None
 result={'schema':'elliptic-curves.kihara-fresh-seed.v1','status':'PASS' if proof else 'UNRESOLVED_SEED_QUOTIENT','parameter':str(t),'family':'kihara-rank14','curve':list(map(str,model)),
 'points':[[str(c) for c in P] for P in cloud],'generic_points':[[str(c) for c in P] for P in cloud],'finite_mod2_rank':rank,'rank_certificate':proof,'rank_lower_bound':14 if proof else None,
 'model_coefficient_bits':max(abs(c.numerator).bit_length() for c in model),'constant_scale':str(constant),'generic_scale_value':str(scale(t)),'point_scale_from_original':str(u*den),
 'scope':'Fixed new path parameter; fourteen source points exactly mapped. Full mod2 independence required before point-search admission. A finite quotient deficit is not a rational dependence or rank upper bound.'}
 checkpoint(out,result);print('fibre',index,'t',t,'bits',result['model_coefficient_bits'],'mod2',rank,result['status'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);a=p.parse_args();intake(a.index)
