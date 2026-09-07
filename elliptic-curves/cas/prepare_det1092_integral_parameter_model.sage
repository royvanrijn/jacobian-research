#!/usr/bin/env sage-python
"""Bounded source-only integral parameter model and height2 local-scaling support."""
import json,hashlib,sys
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,lcm,gcd,ceil
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from research_runtime.store import checkpoint,digest
D=ROOT/'artifacts/local/elliptic-curves/det1092-integral-parameter-model-v1';SOURCE=ROOT/'artifacts/local/elliptic-curves/det1092-point-pilot-v2/parent-sections.json'
def hashed(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 assert not (D/'result.json').exists();p=json.loads((D/'protocol.json').read_bytes());assert p['source_sha256']==hashed(SOURCE) and p['script_sha256']==hashed(Path(__file__))
 parent=json.loads(SOURCE.read_bytes());R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()
 def f(r):return K(R(r['numerator'])/R(r['denominator']))
 E=EllipticCurve(K,[f(r) for r in parent['a_invariants']]);A=R(-E.c4()/48);B=R(-E.c6()/864);assert A.degree()==8 and B.degree()==12
 DA=lcm([a.denominator() for a in A]);DB=lcm([b.denominator() for b in B]);den=lcm(DA,DB)
 data=dict(status='FACTORING_COEFFICIENT_DENOMINATORS',protocol_hash=digest(p),coefficient_denominator=str(den),rows=[]);checkpoint(D/'result.json',data)
 factors=list(ZZ(den).factor(proof=True));assert all(prime.is_prime(proof=True) for prime,e in factors)
 L=ZZ(1)
 for prime,e in factors:L*=prime**max(ceil(QQ(DA.valuation(prime))/4),ceil(QQ(DB.valuation(prime))/6))
 F=R(A*L**4);G=R(B*L**6);assert all(c.denominator()==1 for c in list(F)+list(G))
 data.update(status='FACTORING_FIXED_HEIGHT2_CONTENTS',denominator_factors=[[str(q),int(e)] for q,e in factors],scale_L=str(L),A_coefficients=list(map(str,F)),B_coefficients=list(map(str,G)));checkpoint(D/'result.json',data)
 allprimes={int(q) for q,e in factors}
 for m,n in p['projective_parameters']:
  aa=sum(ZZ(F[i])*ZZ(m)**i*ZZ(n)**(8-i) for i in range(9));bb=sum(ZZ(G[i])*ZZ(m)**i*ZZ(n)**(12-i) for i in range(13));c=gcd(aa,bb)
  row=dict(parameter=[m,n],coefficient_gcd=str(c),status='FACTORING');data['rows'].append(row);checkpoint(D/'result.json',data)
  fac=list(c.factor(proof=True));assert all(q.is_prime(proof=True) for q,e in fac)
  row.update(status='PASS',factors=[[str(q),int(e)] for q,e in fac],removable_short_scalings=[[str(q),int(min(aa.valuation(q)//4,bb.valuation(q)//6))] for q,e in fac if min(aa.valuation(q)//4,bb.valuation(q)//6)>0]);allprimes.update(int(q) for q,e in fac);checkpoint(D/'result.json',data)
  print('HEIGHT2 SUPPORT',m,n,'gcd bits',c.nbits(),'prime bits',[q.nbits() for q,e in fac],flush=True)
 data.update(status='PASS',candidate_primes=sorted(allprimes),boundary='Only coefficient denominators and eight fixed projective-height2 specialization gcds factored. No whole discriminant/resultant factorization, score, point search, parameter promotion or claim of complete global scaling support.');checkpoint(D/'result.json',data);print('PASS integral model and fixed support',len(allprimes),flush=True)
if __name__=='__main__':main()
