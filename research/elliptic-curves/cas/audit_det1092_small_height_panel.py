#!/usr/bin/env python3
"""Fixed projective-height2 diagnostic, before any new candidate selection."""
from pathlib import Path
from fractions import Fraction as F
import math
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/det1092-height2-panel-v1';SOURCE=ROOT/'artifacts/local/elliptic-curves/det1092-point-pilot-v2/parent-sections.json'
def main():
 rows=sorted({F(m,n) for m in range(-2,3) for n in range(1,3) if math.gcd(m,n)==1});parameters=list(map(str,rows))+['infinity'];assert len(parameters)==8
 p=dict(source_sha256=cert.hashed(SOURCE),script_sha256=cert.hashed(Path(__file__)),parameters=parameters,scope='Complete projective-height2 j-height diagnostic, including the excluded anchor0 and previously tested1 as labelled controls. No point search, score, target point or promotion. All parameters fixed before evaluation.')
 assert not (D/'protocol.json').exists();checkpoint(D/'protocol.json',p)
 parent=cert.read(SOURCE);data=[]
 def ev(c,t):
  x=F(0)
  for a in reversed(c):x=x*t+F(a)
  return x
 for t in parameters:
  a=[]
  for r,weight in zip(parent['a_invariants'],[2,4,6,8,12]):
   if t=='infinity':
    assert r['denominator']==['1'];a.append(F(r['numerator'][weight]) if len(r['numerator'])>weight else F(0))
   else:a.append(ev(r['numerator'],F(t))/ev(r['denominator'],F(t)))
  v=cert.weierstrass_invariants(a);j=v['c4']**3/v['discriminant'];bn=abs(j.numerator).bit_length();bd=j.denominator.bit_length();floor=max(1,(bn-25+5)//6,(bd-16+5)//6)
  data.append(dict(parameter=t,j_numerator_bits=bn,j_denominator_bits=bd,quartic_coefficient_bits_floor=floor,excluded_anchor=t=='0',already_point_searched=t=='1',curve=list(map(str,a))))
 checkpoint(D/'result.json',dict(status='PASS',protocol_hash=digest(p),rows=data));print([(r['parameter'],r['j_numerator_bits'],r['quartic_coefficient_bits_floor']) for r in data],flush=True)
if __name__=='__main__':main()
