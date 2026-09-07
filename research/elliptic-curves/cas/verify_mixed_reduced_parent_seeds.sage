#!/usr/bin/env sage-python
"""Standalone specialization, retained-source, and rational-span equality replay."""
import argparse,json,hashlib
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(d):
 inp=json.loads((d/'intake-input.json').read_bytes());rows=json.loads((d/'intake-result.json').read_bytes())['rows'];R=PolynomialRing(QQ,'s');K=R.fraction_field()
 def f(v):return K(R(v['numerator'])/R(v['denominator']))
 a=list(map(f,inp['parent']['a_invariants']));P=[[f(v) for v in Q] for Q in inp['parent']['basis_weierstrass_coordinates']]
 for row in rows:
  path=d/(row['id']+'-seed.json');assert sha(path)==row['seed_sha256'];s=json.loads(path.read_bytes());E=EllipticCurve(QQ,s['curve']);basis=[E(Q) for Q in s['points']];assert len(basis)==row['initial_rank']
  if row['parent']=='det1092':
   t=QQ(row['parameter']);assert list(E.a_invariants())==[v(t) for v in a];original=[]
   for x,y in P:
    if x.denominator()(t)==0 or y.denominator()(t)==0:assert x.denominator()(t)==y.denominator()(t)==0;original.append(E(0))
    else:original.append(E([x(t),y(t)]))
   assert [None if Q.is_zero() else list(map(str,Q.xy())) for Q in original]==s['source_section_points']
   assert all(Q in original for Q in basis) and s['generic_points']==s['points']
   n=t.denominator();assert list(map(QQ,row['score_model']))==[0,0,0,E.a4()*n**8,E.a6()*n**12]
   aa,bb,cc,dd=map(QQ,inp['parameter_matrix']);assert QQ(row['original_parameter'])==(aa*t+bb)/(cc*t+dd)
   assert len(s['section_relations'])==17
   for rel in s['section_relations']:assert int(rel['denominator'])*original[rel['section_index']]==sum((int(v)*Q for v,Q in zip(rel['word'],basis)),E(0))
  else:
   source=inp['retained'][row['id']];assert source['curve']==s['curve'] and source['points']==s['points'] and source['generic_points']==s['generic_points']
  print('PASS seed',row['id'],len(basis),flush=True)
 out=d/'result.json';assert not out.exists();out.write_text(json.dumps(dict(status='PASS',rows=len(rows),input_sha256=sha(d/'intake-input.json'),intake_sha256=sha(d/'intake-result.json'),boundary='Standalone exact seed memberships, generic specializations and rational-span relations. Finite independence is certified separately.'),indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--directory',type=Path,required=True);main(p.parse_args().directory)
