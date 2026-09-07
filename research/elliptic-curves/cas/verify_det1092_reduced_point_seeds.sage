#!/usr/bin/env sage-python
"""Independent exact specialization and short-coordinate replay of two fixed fibres."""
import argparse,json
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve

def main(parent,paths,output):
 d=json.loads(parent.read_bytes());R=PolynomialRing(QQ,'t');K=R.fraction_field()
 def value(r):return K(R(r['numerator'])/R(r['denominator']))
 a=[value(r) for r in d['a_invariants']];points=[[value(r) for r in P] for P in d['basis_weierstrass_coordinates']];rows=[]
 assert len(points)==17 and set(d)=={'a_invariants','basis_weierstrass_coordinates','source_sha256'}
 for path in paths:
  seed=json.loads(path.read_bytes());t=QQ(seed['parameter']);assert t in [QQ(1),QQ(1009)/101]
  invariants=[v(t) for v in a];E=EllipticCurve(QQ,invariants);P=[E([v(t) for v in Q]) for Q in points]
  assert list(map(str,invariants))==seed['long_curve'] and [list(map(str,Q.xy())) for Q in P]==seed['long_points']
  S=EllipticCurve(QQ,seed['curve']);assert list(S.a_invariants())==[0,0,0,-E.c4()/48,-E.c6()/864]
  actual=[]
  for Q in P:
   x,y=Q.xy();actual.append(S([x+E.b2()/12,y+(E.a1()*x+E.a3())/2]))
  assert [list(map(str,Q.xy())) for Q in actual]==seed['original_generic_points']
  subset=[actual[i] for i in seed['independent_indices']]
  assert [list(map(str,Q.xy())) for Q in subset]==seed['points']==seed['generic_points']
  if t==1:
   assert seed['independent_indices']==list(range(14))+[16]
   assert actual[14]==actual[2]+actual[5]-actual[6]+actual[9]-actual[13]
   assert actual[15]==-actual[3]-actual[5]+actual[7]+actual[8]-actual[12]+2*actual[13]
  else: assert seed['independent_indices']==list(range(17))
  rows.append(dict(parameter=str(t),point_count=len(subset),original_section_count=17,model_bits=seed['model_coefficient_bits']))
 assert len(rows)==2 and rows[0]['parameter']!=rows[1]['parameter']
 assert not output.exists();output.write_text(json.dumps(dict(status='PASS',rows=rows,boundary='Independent generic-section specialization and exact coordinate transport. Point independence is verified separately on the complete clouds.'),indent=2)+'\n')
 print('PASS two determinant1092 exact seed transports and unit relations')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--parent',type=Path,required=True);p.add_argument('--seeds',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.parent,a.seeds,a.output)
