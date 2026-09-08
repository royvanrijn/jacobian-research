#!/usr/bin/env sage-python
"""Exact raw slope identity and new degree-two centre classes, with finite-group replay."""
import argparse,json,hashlib
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import QQ,PolynomialRing,EllipticCurve
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(args):
 p=json.loads(args.protocol.read_bytes());assert p['verifier_sha256']==sha(Path(__file__)) and p['input_sha256']==sha(args.input) and p['rank_verifier_sha256']==sha(args.rank_verifier)
 data=json.loads(args.input.read_bytes());rankcheck=SourceFileLoader('independent_finite_rank',str(args.rank_verifier)).load_module()
 R=PolynomialRing(QQ,names=['A','c','d','m','X']);A,c,d,m,X=R.gens();B=d*d-c**3-A*c
 q=X*X+(c-m*m)*X+(c*c+A+m*m*c+2*m*d)
 assert X**3+A*X+B-(m*(X-c)-d)**2==(X-c)*q
 discriminant=(c-m*m)**2-4*(c*c+A+m*m*c+2*m*d)
 assert discriminant==m**4-6*c*m*m-8*d*m-3*c*c-4*A
 rows=[]
 for row in data['rows']:
  seed=row['seed'];E=EllipticCurve(QQ,seed['curve']);points=[E(P) for P in seed['points']];r=len(points);old=row['old_dimension'];assert 0<old<r and seed['points'][:old]==seed['generic_points']
  proof=seed['rank_certificate'];cloud=dict(curve=seed['curve'],points=seed['points'],signatures=proof['signatures'],rank_certificate=proof,rank_lower_bound=r,independent_column_indices=list(range(r)))
  cf=args.output.parent/(row['id']+'-rank-input.json');assert not cf.exists();cf.write_text(json.dumps(cloud,indent=2)+'\n');finite=rankcheck.check_one(cf);assert finite['union_finite_rank']==r
  masks=[]
  for chart in row['maps']:
   centre=chart['centre'];word=list(map(int,centre['representative']));assert len(word)==r;mask=sum((v%2)<<i for i,v in enumerate(word));assert mask==centre['parity'] and mask>>old and mask not in masks;masks.append(mask)
   C=sum((v*P for v,P in zip(word,points)),E(0));assert not C.is_zero();x,y=C.xy()
   assert list(map(QQ,chart['raw_coefficients']))==[-3*x*x-4*E.a4(),-8*y,-6*x,0,1]
  assert len(masks)==49
  rows.append(dict(id=row['id'],certified_dimension=r,old_dimension=old,selected_distinct_new_centre_classes=49,available_classes_outside_old_span=2**r-2**old,parity_masks=masks,finite_groups=finite['groups'],no_two_torsion_prime=finite['no_two_torsion_prime']))
 out=dict(status='PASS',input_sha256=sha(args.input),protocol_sha256=sha(args.protocol),symbolic_raw_quartic_identity=True,rows=rows,total_distinct_new_centre_classes=sum(r['selected_distinct_new_centre_classes'] for r in rows),boundary='Classical degree-two quotient classification applied to exact finite-independent point subgroups. Every selected centre has a class in E(Q)/2E(Q) outside the old seed image; the associated degree-two map cannot be equivalent to an old-subgroup-centre map by a rational curve automorphism and PGL2 change. These are presentations of the same elliptic curve, not new curves, Selmer directions, or a proof of finite-box point non-overlap or improved detection.')
 assert not args.output.exists();args.output.write_text(json.dumps(out,indent=2)+'\n');print('PASS symbolic slope identity and',out['total_distinct_new_centre_classes'],'new centre classes',flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser()
 for n in ['protocol','input','rank-verifier','output']:ap.add_argument('--'+n,type=Path,required=True)
 main(ap.parse_args())
