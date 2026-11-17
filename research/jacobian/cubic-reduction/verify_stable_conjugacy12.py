#!/usr/bin/env python3
"""Independent low-jet certificate and finite coupled-translation replay."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from collections import Counter
import hashlib,json
from verify import add,mul,scale,variable,derivative,decode,compose,evaluate
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def determinant(A):
 A=[list(row) for row in A];out=Q(1)
 for j in range(len(A)):
  pivot=next((i for i in range(j,len(A)) if A[i][j]),None)
  if pivot is None:return Q(0)
  if pivot!=j:A[j],A[pivot]=A[pivot],A[j];out=-out
  d=A[j][j];out*=d
  for i in range(j+1,len(A)):
   c=A[i][j]/d
   if c:A[i]=[a-c*b for a,b in zip(A[i],A[j])]
 return out
def inverse(A):
 n=len(A);B=[list(row)+[Q(i==j) for j in range(n)] for i,row in enumerate(A)]
 for j in range(n):
  pivot=next(i for i in range(j,n) if B[i][j]);B[pivot],B[j]=B[j],B[pivot]
  B[j]=[x/B[j][j] for x in B[j]]
  for i in range(n):
   if i!=j and B[i][j]:
    c=B[i][j];B[i]=[x-c*y for x,y in zip(B[i],B[j])]
 inv=[row[n:] for row in B]
 assert all(sum(A[i][k]*inv[k][j] for k in range(n))==Q(i==j) for i in range(n) for j in range(n))
 return inv
def lowmul(p,q,cap=3):
 return {e:c for e,c in mul(p,q).items() if sum(e)<=cap}
def monomial(indices,n):
 ex=[0]*n
 for i in indices:ex[i]+=1
 return {tuple(ex):Q(1)}
def lowpull(indices,F):
 out={(0,)*len(next(iter(F[0]))):Q(1)}
 for i in indices:out=lowmul(out,F[i])
 return out
def main():
 folder=HERE/'stable-conjugacy12';cert=json.loads((folder/'certificate.json').read_text());path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());assert cert['input_sha256']==sha(path)
 F=[decode(row,12) for row in raw['coefficients']];X=[variable(i,12) for i in range(12)]
 H2=[{e:c for e,c in f.items() if sum(e)==2} for f in F];H3=[{e:c for e,c in f.items() if sum(e)==3} for f in F]
 assert all(f==add(add(X[i],H2[i]),H3[i]) for i,f in enumerate(F))
 rows=cert['quadratic_minor_rows'];cols=list(map(tuple,cert['quadratic_minor_columns']))
 assert rows==[1,2,3,4,6,7,8,9,10,11,12] and len(set(cols))==11
 mat=[[H2[i-1].get(e,Q(0)) for e in cols] for i in rows];det=determinant(mat)
 assert det==Q(cert['quadratic_minor_determinant'])==-288 and not H2[4]
 weights={tuple(row['exponent']):Q(row['weight']) for row in cert['functional']}
 def lam(p):return sum(w*p.get(e,Q(0)) for e,w in weights.items())
 count=0
 for ij in combinations_with_replacement(range(12),2):
  q=monomial(ij,12);d2={}
  for i in range(12):d2=add(d2,mul(derivative(q,i),H2[i]))
  assert lam(d2)==0,ij
  assert add(lowpull(ij,F),q,Q(-1))==d2
  count+=1
 assert count==78 and lam(H3[4])==1 and H3[4]==monomial([1,1,1],12)
 # Entire degree<=3 pullback basis, with no assumption on higher formal coefficients.
 checks=0
 for degree in range(4):
  for indices in combinations_with_replacement(range(12),degree):
   q=monomial(indices,12);diff=add(lowpull(indices,F),q,Q(-1))
   if degree in (0,3):assert not diff
   if degree==1:assert diff==add(H2[indices[0]],H3[indices[0]])
   checks+=1
 assert checks==455
 # Two stabilization variables: all cross-quadratic contributions disappear at Z=0.
 F14=[{e+(0,0):c for e,c in f.items()} for f in F]+[variable(12,14),variable(13,14)]
 cross=0
 for ij in combinations_with_replacement(range(14),2):
  if max(ij)<12:continue
  diff=add(lowpull(ij,F14),monomial(ij,14),Q(-1))
  assert not {e[:12]:c for e,c in diff.items() if not any(e[12:])}
  cross+=1
 assert cross==27
 for i in (12,13):assert lowpull([i],F14)==variable(i,14)
 # Independently reconstruct the finite affine experiment with generic composition,
 # not the primary binomial-shift function.
 data_path=HERE/'joint-translation12/checkpoint.json';data=json.loads(data_path.read_text())
 assert data['status']=='finite-centres-complete' and data['source_sha256']==sha(HERE/'joint_translation12.py') and data['input_sha256']==sha(path)
 J=[[derivative(f,j) for j in range(12)] for f in F];selected={tuple(row['centre']):row for row in data['selected']};hist=Counter();seen=set()
 for row in data['centres']:
  centre=tuple(map(Q,row['centre']));assert centre not in seen;seen.add(centre)
  assert sum(bool(v) for v in centre)>1 and all(evaluate(F[i],centre)==0 for i in range(3,12))
  A=[[evaluate(f,centre) for f in jr] for jr in J];assert determinant(A)==1;inv=inverse(A)
  shifted=[add(X[i],{(0,)*12:centre[i]}) for i in range(12)]
  pulled=[compose(f,shifted,12) for f in F]
  pulled=[add(p,{(0,)*12:evaluate(f,centre)},Q(-1)) for p,f in zip(pulled,F)]
  G=[]
  for rr in inv:
   out={}
   for c,p in zip(rr,pulled):out=add(out,p,c)
   G.append(out)
  for i,p in enumerate(G):
   assert not p.get((0,)*12,0)
   assert {e:c for e,c in p.items() if sum(e)==1}==X[i]
  score=sum(sum(e)>1 for p in G for e in p);assert score==row['nonlinear_total'];hist[score]+=1
  key=tuple(map(str,centre))
  if key in selected:
   rr=selected[key];assert G==[decode(p,12) for p in rr['coefficients']]
   assert inv==[[Q(c) for c in r] for r in rr['matrix']]
   pts=[[Q(v)-c for v,c in zip(p,centre)] for p in raw['points']]
   assert pts[0]!=pts[1] and [evaluate(p,pts[0]) for p in G]==[evaluate(p,pts[1]) for p in G]
 assert len(seen)==97 and min(hist)==145 and dict(sorted(hist.items()))=={int(k):v for k,v in data['score_histogram'].items()}
 old=json.loads((HERE/'integer-validation12.json').read_text());assert old['witness_sha256']==sha(path) and old['determinant_one'] and old['nonlinear_total']==36
 receipt={'status':'pass','claim':'KCUB25','quadratic_rank':11,'quadratic_minor_determinant':str(det),'sole_surviving_linear_jet_before_cubic':5,'functional_annihilates_quadratic_lie_image':78,'functional_on_cubic_row5':'1','pullback_monomials_through_degree3':checks,'two_identity_variable_cross_controls':cross,
 'all_formal_invariant_functions_have_zero_linear_X_jet':True,'arbitrary_stabilization_proof_written':True,'minimum_active_dimension_in_stable_conjugacy_class':12,'global_cubic_dimension_minimum_proved':False,'formal_verification':False,
 'finite_translation_replay':{'centres':97,'best':145,'stored_maps_reconstructed':len(selected),'stored_collision_checks':len(selected),'scope':'Finite coupled carrier-zero affine centres only; unequal left/right normalizations are outside the universal conjugacy theorem. No global translation/affine optimum.'},
 'checker_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'certificate_sha256':sha(folder/'certificate.json'),'proof_sha256':sha(HERE/'STABLE_CONJUGACY_OBSTRUCTION.md'),'input_sha256':sha(path),'base_validation_sha256':sha(HERE/'integer-validation12.json'),'translation_source_sha256':sha(HERE/'joint_translation12.py'),'translation_checkpoint_sha256':sha(data_path),'positive_bounds_unchanged':True,'goal_complete':False}
 (folder/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print('PASS:rank11/minor-288;78 Lie controls;455 pullback controls;stable conjugacy needs12 active dimensions;97 affine centres replayed(best145).')
if __name__=='__main__':main()
