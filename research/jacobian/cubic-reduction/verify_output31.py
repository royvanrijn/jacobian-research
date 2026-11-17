#!/usr/bin/env python3
"""Independent exact check of output-frame cubic and sharp count tradeoffs."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,permutations
from collections import Counter
import json,hashlib
from verify import decode,add,scale,derivative,evaluate,identity,madd,mmul,determinant3,compose,variable
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def eye(n):return [[Q(i==j) for j in range(n)] for i in range(n)]
def numeric_det(A):
 A=[row[:] for row in A];d=Q(1)
 for j in range(len(A)):
  k=next(i for i in range(j,len(A)) if A[i][j])
  if k!=j:A[j],A[k]=A[k],A[j];d=-d
  p=A[j][j];d*=p
  for i in range(j+1,len(A)):
   c=A[i][j]/p
   for k in range(j,len(A)):A[i][k]-=c*A[j][k]
 return d

def main():
 path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(p,12) for p in raw['coefficients']];wp=HERE/'output-frame31/witness.json';w=json.loads(wp.read_text());assert w['source_sha256']==sha(HERE/'construct_output31.py') and w['input_sha256']==sha(path)
 oldpath=HERE/'integer-validation12.json';old=json.loads(oldpath.read_text());assert old['status']=='pass' and old['source_sha256']==sha(HERE/'verify_integer12.py') and old['arithmetic_sha256']==sha(HERE/'verify.py') and old['witness_sha256']==sha(path) and old['full_polynomial_determinant'] and old['determinant_one']
 operations=[(1,6,Q(-3)),(1,9,Q(-12)),(2,7,Q(-8)),(8,7,Q(3)),(2,8,Q(-4))];assert w['operations']==[{'target':i+1,'from':j+1,'coefficient':str(c)} for i,j,c in operations]
 M=eye(12);points=[[Q(v) for v in p] for p in raw['points']];assert points[0]!=points[1];trades=[]
 for k,stage in enumerate(w['stages']):
  if k:
   i,j,c=operations[k-1];M[i]=[a+c*b for a,b in zip(M[i],M[j])]
  assert [[Q(v) for v in row] for row in stage['linear_matrix']]==M and numeric_det(M)==1
  G=[]
  for row in M:
   p={}
   for c,f in zip(row,F):p=add(p,f,c)
   G.append(p)
  assert G==[decode(p,12) for p in stage['coefficients']]
  assert all(c.denominator==1 and sum(e)>=1 for p in G for e,c in p.items())
  J0=[[derivative(f,j).get((0,)*12,Q(0)) for j in range(12)] for f in G];assert J0==M
  linear=[sum(sum(e)==1 for e in p) for p in G];nonlinear=[sum(sum(e)>=2 for e in p) for p in G]
  assert linear==stage['linear_counts'] and nonlinear==stage['nonlinear_counts']
  assert sum(linear)==stage['linear_total']==12+k and sum(nonlinear)==stage['nonlinear_total']==36-k and stage['total']==48
  images=[[evaluate(p,point) for p in G] for point in points];assert images[0]==images[1]==list(map(Q,w['common_image']))
  assert max(sum(e) for p in G for e in p)==3
  trades.append([sum(linear),sum(nonlinear)])
 # Independent simultaneous description, using only the original raw map.
 expected=list(F);expected[1]=add(add(F[1],F[6],Q(-3)),F[9],Q(-12));expected[2]=add(add(F[2],F[8],Q(-4)),F[7],Q(-20));expected[8]=add(F[8],F[7],Q(3));assert G==expected
 explicit=eye(12)
 for i,j,c in [(1,6,-3),(1,9,-12),(2,8,-4),(2,7,-20),(8,7,3)]:explicit[i][j]=Q(c)
 assert M==explicit and M!=eye(12)
 A=[[M[i][j]-Q(i==j) for j in range(12)] for i in range(12)]
 def numeric_mul(A,B):return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
 A2=numeric_mul(A,A);assert not any(v for row in numeric_mul(A2,A) for v in row)
 Minv=[[Q(i==j)-A[i][j]+A2[i][j] for j in range(12)] for i in range(12)];assert numeric_mul(M,Minv)==eye(12) and numeric_mul(Minv,M)==eye(12)
 con=w['identity_normalized_conjugate'];assert [[Q(v) for v in row] for row in con['inverse_linear_matrix']]==Minv
 X=[variable(j,12) for j in range(12)];sub=[]
 for row in Minv:
  p={}
  for a,x in zip(row,X):p=add(p,x,a)
  sub.append(p)
 C=[compose(p,sub,12) for p in G];assert C==[decode(p,12) for p in con['coefficients']]
 assert [[derivative(p,j).get((0,)*12,Q(0)) for j in range(12)] for p in C]==eye(12)
 conpoints=[[sum(a*b for a,b in zip(row,point)) for row in M] for point in points];assert conpoints==[list(map(Q,p)) for p in con['points']] and conpoints[0]!=conpoints[1]
 assert all([evaluate(p,point) for p in C]==list(map(Q,con['common_image'])) for point in conpoints)
 concounts=[sum(sum(e)>=2 for e in p) for p in C];assert concounts==con['nonlinear_counts']==[7,12,22,6,10,1,3,6,1,1,2,3] and sum(concounts)==con['nonlinear_total']==74
 # Full all-variable polynomial determinant of the final map; no inherited determinant assertion is trusted here.
 J=[[derivative(f,j) for j in range(12)] for f in G];D=[row[3:] for row in J[3:]];I=identity(9,12);N=madd(D,I,Q(-1));assert all(not N[i][j] for i in range(9) for j in range(i,9))
 inv=I;power=I
 for k in range(1,9):power=mmul(power,N);inv=madd(inv,power,Q((-1)**k))
 assert mmul(power,N)==[[{} for j in range(9)] for i in range(9)] and mmul(D,inv)==I and mmul(inv,D)==I
 schur=madd([row[:3] for row in J[:3]],mmul(mmul([row[3:] for row in J[:3]],inv),[row[:3] for row in J[3:]]),Q(-1));assert determinant3(schur,12)=={(0,)*12:Q(1)}
 # Reconstruct every source-column vector used by the written universal lower bounds.
 exponents=set().union(*(set(p) for p in F));columns={e:tuple(p.get(e,Q(0)) for p in F) for e in sorted(exponents)}
 nonlinear={e:v for e,v in columns.items() if sum(e)>=2};assert len(columns)==42 and len(nonlinear)==30
 def exponent(**kwargs):
  e=[0]*12
  for j,k in kwargs.items():e[int(j[1:])-1]=k
  return tuple(e)
 E=[tuple(Q(i==j) for i in range(12)) for j in range(12)]
 def vc(*terms):return tuple(sum(Q(c)*E[j][i] for j,c in terms) for i in range(12))
 exceptional=[exponent(x1=1,x2=1),exponent(x1=2,x3=1),exponent(x1=1,x10=1)]
 assert [nonlinear[e] for e in exceptional]==[vc((3,3),(11,1)),E[3],E[11]]
 for a,b in combinations([nonlinear[e] for e in exceptional],2):assert a[3]*b[11]-a[11]*b[3]
 blockrows=[{0},{1,6,9},{2,7,8},{3,11},{4},{5},{10}];blocks=[]
 for rows in blockrows:
  block={e:v for e,v in columns.items() if any(v[j] for j in rows)};assert all({j for j,a in enumerate(v) if a}<=rows for v in block.values());blocks.append(block)
 assert sum(len(b) for b in blocks)==42 and [len(b) for b in blocks]==[5,13,13,5,2,2,2]
 # Rank-one blocks: every nonzero column survives any invertible output frame.
 for b,axis in [(blocks[0],0),(blocks[4],4),(blocks[5],5),(blocks[6],10)]:assert all(v[axis] and all(not a for j,a in enumerate(v) if j!=axis) for v in b.values())
 heavy=lambda b,axis:[v for v in b.values() if v[axis] and all(not a for j,a in enumerate(v) if j!=axis)]
 assert len(heavy(blocks[1],1))==9 and len(heavy(blocks[2],2))==9
 lights1={exponent(x7=1):E[6],exponent(x10=1):E[9],exponent(x1=1,x3=1):vc((1,3),(6,1)),exponent(x1=1,x8=1):vc((1,12),(9,1))}
 lights2={exponent(x8=1):E[7],exponent(x9=1):E[8],exponent(x2=2):vc((2,8),(7,1),(8,-3)),exponent(x1=1,x5=1):vc((2,8),(8,2))}
 for b,axis,lights in [(blocks[1],1,lights1),(blocks[2],2,lights2)]:
  assert all(b[e]==v for e,v in lights.items()) and len(b)==len(heavy(b,axis))+len(lights)
 # Any three second-block light vectors span its3-space, so a functional nonzero on the heavy axis cannot kill three.
 triples=[];lights=list(lights2.values())
 for indices in combinations(range(4),3):
  matrix=[[lights[j][i] for j in indices] for i in [2,7,8]];d=numeric_det(matrix);assert d;triples.append(str(d))
 assert [tuple(v[j] for j in [7,8]) for v in lights]==[(Q(1),Q(0)),(Q(0),Q(1)),(Q(1),Q(-3)),(Q(0),Q(2))]
 assert Counter(blocks[3].values())==Counter({E[3]:2,E[11]:2,vc((3,3),(11,1)):1})
 assert w['output_block_rows']==[[i+1 for i in sorted(rows)] for rows in blockrows] and w['total_block_lower_bounds']==[5,15,16,6,2,2,2] and sum(w['total_block_lower_bounds'])==48
 receipt={'status':'pass','claim':'KCUB30','checker_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'input_sha256':sha(path),'witness_sha256':sha(wp),'producer_sha256':sha(HERE/'construct_output31.py'),'base_validation_sha256':sha(oldpath),'dimension':12,'degree':3,'integral_coefficients_and_collision':True,'full_symbolic_determinant_one':True,'linear_part_unipotent_nonidentity':True,'linear_occurrences':17,'nonlinear_occurrences':31,'total_occurrences':48,'nonlinear_source_columns':30,'full_source_columns':42,'exceptional_projective_columns':3,'block_column_sizes':[len(b) for b in blocks],'block_total_lower_bounds':[5,15,16,6,2,2,2],'second_block_light_triple_minors':triples,'pareto_pairs_linear_nonlinear':trades,'universal_output_only_lower_bounds':{'linear':12,'nonlinear':31,'total':48},'identity_linear_part_bound_unchanged':36,'checked_identity_normalized_conjugate_nonlinear_total':74,'checked_identity_normalized_conjugate_counts':concounts,'source_frame_fixed':True,'goal_complete':False,'scope':'Finite exact premises for a written all-output-GL counting theorem. No arbitrary source-frame, identity-normalized improvement, smaller dimension, total-term improvement or global optimality claim.'}
 proof=HERE/'OUTPUT_FRAME_TRADEOFF.md'
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (HERE/'output-frame31/validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS integral cubic31 nonlinear+17 linear=48; full determinant and collision; six Pareto witnesses; exact output-only lower-bound premises',flush=True)
if __name__=='__main__':main()
