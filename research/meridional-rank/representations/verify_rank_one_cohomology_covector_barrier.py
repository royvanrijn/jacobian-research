#!/usr/bin/env python3
"""Exact controls for the all-field universal rank-one cohomology barrier."""
if not __debug__:raise RuntimeError('Assertions must be enabled')
from pathlib import Path
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=HERE/'rank_one_cohomology_covector_barrier_verification.json')
args=parser.parse_args()
ledger_path=HERE/'rank_one_cohomology_covector_barrier_source_ledger.json'
ledger=json.loads(ledger_path.read_text())
for source in ledger['sources']:
 p=HERE/source['file']
 assert p.stat().st_size==source['bytes']
 assert hashlib.sha256(p.read_bytes()).hexdigest()==source['sha256']

def eye(n):return[[int(i==j)for j in range(n)]for i in range(n)]
def add(A,B,p,sign=1):return[[(a+sign*b)%p for a,b in zip(x,y)]for x,y in zip(A,B)]
def mm(A,B,p):return[[sum(a*b for a,b in zip(row,col))%p for col in zip(*B)]for row in A]
def rref(A,p):
 A=[row.copy()for row in A];piv=[];r=0
 for j in range(len(A[0])):
  k=next((i for i in range(r,len(A))if A[i][j]%p),None)
  if k is None:continue
  A[r],A[k]=A[k],A[r];iv=pow(A[r][j]%p,-1,p);A[r]=[(iv*a)%p for a in A[r]]
  for i in range(len(A)):
   if i==r or not A[i][j]%p:continue
   q=A[i][j]%p;A[i]=[(a-q*b)%p for a,b in zip(A[i],A[r])]
  piv.append(j);r+=1
  if r==len(A):break
 return A[:r],piv

def rank(A,p):return len(rref(A,p)[1])
def nullspace(A,p):
 R,piv=rref(A,p);free=[j for j in range(len(A[0]))if j not in piv];basis=[]
 for j in free:
  v=[0]*len(A[0]);v[j]=1
  for row,col in zip(R,piv):v[col]=-row[j]%p
  basis.append(v)
 return basis

def roots(M,p):
 D=add(eye(len(M)),M,p,-1);j=next(j for j in range(len(M))if any(row[j]for row in D));u=[row[j]for row in D];i=next(i for i,a in enumerate(u)if a);f=[a*pow(u[i],-1,p)%p for a in D[i]]
 assert D==[[(a*b)%p for b in f]for a in u]
 return u,f,i

def jacobian(xs,word,p):
 n=len(xs);d=len(xs[0]);I=eye(d);big=eye(n*d);t=[[row.copy()for row in M]for M in xs];co=[big[j*d:(j+1)*d]for j in range(n)]
 for l in word:
  j=abs(l)-1;A,B=t[j:j+2];f,g=co[j:j+2]
  if l>0:
   C=mm(mm(B,A,p),B,p);ng=add(mm(B,f,p),mm(mm(B,add(A,I,p,-1),p),g,p),p);t[j:j+2]=[B,C];co[j:j+2]=[g,ng]
  else:
   C=mm(mm(A,B,p),A,p);nf=add(mm(add(I,C,p,-1),f,p),mm(A,g,p),p);t[j:j+2]=[C,A];co[j:j+2]=[nf,f]
 assert t==xs
 D=[row for block in co for row in block]
 return add(D,big,p,-1)

def braid_matrices(xs,word,p):
 t=xs.copy()
 for l in word:
  j=abs(l)-1;A,B=t[j:j+2]
  t[j:j+2]=[B,mm(mm(B,A,p),B,p)]if l>0 else[mm(mm(A,B,p),A,p),A]
 return t

def red(w):
 a=[]
 for i in w:
  if a and a[-1]==-i:a.pop()
  else:a.append(i)
 return tuple(a)
def wi(w):return tuple(-a for a in w[::-1])
def wm(*ws):return red(a for w in ws for a in w)
def free_action(values,word):
 t=list(values)
 for l in word:
  j=abs(l)-1;a,b=t[j:j+2];t[j:j+2]=[b,wm(wi(b),a,b)]if l>0 else[wm(a,b,wi(a)),a]
 return t

inputs=HERE/'rank_one_cohomology_barrier_inputs.json';data=json.loads(inputs.read_text());p=data['prime'];xs=data['input_matrices'];word=data['base_braid'];d=6;b=8;block_count=2;c=b-block_count
assert b>=4*block_count
U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3]
assert word[:28]==U+[i+4 if i>0 else i-4 for i in U]
free=[(i,)for i in range(1,9)];tail=free_action(free,word[28:]);assert all(tail[i]==free[i]for i in[0,3,4,7])
perm=list(range(8))
for l in word:j=abs(l)-1;perm[j],perm[j+1]=perm[j+1],perm[j]
seen={0};j=perm[0]
while j not in seen:seen.add(j);j=perm[j]
assert len(seen)==8
for M in xs:
 assert mm(M,M,p)==eye(d)and rank(add(eye(d),M,p,-1),p)==1
A=jacobian(xs,word,p);Z=nullspace(A,p)
u,f,pivot=roots(xs[0],p);constraint=[]
for i in range(d):
 if i==pivot:continue
 row=[0]*(d*b);row[i]=u[pivot];row[pivot]=-u[i]%p;constraint.append(row)
ZR=nullspace(A+constraint,p);r=len(ZR)
assert len(Z)==11 and r==6 and len(Z)-d+1==r==c
original_covectors=[roots(M,p)[1]for M in xs]
coboundary_matrix=[row for M in xs for row in add(M,eye(d),p,-1)]
original_covector_rank=rank(original_covectors,p)
coboundary_rank=rank(coboundary_matrix,p)
global_invariant_dimension=len(nullspace(coboundary_matrix,p))
assert original_covector_rank==coboundary_rank==d-global_invariant_dimension==5
assert global_invariant_dimension==1
assert mm(A,coboundary_matrix,p)==[[0]*d for _ in A]
assert coboundary_rank<=r
aug=[];cov=[];coeff=[]
for k,M in enumerate(xs):
 u,f,pivot=roots(M,p);ai=[z[k*d+pivot]*pow(u[pivot],-1,p)%p for z in ZR]
 E=[[z[k*d+i]for z in ZR]for i in range(d)]
 assert E==[[(a*v)%p for v in ai]for a in u]
 hat=[M[i]+E[i]for i in range(d)]+[[0]*d+row for row in eye(r)]
 assert rank(add(eye(d+r),hat,p,-1),p)==1 and mm(hat,hat,p)==eye(d+r)
 hatcov=f+[-a%p for a in ai];hatu=u+[0]*r
 assert add(eye(d+r),hat,p,-1)==[[a*v%p for v in hatcov]for a in hatu]
 aug.append(hat);cov.append(hatcov);coeff.append(ai)
assert braid_matrices(aug,word,p)==aug
assert rank(coeff,p)==r and rank(cov,p)==6
# mu=1 / characteristic-two scope has no diagonalizability requirement.
M=[[1,1],[0,1]];H=[[1,1,1],[0,1,0],[0,0,1]]
assert rank(add(eye(2),M,2,-1),2)==1 and rank(add(eye(3),H,2,-1),2)==1
assert mm(H,H,2)==eye(3)
# An arbitrary cocycle translation outside the root line breaks rank one.
bad=[[p-1,0,0],[0,1,1],[0,0,1]]
assert rank(add(eye(3),bad,p,-1),p)==2
proof=HERE/'RANK_ONE_COHOMOLOGY_COVECTOR_BARRIER.md'
receipt={'status':'passed','arithmetic':'exact finite-field matrices, complete nullspaces and freely reduced words',
 'case':{'field_characteristic':p,'V_dimension':d,'Z1_dimension':len(Z),'Zroot_dimension':r,'extended_dimension':d+r,'extended_covector_span':rank(cov,p),'projected_cocycle_coefficient_rank':rank(coeff,p),'original_covector_span':original_covector_rank,'coboundary_map_rank':coboundary_rank,'global_invariant_dimension':global_invariant_dimension,'b':b,'k':block_count},
 'controls':{'whole_augmented_braid_closure':True,'exact_protected_free_words':True,'one_knot_permutation_cycle':True,'mu_one_characteristic_two':True,'nonroot_translation_rejected':True,'all_imported_source_pins':True,'reverse_bridge_covectors_coboundaries_invariants':True},
 'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),inputs,proof,ledger_path]+[HERE/s['file']for s in ledger['sources']]},
 'scope':'finite controls of the universal affine-extension proof; all-field/mu1 theorem follows from written argument, not sampling',
 'external_review':False,'formal_verification':False}
args.output.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'passed','Z1':len(Z),'Zroot':r,'extended_covector_span':rank(cov,p),'mu1_char2_control':True}))
