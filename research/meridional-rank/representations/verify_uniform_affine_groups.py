#!/usr/bin/env python3
"""Exact integral certificate: no affine-group enumeration or CAS dependency."""
import hashlib,json
from pathlib import Path
BASE=Path(__file__).resolve().parent
def compose(a,b): return tuple(a[b[i]] for i in range(7))
def permute(a,v):
 out=[0]*7
 for i,x in enumerate(v): out[a[i]]=x
 return tuple(out)
def multiply(a,b):
 return compose(a[0],b[0]),tuple(x+y for x,y in zip(a[1],permute(a[0],b[1])))
def inverse(a):
 p=tuple(a[0].index(i) for i in range(7))
 return p,tuple(-x for x in permute(p,a[1]))
UNIT=(tuple(range(7)),(0,)*7)
def evaluate(word,generators):
 result=UNIT
 for i in word:
  g=generators[abs(i)-1]
  result=multiply(result,g if i>0 else inverse(g))
 return result

def determinant(rows):
 a=[list(row) for row in rows];sgn=1;den=1
 for k in range(len(a)-1):
  if not a[k][k]:
   j=next((j for j in range(k+1,len(a)) if a[j][k]),None)
   if j is None: return 0
   a[k],a[j]=a[j],a[k];sgn=-sgn
  pivot=a[k][k]
  for i in range(k+1,len(a)):
   for j in range(k+1,len(a)):
    numerator=a[i][j]*pivot-a[i][k]*a[k][j]
    assert numerator%den==0
    a[i][j]=numerator//den
  den=pivot
 return sgn*a[-1][-1]

cert=json.loads((BASE/'uniform_affine_group_inputs.json').read_text())
x=[(tuple(p),tuple(v)) for p,v in zip(cert['quartet_permutations'],cert['quartet_translations'])]
assert len(x)==4
for p,v in x:
 assert compose(p,p)==UNIT[0] and sum(v)==0
 assert tuple(i+j for i,j in zip(v,permute(p,v)))==UNIT[1]
# The nontrivial translate is explicitly conjugate to the split involution.
potential=(0,1,-1,0,0,0,0)
assert x[0][1]==tuple(a-b for a,b in zip(potential,permute(x[0][0],potential)))
assert all(not any(z[1]) for z in x[1:])
a,b=[evaluate([i+1 for i in w],x) for w in cert['ordinary_generators_as_quartet_words']]
seen={UNIT[0]};queue=[UNIT[0]]
for p in queue:
 for q in (a[0],b[0]):
  r=compose(p,q)
  if r not in seen:seen.add(r);queue.append(r)
assert len(seen)==168
# Every enumerated permutation is induced by a binary linear map.
for p in queue:
 for i in range(1,8):
  for j in range(1,8):
   if i!=j:assert (p[i-1]+1)^(p[j-1]+1)==p[(i^j)-1]+1
translation=evaluate(cert['translation_word'],[a,b])
assert translation==(UNIT[0],tuple(cert['translation']))
rows=[]
for item in cert['unimodular_conjugates']:
 g=evaluate(item['conjugator'],[a,b])
 t=multiply(multiply(g,translation),inverse(g))
 assert t==(UNIT[0],tuple(item['translation']))
 assert sum(t[1])==0
 rows.append(t[1][:6])
det=determinant(rows)
assert abs(det)==1
corrupt=[list(r) for r in rows];corrupt[0]=corrupt[1][:]
assert determinant(corrupt)==0,'duplicated translation basis not detected'
input_names=['TERNARY_AFFINE_QUOTIENTS.tex','verify_uniform_affine_groups.py','uniform_affine_group_inputs.json']
receipt={'input_sha256':{n:hashlib.sha256((BASE/n).read_bytes()).hexdigest() for n in input_names},'status':'PASS','arithmetic':'integers','projected_order':len(queue),'translation':translation[1],'basis_determinant':det,'ordinary_generating_rank':2,'distinguished_involution_class_rank':4,'scope':'All primes p and all integer moduli m>=2; class-rank lower bound additionally uses the written fixed-space argument.'}
(BASE/'uniform_affine_group_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
