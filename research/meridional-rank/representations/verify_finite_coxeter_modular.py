#!/usr/bin/env python3
"""Independent prime-field replay of the H4 reflection-coloring exclusion.

No exact quadratic-field producer, diagram parser, saved propagation trace,
or retired campaign is imported. Root coordinates use only a standard
signed-permutation descriptor formula. The norm bound in the source note
connects two-prime rank deficiency to exact characteristic-zero deficiency.
"""
import hashlib,itertools,json,time
from pathlib import Path
D=Path(__file__).resolve().parent
# Root descriptors are only added/permuted/sign-changed. All arithmetic
# for reflections and determinants is done independently in prime fields.
z=(0,0)
base=[z,(1,0),(0,1),(-1,1)]
def minus(v):return tuple((-a,-b) for a,b in v)
def canonical(v):return min(v,minus(v))
formal=set()
for i in range(4):
 v=[z]*4;v[i]=(2,0);formal.add(canonical(tuple(v)))
for s in itertools.product((1,-1),repeat=4):formal.add(canonical(tuple((a,0) for a in s)))
for perm in itertools.permutations(range(4)):
 inversions=sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
 if inversions%2:continue
 for s in itertools.product((1,-1),repeat=3):
  v=[z,(s[0],0),(0,s[1]),(-s[2],s[2])]
  formal.add(canonical(tuple(v[i] for i in perm)))
formal=sorted(formal)
assert len(formal)==60

def project(v,p):
 first=next(a for a in v if a%p)
 inv=pow(first,-1,p)
 return tuple(a*inv%p for a in v)
def model(p,phi):
 assert (phi*phi-phi-1)%p==0
 vectors=[tuple((a+b*phi)%p for a,b in r) for r in formal]
 assert all(sum(a*a for a in v)%p==4 for v in vectors)
 ix={project(v,p):i for i,v in enumerate(vectors)}
 assert len(ix)==60
 q=[]
 half=pow(2,-1,p)
 for u in vectors:
  row=[]
  for v in vectors:
   coefficient=sum(a*b for a,b in zip(u,v))*half%p
   w=tuple((b-coefficient*a)%p for a,b in zip(u,v))
   row.append(ix[project(w,p)])
  q.append(row)
 assert all(q[a][q[a][b]]==b for a in range(60) for b in range(60))
 return vectors,q
# Dot products divided by2 remain in Z[phi]; this independent parity
# check uses the two coefficient formulas only modulo2, not exact producer
# arithmetic in the quadratic integer ring.
for u in formal:
 for v in formal:
  constant=sum(a*c+b*d for (a,b),(c,d) in zip(u,v))
  phi_coefficient=sum(a*d+b*c+b*d for (a,b),(c,d) in zip(u,v))
  assert constant%2==phi_coefficient%2==0
v11,q11=model(11,4);v31,q31=model(31,13)
assert q11==q31
try:
 model(11,5)
except AssertionError:
 invalid_phi_rejected=True
else:
 raise AssertionError('invalid quadratic specialization accepted')
Q=q11
# Normalize first seed to the axis-root reflection, with fresh label choice.
first=formal.index(canonical(((2,0),z,z,z)))
seen={first};queue=[first]
for a in queue:
 for s in range(60):
  b=Q[s][a]
  if b not in seen:seen.add(b);queue.append(b)
assert len(seen)==60

permutations=list(itertools.permutations(range(4)))
signs=[1 if sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))%2==0 else -1 for p in permutations]
def detmod(labels,vectors,p):
 total=0
 for perm,sign in zip(permutations,signs):
  term=sign
  for i,j in enumerate(perm):term=term*vectors[labels[i]][j]%p
  total=(total+term)%p
 return total

def diagram(pd,seeds):
 # Independent explicit component unions; no producer parser or trace.
 sets=[{v} for v in range(1,27)]
 for row in reversed(pd):
  a=next(c for c in sets if row[1] in c);b=next(c for c in sets if row[3] in c)
  if a is not b:sets.remove(a);sets.remove(b);sets.append(a|b)
 sets=sorted(sorted(c) for c in sets)
 arc={v:i for i,s in enumerate(sets) for v in s}
 crossings=[(arc[r[1]],arc[r[0]],arc[r[2]]) for r in pd]
 known=set(seeds);plan=[]
 while len(known)<13:
  step=None
  for o,a,b in reversed(crossings):
   if o in known and b in known and a not in known:step=(o,b,a);break
   if o in known and a in known and b not in known:step=(o,a,b);break
  assert step
  plan.append(step);known.add(step[2])
 return crossings,plan

# Positive control: the four standard simple roots have independent
# images and their reflection permutations have the H4 Coxeter orders.
simple=[((-2,0),z,z,z),((0,1),(-1,0),(-1,1),z),
        (z,(0,1),(1,-1),(-1,0)),(z,z,z,(2,0))]
simple_labels=[formal.index(canonical(r)) for r in simple]
assert detmod(simple_labels,v11,11) and detmod(simple_labels,v31,31)
def permutation_order(a,b):
 p=tuple(Q[a][Q[b][i]] for i in range(60))
 seen=set();order=1
 from math import lcm
 for start in range(60):
  if start in seen:continue
  v=start;length=0
  while v not in seen:seen.add(v);length+=1;v=p[v]
  order=lcm(order,length)
 return order
orders=[permutation_order(simple_labels[i],simple_labels[j])
        for i in range(4) for j in range(i+1,4)]
assert orders==[5,2,2,3,2,3]
seen={simple_labels[0]};queue=[simple_labels[0]]
for a in queue:
 for s in simple_labels:
  b=Q[s][a]
  if b not in seen:seen.add(b);queue.append(b)
assert len(seen)==60
def symmetric_quandle(n):
 # Independent full permutation multiplication, rather than endpoint swaps.
 transpositions=[]
 for a,b in itertools.combinations(range(n),2):
  p=list(range(n));p[a],p[b]=p[b],p[a];transpositions.append(tuple(p))
 index={p:i for i,p in enumerate(transpositions)}
 return [[index[tuple(a[b[a[v]]] for v in range(n))]
          for b in transpositions] for a in transpositions]

def constant_only(seeds,crossings,plan,q):
 count=0
 for tail in itertools.product(range(len(q)),repeat=3):
  labels=[-1]*13
  for a,b in zip(seeds,(0,)+tail):labels[a]=b
  for o,a,b in plan:labels[b]=q[labels[o]][labels[a]]
  if all(q[labels[o]][labels[a]]==labels[b] for o,a,b in crossings):
   count+=1
   assert all(a==0 for a in labels)
 return count

symmetric={n:symmetric_quandle(n) for n in (3,4,5)}
inputs=json.loads((D/'seven_sheet_cover_inputs.json').read_text())
rows=[];start=time.monotonic()
for case in inputs['cases']:
 if not case['actions'][0]['expected_H1_F3_dimension']:continue
 seeds=case['seed_arcs'];crossings,plan=diagram(case['pd'],seeds)
 accepted=mod11_full=mod31_full=0
 symmetric_counts={n:constant_only(seeds,crossings,plan,q) for n,q in symmetric.items()}
 assert all(count==1 for count in symmetric_counts.values())
 for tail in itertools.product(range(60),repeat=3):
  labels=[-1]*13
  for a,b in zip(seeds,(first,)+tail):labels[a]=b
  for o,a,b in plan:labels[b]=Q[labels[o]][labels[a]]
  if any(Q[labels[o]][labels[a]]!=labels[b] for o,a,b in crossings):continue
  accepted+=1
  basis=(first,)+tail
  mod11_full+=detmod(basis,v11,11)!=0
  mod31_full+=detmod(basis,v31,31)!=0
 assert mod11_full==mod31_full==0
 row={'name':case['name'],'trials':60**3,'accepted':accepted,'rank4_mod11':mod11_full,'rank4_mod31':mod31_full,'derived_plan':plan,'symmetric_constant_only_counts':symmetric_counts}
 rows.append(row)
 print(row['name'],accepted,'no rank4 at11/31','seconds',round(time.monotonic()-start,2),flush=True)
report={'status':'PASS','arithmetic':'Independent prime-field reflections/determinants at11(phi4),31(phi13). Root descriptor generation shared mathematical formula only. No producer arithmetic, parser or propagation trace imported.','normalization_first_root':first,'results':rows,'exact_zero_argument':'Both real embeddings of every scaled root have Euclidean norm2. Hadamard gives determinant embeddings abs≤16 and integer norm abs≤256. A determinant zero modulo primes above11 and31 has integer norm divisible341>256, forcing determinant0.','simple_root_positive_control':{'labels':simple_labels,'Coxeter_orders':orders,'root_line_orbit':60},'root_line_injectivity':True,'quandles_agree_at_both_primes':True,'invalid_phi_rejected':invalid_phi_rejected}
assert len(rows)==24
paths=['verify_finite_coxeter_modular.py','seven_sheet_cover_inputs.json']
report['input_sha256']={p:hashlib.sha256((D/p).read_bytes()).hexdigest() for p in paths}
(D/'finite_coxeter_modular_verification.json').write_text(json.dumps(report,indent=2)+'\n')
