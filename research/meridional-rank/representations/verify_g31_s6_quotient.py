#!/usr/bin/env python3
"""Exact stdlib replay of G31 -> S6 from rational Gaussian matrices."""
from pathlib import Path
from fractions import Fraction
from itertools import combinations,product
from collections import deque
import json,hashlib,time,sys
HERE=Path(__file__).resolve().parent
INPUT=HERE/"g31_s6_quotient_inputs.json"
RECEIPT=HERE/"g31_s6_quotient_verification.json"
def gm(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def gs(a,b):return (a[0]-b[0],a[1]-b[1])
def ga(a,b):return (a[0]+b[0],a[1]+b[1])
def conj(a):return (a[0],-a[1])
def flatten(a):return tuple(x for row in a for z in row for x in z)
def matrows(a):return [[(a[8*i+2*j],a[8*i+2*j+1]) for j in range(4)] for i in range(4)]
ID=tuple(x for i in range(4) for j in range(4) for x in (2*int(i==j),0))
def mm(a,b):
 out=[]
 for i in range(4):
  for j in range(4):
   r=s=0
   for k in range(4):
    ar,ai=a[8*i+2*k],a[8*i+2*k+1];br,bi=b[8*k+2*j],b[8*k+2*j+1]
    r+=ar*br-ai*bi;s+=ar*bi+ai*br
   assert r%2==s%2==0,"Product outside retained half-integral Gaussian encoding"
   out.extend((r//2,s//2))
 return tuple(out)
def star(a):return tuple(x for i in range(4) for j in range(4) for x in (a[8*j+2*i],-a[8*j+2*i+1]))
def reflection(u):
 norm=sum(a*a+b*b for a,b in u);assert norm
 out=[]
 for i in range(4):
  row=[]
  for j in range(4):
   a,b=gm(u[i],conj(u[j]));assert (4*a)%norm==(4*b)%norm==0
   row.append((2*int(i==j)-4*a//norm,-4*b//norm))
  out.append(row)
 return flatten(out)
def normalized(u):
 a=next(z for z in u if z!=(0,0));norm=a[0]*a[0]+a[1]*a[1]
 return tuple((Fraction(z[0]*a[0]+z[1]*a[1],norm),Fraction(z[1]*a[0]-z[0]*a[1],norm)) for z in u)
def mv_twice(a,u):
 return tuple((sum(a[8*i+2*j]*u[j][0]-a[8*i+2*j+1]*u[j][1] for j in range(4)),
 sum(a[8*i+2*j]*u[j][1]+a[8*i+2*j+1]*u[j][0] for j in range(4))) for i in range(4))
def perm_closure(gens,n):
 identity=tuple(range(n));seen={identity};queue=[identity]
 for a in queue:
  for b in gens:
   c=tuple(a[b[i]] for i in range(n))
   if c not in seen:seen.add(c);queue.append(c)
 return seen
def rankone_displacement(a):
 d=list(a)
 for i in range(4):d[8*i+2*i]-=2
 d=matrows(d)
 pos=next(((i,j) for i in range(4) for j in range(4) if d[i][j]!=(0,0)),None)
 if pos is None:return False
 i,j=pos;v=d[i][j]
 return all(gm(v,d[k][l])==gm(d[k][j],d[i][l]) for k in range(4) for l in range(4))
def mul_right_special(a,columns):
 out=[]
 for i in range(4):
  for col in columns:
   r=s=0
   for k,br,bi in col:
    ar,ai=a[8*i+2*k],a[8*i+2*k+1]
    r+=ar*br-ai*bi;s+=ar*bi+ai*br
   assert r%2==s%2==0
   out.extend((r//2,s//2))
 return tuple(out)
def full_closure(gens):
 columns=[[[ (k,g[8*k+2*j],g[8*k+2*j+1]) for k in range(4)
             if g[8*k+2*j] or g[8*k+2*j+1]] for j in range(4)] for g in gens]
 seen={ID};queue=[ID];level=[1];depth=[0]
 for k,a in enumerate(queue):
  for cols in columns:
   b=mul_right_special(a,cols)
   if b not in seen:
    seen.add(b);queue.append(b);d=depth[k]+1;depth.append(d)
    if d==len(level):level.append(0)
    level[d]+=1
  assert len(seen)<=46080
 return seen,level
def main():
 start=time.monotonic();data=json.loads(INPUT.read_text())
 sources={}
 for item in data["primary_sources"]:
  p=HERE/item["local_path"];raw=p.read_bytes()
  assert hashlib.sha256(raw).hexdigest()==item["sha256"] and len(raw)==item["bytes"]
  sources[item["local_path"]]=item["sha256"]
 roots=[tuple(tuple(z) for z in u) for u in data["roots_gaussian_integer_pairs"]]
 assert len(roots)==60 and len(set(map(normalized,roots)))==60
 # Independently regenerate the displayed root-line family in its exact order.
 expected=[];units=[(1,0),(0,1),(-1,0),(0,-1)]
 for i in range(4):expected.append(tuple((int(j==i),0) for j in range(4)))
 for i in range(4):
  for j in range(i+1,4):
   for k in range(4):
    u=[(0,0)]*4;u[i]=(1,0);u[j]=(-units[k][0],-units[k][1]);expected.append(tuple(u))
 for a in range(4):
  for b in range(4):
   for c in range(4):
    if (a+b+c)%2==0:expected.append(((1,0),units[a],units[b],units[c]))
 assert roots==expected
 monomial_orbit={tuple(units[k] for k in exponents) for exponents in product(range(4),repeat=4) if sum(exponents)%2==0}
 assert len(monomial_orbit)==128
 assert {normalized(u) for u in monomial_orbit}=={normalized(u) for u in roots[28:]}
 assert len({normalized(u) for u in monomial_orbit})==32

 refs=[reflection(u) for u in roots]
 assert refs==[flatten(a) for a in data["twice_reflection_matrices"]]
 assert len(set(refs))==60
 for a in refs:assert mm(a,a)==ID and mm(star(a),a)==ID and rankone_displacement(a)
 rootid={normalized(u):i for i,u in enumerate(roots)}
 root_actions=[[rootid[normalized(mv_twice(a,u))] for u in roots] for a in refs]
 Q=data["conjugation_table"];assert Q==root_actions
 for a in Q:assert sorted(a)==list(range(60))
 generators=data["generator_indices"];assert generators==[3,4,14,37,28]
 assert [roots[i] for i in generators]==[tuple(tuple(z) for z in u) for u in data["published_generator_roots"]]
 gens=[refs[i] for i in generators]
 assert data["generator_root_line_actions"]==[Q[i] for i in generators]
 blocks=data["blocks"];assert len(blocks)==15 and all(len(b)==4 for b in blocks)
 assert sorted(sum(blocks,[]))==list(range(60))
 blockof={x:i for i,b in enumerate(blocks) for x in b}
 labels=[tuple(p) for p in data["block_labels_as_pairs"]]
 pairs=list(combinations(range(1,7),2));assert sorted(labels)==pairs
 f=[labels[blockof[x]] for x in range(60)]
 assert set(f)==set(pairs) and all(f.count(p)==4 for p in pairs)
 def swap(pair,x):
  a,b=pair;return b if x==a else a if x==b else x
 def target_action(a,b):return tuple(sorted(swap(a,x) for x in b))
 blockactions=[]
 for i in range(60):
  action=[]
  for b in blocks:
   image={blockof[Q[i][x]] for x in b};assert len(image)==1;action.append(image.pop())
  assert all(labels[action[j]]==target_action(f[i],labels[j]) for j in range(15))
  blockactions.append(tuple(action))
 assert all(f[Q[i][j]]==target_action(f[i],f[j]) for i in range(60) for j in range(60))
 targetgens=[tuple(swap(f[i],x)-1 for x in range(1,7)) for i in generators]
 target=perm_closure(targetgens,6);assert len(target)==720
 pairkernel=[g for g in target if all(tuple(sorted(g[x-1]+1 for x in p))==p for p in pairs)]
 assert pairkernel==[tuple(range(6))]
 assert len(perm_closure([blockactions[i] for i in generators],15))==720
 projective=perm_closure([tuple(Q[i]) for i in generators],60);assert len(projective)==11520
 matrices,levels=full_closure(gens);assert len(matrices)==46080 and all(a in matrices for a in refs)
 scalar_matrices={tuple(x for i in range(4) for j in range(4) for x in ((a if i==j else 0),(b if i==j else 0))) for a,b in [(2,0),(-2,0),(0,2),(0,-2)]}
 assert len(scalar_matrices)==4 and scalar_matrices<=matrices
 assert len(matrices)//len(projective)==len(scalar_matrices)

 allrefs={a for a in matrices if rankone_displacement(a)};assert allrefs==set(refs)
 seen={0};queue=[0]
 for a in queue:
  for i in generators:
   b=Q[i][a]
   if b not in seen:seen.add(b);queue.append(b)
 assert len(seen)==60
 assert len(matrices)//len(target)==64
 # Full closure checks the displayed five matrices, not only their projective action.
 out={"status":"PASS exact G31 matrix-group quotient onto S6",
 "scope":"Genuine group homomorphism from a verified exact matrix action on60 root lines and15 invariant blocks; no knot count or diagram identification asserted.",
 "matrix_group_order":len(matrices),"root_line_action_order":len(projective),"scalar_kernel_order":len(matrices)//len(projective),
 "all_group_reflections":len(allrefs),"G424_vector_orbit_size":len(monomial_orbit),"G424_projective_orbit_size":32,"single_reflection_conjugacy_class":len(seen),"target_order":len(target),
 "target_pair_action_kernel_order":len(pairkernel),"full_quotient_kernel_order":len(matrices)//len(target),
 "reflection_preimages_per_transposition":4,"all_reflection_block_action_checks":60*15,"all_quandle_map_checks":3600,
 "generator_images_as_transpositions":[f[i] for i in generators],"minimum_reflection_class_generating_number":5,
 "matrix_BFS_level_counts":levels,"elapsed_seconds":round(time.monotonic()-start,3),
 "python_version":sys.version.split()[0],"primary_source_sha256":sources,
 "input_sha256":hashlib.sha256(INPUT.read_bytes()).hexdigest(),"checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 "proof_source_sha256":hashlib.sha256((HERE/"G31_COXETER_QUOTIENT.tex").read_bytes()).hexdigest(),
 "assurance":"Exact internal standard-library replay; written homomorphism, class-rank and KCH consequences in the note; not formal verification or external review."}
 # Keep dated review records attached to their own source/PDF hashes.
 # A replay never upgrades an old review to the current proof edition.
 for attempt in range(5):
  previous=RECEIPT.read_bytes() if RECEIPT.exists() else None
  saved=json.loads(previous) if previous is not None else {}
  merged=dict(out)
  if "internal_pdf_reviews" in saved:merged["internal_pdf_reviews"]=saved["internal_pdf_reviews"]
  temporary=RECEIPT.with_name(RECEIPT.name+".tmp")
  temporary.write_text(json.dumps(merged,indent=2)+"\n")
  current=RECEIPT.read_bytes() if RECEIPT.exists() else None
  if current==previous:temporary.replace(RECEIPT);break
 else:raise RuntimeError("Concurrent receipt changes; no review metadata overwritten")
 print(json.dumps(out,indent=2),flush=True)
if __name__=="__main__":main()
