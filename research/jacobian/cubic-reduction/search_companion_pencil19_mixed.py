#!/usr/bin/env python3
"""Exact shared-cube basis exchanges optimizing a nonzero rational pencil parameter."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import hashlib,json,random,time
from flint import fmpq_mat
from prepare_companion_pencil19 import cubes
HERE=Path(__file__).resolve().parent
OUT=HERE/'companion-pencil19-mixed'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def plus(a,b,c):
 out=dict(a)
 for j,v in b.items():out[j]=out.get(j,Q(0))+c*v
 return {j:v for j,v in out.items() if v}
def pivot(A,B,basis,p,q):
 ap=A[p][q];ar={j:c/ap for j,c in A[p].items()};br={j:c/ap for j,c in B[p].items()}
 for k in range(len(A)):
  if k!=p and A[k].get(q):
   value=A[k][q];A[k]=plus(A[k],ar,-value)
   if br:B[k]=plus(B[k],br,-value)
 A[p]=ar;B[p]=br;basis[p]=q

def signature_groups(B):
 groups={};scalars={}
 for p,b in enumerate(B):
  if b:
   scalar=b[min(b)];groups.setdefault(tuple((j,c/scalar) for j,c in sorted(b.items())),[]).append(p);scalars[p]=scalar
 return {p:g for g in groups.values() for p in g},scalars

def best_parameter(B):
 roots=Counter();generic=sum(bool(b) for b in B)
 for b in B:
  ids=[j for j in range(19) if b.get(j+19)]
  if not ids:continue
  j=ids[0];root=-b.get(j,Q(0))/b[j+19]
  if root and all(b.get(j,Q(0))+root*b.get(j+19,Q(0))==0 for j in range(19)):roots[root]+=1
 if not roots:return generic,Q(-3,2),generic,{}
 root=min(roots,key=lambda r:(-roots[r],abs(r.numerator)+r.denominator,r));return generic-roots[root],root,generic,roots

def main():
 ip=OUT/'input.json';bp=OUT/'basis.json';d=json.loads(ip.read_text());bd=json.loads(bp.read_text());assert bd['input_sha256']==sha(ip);initialA=[{int(k):Q(c) for k,c in row.items()} for row in bd['A']];initialB=[{int(k):Q(c) for k,c in row.items()} for row in bd['B']];initialbasis=bd['basis'];forms=list(map(tuple,d['forms']));rank=len(initialbasis);m=len(forms);padding=d['padding_index'];H0=[{tuple(e):Q(c) for e,c in row} for row in d['H0']];H1=[{tuple(e):Q(c) for e,c in row} for row in d['H1']];mons=sorted(set().union(*(set(p) for p in H0+H1)));polys=[cubes(f) for f in forms];mons=sorted(set(mons).union(*(set(p) for p in polys)));C=fmpq_mat([[str(p.get(e,0)) for p in polys] for e in mons]);best,sigma,_,_=best_parameter(initialB);rng=random.Random(261121);started=time.monotonic();log=[];steps=0;seen=set();warm=(initialA,initialB,initialbasis,[])
 def save(B,basis,sigma,history):
  rows={basis[k]:{j:v for j in range(19) if (v:=B[k].get(j,Q(0))+sigma*B[k].get(j+19,Q(0)))} for k in range(rank)};rows={k:v for k,v in rows.items() if v};selected=sorted(rows);assert padding in rows and [rows[k].get(18,Q(0)) for k in selected]==[Q(k==padding) for k in selected];mat=fmpq_mat([[str(rows[k].get(j,Q(0))) for j in range(19)] for k in selected]);cc=fmpq_mat([[str(C[i,k]) for k in selected] for i in range(len(mons))]);T=fmpq_mat([[str(a.get(e,Q(0))+sigma*b.get(e,Q(0))) for a,b in zip(H0,H1)] for e in mons]);assert cc*mat==T
  candidate=dict(status='exact-primary-shared-cube-identity',source_sha256=sha(Path(__file__)),input_sha256=sha(ip),sigma=str(sigma),size=len(selected),forms=[forms[k] for k in selected],B=[[str(rows[k].get(j,Q(0))) for k in selected] for j in range(19)],padding_index=selected.index(padding),pivots=history,scope='Exact identity for the changed19-variable companion at nonzero rational sigma. This is not an improved Keller witness until independent base/shear/companion/collision and full cubic-linear matrix verification.');(OUT/'best-identity.json').write_text(json.dumps(candidate,indent=2,sort_keys=True)+'\n');return len(selected)
 assert save(initialB,initialbasis,sigma,[])==best;print('INITIAL',best,'sigma',sigma,flush=True);stop=False
 for restart in range(12):
  initialA,initialB,initialbasis,history0=warm;A=[dict(row) for row in initialA];B=[dict(row) for row in initialB];basis=list(initialbasis);history=list(history0)
  for step in range(1500):
   if time.monotonic()-started>40:stop=True;break
   active=[k for k,b in enumerate(B) if b];zeros=[k for k,b in enumerate(B) if not b];generic=len(active)
   if step%3==0 and zeros:
    possible=[(p,q) for p in zeros for q in A[p] if q not in basis]
    if possible:
     p,q=rng.choice(possible);pivot(A,B,basis,p,q);history.append([p,q])
   groups,scalars=signature_groups(B);basis_set=set(basis);current={basis[k] for k in active};options=[];mindelta=999
   for q in range(m):
    if q in basis_set:continue
    added=[k for k in zeros if A[k].get(q)]
    for p in active:
     ap=A[p].get(q)
     if not ap or basis[p]==padding:continue
     removed=[k for k in groups[p] if k!=p and A[k].get(q,Q(0))*scalars[p]==ap*scalars[k]];delta=len(added)-len(removed)
     if generic+delta>best+7:continue
     key=tuple(sorted((current-{basis[p]}-{basis[k] for k in removed})|{q}|{basis[k] for k in added}))
     if key in seen:continue
     options.append((delta,p,q,key));mindelta=min(mindelta,delta)
   if not options:break
   delta_group=mindelta if rng.random()<0.75 else rng.choice(sorted({o[0] for o in options if o[0]<=mindelta+2}));delta,p,q,key=rng.choice([o for o in options if o[0]==delta_group]);pivot(A,B,basis,p,q);history.append([p,q]);steps+=1;seen.add(key);value,sigma,newgeneric,roots=best_parameter(B);assert newgeneric==generic+delta
   if value<best:
    best=value;assert save(B,basis,sigma,history)==best;warm=([dict(row) for row in A],[dict(row) for row in B],list(basis),list(history));print('BEST',best,'sigma',sigma,'step',steps,'seconds',round(time.monotonic()-started,2),flush=True)
    if best<89:stop=True;break
   if steps%10==0:
    log.append(dict(restart=restart,step=step,generic=newgeneric,pencil_support=value,sigma=str(sigma),best=best,seconds=round(time.monotonic()-started,3)));receipt=dict(status='running',source_sha256=sha(Path(__file__)),input_sha256=sha(ip),basis_sha256=sha(bp),steps=steps,best=best,seed=261121,time_budget_seconds=40,log=log,scope='Bounded exact basis exchanges in427 cubes, simultaneous optimization over nonzero rational pencil parameters. No finite-pool or global minimum.');(OUT/'search-checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
  if stop:break
 receipt=dict(status='candidate-below89-needs-independent-witness' if best<89 else 'bounded-no-improvement-below89',source_sha256=sha(Path(__file__)),input_sha256=sha(ip),basis_sha256=sha(bp),steps=steps,best=best,seed=261121,time_budget_seconds=40,seconds=round(time.monotonic()-started,3),log=log,best_identity_sha256=sha(OUT/'best-identity.json'),scope='Bounded exact basis exchanges in427 cubes, simultaneous optimization over nonzero rational pencil parameters. Exact saved primary cube identity but no finite-pool/global minimum or smaller Keller witness unless separately verified.');(OUT/'search-checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print('DONE',best,'STEPS',steps,'SECONDS',receipt['seconds'],flush=True)
if __name__=='__main__':main()
