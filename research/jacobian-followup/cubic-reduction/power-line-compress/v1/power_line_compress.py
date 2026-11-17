#!/usr/bin/env python3
"""Exact two-to-one cube compression along lines through the saved89 forms."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,combinations_with_replacement
from collections import defaultdict
import hashlib,json,time
from flint import fmpq_mat
import sympy as sp
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cubic(f):
 out={}
 for i,j,k in combinations_with_replacement([i for i,c in enumerate(f) if c],3):
  e=[0]*18
  for a in (i,j,k):e[a]+=1
  out[tuple(e)]=(1 if i==k else (3 if i==j or j==k else 6))*f[i]*f[j]*f[k]
 return out
def cross(f,g):
 a={};b={}
 for i,j,k in combinations_with_replacement([i for i in range(18) if f[i] or g[i]],3):
  e=[0]*18
  for z in (i,j,k):e[z]+=1
  mult=1 if i==k else (3 if i==j or j==k else 6)
  x=mult*(g[i]*f[j]*f[k]+f[i]*g[j]*f[k]+f[i]*f[j]*g[k])
  y=mult*(f[i]*g[j]*g[k]+g[i]*f[j]*g[k]+g[i]*g[j]*f[k])
  if x:a[tuple(e)]=x
  if y:b[tuple(e)]=y
 return a,b
def constraint(a,b,indices,t0=None):
 for e in indices:
  x=a.get(e,Q(0));y=b.get(e,Q(0))
  if not x and not y:continue
  if not x or not y:return False,None
  ratio=-Q(x)/Q(y)
  if t0 is not None and t0!=ratio:return False,None
  t0=ratio
 return True,t0
def main():
 out=HERE/"power-line-compress";out.mkdir(exist_ok=True);source=HERE/"power-basis-anneal/candidate.json";d=json.loads(source.read_text())
 forms=d["forms"];B=[[Q(c) for c in row] for row in d["B"]];n=len(forms);assert n==89
 cubes=[cubic(f) for f in forms];mons=sorted(set().union(*(set(p) for p in cubes)));mid={e:i for i,e in enumerate(mons)}
 C=fmpq_mat([[str(p.get(e,0)) for p in cubes] for e in mons])
 joint=fmpq_mat([[str(C[i,j]) for j in range(n)]+[str(int(i==j)) for j in range(len(mons))] for i in range(len(mons))])
 rr,rank=joint.rref();assert rank==len(mons)
 for i in range(n):assert rr[i,i]==1
 R=fmpq_mat([[str(rr[i,n+j]) for j in range(len(mons))] for i in range(len(mons))])
 groups=defaultdict(list)
 for j in range(n):
  v=tuple(row[j] for row in B);s=next(c for c in v if c)
  groups[tuple(c/s for c in v)].append(j)
 eligible=[(j,k) for key,ids in groups.items() if key[17]==0 for j,k in combinations(ids,2)]
 t=sp.Symbol("t");root_cache={};seen_lines=set();stats=dict(pairs=0,outside_rejected=0,residual_rejected=0,fixed_points=0,full_binary_lines=0,rational_trials=0,algebraic_trials=0,best=89)
 started=time.monotonic();answer=None
 def attempt(polys,modulus,f,g):
  nonlocal answer
  # Rational evaluation uses degree-one moduli; irreducible higher degrees represent algebraic roots.
  reduced=[sp.rem(p,modulus,t) for p in polys]
  for j,k in eligible:
   if reduced[j]==0:continue
   if sp.rem(reduced[j]*sp.Rational(str(B[0][k]))-reduced[k]*sp.Rational(str(B[0][j])),modulus,t)!=0 and B[0][j]:continue
   bj=[sp.Rational(str(row[j])) for row in B]
   canceled=[]
   for h in range(n):
    if all(sp.rem(reduced[j]*sp.Rational(str(B[o][h]))-bj[o]*reduced[h],modulus,t)==0 for o in range(18)):canceled.append(h)
   if len(canceled)<2 or d["padding_index"] in canceled:continue
   size=n-len(canceled)+1
   if size>=stats["best"]:continue
   stats["best"]=size
   answer=dict(status="compression-candidate",source_sha256=sha(Path(__file__)),input_sha256=sha(source),
     size=size,modulus=[str(c) for c in reversed(sp.Poly(modulus,t).all_coeffs())],
     first_form=f,second_form=g,reference_column=j,canceled=canceled,
     cube_coordinates=[[str(c) for c in reversed(sp.Poly(p,t).all_coeffs())] for p in reduced],
     scope="Exact primary two-to-one cube compression candidate. Requires full cube identities, Keller matrix and collision verification.")
   (out/"candidate.json").write_text(json.dumps(answer,indent=2,sort_keys=True)+"\n");print("BEST",size,"field degree",sp.degree(modulus,t),flush=True)
   return
 for p,q in combinations(range(n),2):
  if time.monotonic()-started>90:break
  stats["pairs"]+=1;a,b=cross(forms[p],forms[q]);outside=(set(a)|set(b))-set(mid)
  ok,fixed=constraint(a,b,outside)
  if not ok:stats["outside_rejected"]+=1;continue
  E=fmpq_mat([[str(a.get(e,0)),str(b.get(e,0))] for e in mons]);P=R*E
  aa={i:Q(str(P[i,0])) for i in range(len(mons)) if P[i,0]};bb={i:Q(str(P[i,1])) for i in range(len(mons)) if P[i,1]}
  ok,fixed=constraint(aa,bb,range(n,len(mons)),fixed)
  if not ok:stats["residual_rejected"]+=1;continue
  polys=[sp.Rational(str(aa.get(i,0)))*t+sp.Rational(str(bb.get(i,0)))*t**2+int(i==p)+int(i==q)*t**3 for i in range(n)]
  if fixed is not None:
   stats["fixed_points"]+=1;stats["rational_trials"]+=1
   attempt(polys,t-sp.Rational(str(fixed)),forms[p],forms[q])
  else:
   # The whole binary cubic space is contained in the saved cube span.
   line_key=tuple(tuple(str(c) for c in [rrr]) for rrr in [])
   # Deduplicate by the reduced row-space of the two source linear forms.
   lr,_=fmpq_mat([forms[p],forms[q]]).rref();line_key=str(lr)
   if line_key in seen_lines:continue
   seen_lines.add(line_key);stats["full_binary_lines"]+=1
   factors=set()
   for j,k in eligible:
    v=next(o for o in range(18) if B[o][j])
    ratio=sp.expand(sp.Rational(str(B[v][k]))*polys[j]-sp.Rational(str(B[v][j]))*polys[k])
    if ratio==0:
     factors.update([t-1,t+1,t-2]);continue
    if sp.degree(ratio,t)==0:continue
    for factor,_ in sp.factor_list(ratio,t)[1]:
     if factor==t:continue
     factors.add(sp.Poly(factor,t).monic().as_expr())
   for factor in factors:
    if sp.degree(factor,t)==1:stats["rational_trials"]+=1
    else:stats["algebraic_trials"]+=1
    attempt(polys,factor,forms[p],forms[q])
  if stats["pairs"]%100==0:
   (out/"checkpoint.json").write_text(json.dumps(dict(status="running",**stats),indent=2)+"\n")
 receipt=dict(status="candidate" if answer else ("completed-no-compression" if stats["pairs"]==n*(n-1)//2 else "bounded-no-compression"),
    source_sha256=sha(Path(__file__)),input_sha256=sha(source),**stats,seconds=round(time.monotonic()-started,3),
    scope="Lines through pairs of the89 saved forms; possible replacements lie in the saved89-cube span. Includes irreducible rational parameter factors of degree up to3. No full Veronese intersection or global optimality claim.")
 (out/"checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt),flush=True)
if __name__=="__main__":main()
