#!/usr/bin/env python3
"""Add mixed product-pencil cubes beyond the original changed-tensor pool."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
from flint import fmpq_mat
from prepare_companion_pencil19 import cubes,decompose,plus,scaled
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();parent=HERE/'companion-pencil19';ip=parent/'input.json';d=json.loads(ip.read_text());candidate=json.loads((parent/'best-identity.json').read_text());assert candidate['input_sha256']==sha(ip);oldforms=list(map(tuple,d['forms']));pool=set(oldforms);unit=[tuple(int(i==j) for i in range(19)) for j in range(19)];a,b,c,R,Tau,S=[unit[j] for j in [0,1,2,3,18,6]];products=[]
 for k in [-3,-2,-1,1,2,3]:
  products += [[b,c,plus(Tau,scaled(a,k))],[b,R,plus(Tau,scaled(S,k))],[b,plus(c,scaled(R,k)),Tau]]
 for fs in products:
  for mode in range(5):pool.update(decompose(fs,mode))
 pool.discard((0,)*19);forms=oldforms+sorted(pool-set(oldforms));positions={f:i for i,f in enumerate(forms)};seed=[positions[tuple(f)] for f in candidate['forms']];priority=seed+[i for i in range(len(forms)) if i not in seed];P0=[{tuple(e):Q(c) for e,c in p} for p in d['H0']];P1=[{tuple(e):Q(c) for e,c in p} for p in d['H1']];polys=[cubes(f) for f in forms];mons=sorted(set().union(*(set(p) for p in polys+P0+P1)));m=len(forms);C=fmpq_mat([[str(p.get(e,0)) for p in polys] for e in mons]);T=fmpq_mat([[str(p.get(e,0)) for p in P0+P1] for e in mons]);ordered=fmpq_mat([[str(C[i,j]) for j in priority] for i in range(len(mons))]);rr,rank=ordered.rref();basis=[priority[next(j for j in range(m) if rr[i,j])] for i in range(rank)];assert basis[:len(seed)]==seed;matrix=fmpq_mat([[str(C[i,j]) for j in basis]+[str(C[i,j]) for j in range(m)]+[str(T[i,j]) for j in range(38)] for i in range(len(mons))]);rr,r=matrix.rref();assert r==rank;A=[{str(j):str(rr[i,rank+j]) for j in range(m) if rr[i,rank+j]} for i in range(rank)];B=[{str(j):str(rr[i,rank+m+j]) for j in range(38) if rr[i,rank+m+j]} for i in range(rank)];out=HERE/'companion-pencil19-mixed';out.mkdir(exist_ok=True);payload=dict(d);payload.update(source_sha256=sha(Path(__file__)),parent_input_sha256=sha(ip),seed_identity_sha256=sha(parent/'best-identity.json'),forms=forms,rank=rank,monomials=len(mons),mixed_product_bases=[[list(v) for v in fs] for fs in products],scope='Same changed19-variable companion pencil; old277 forms enlarged with cube polarizations of bc*(tau+k*a), bR*(tau+k*S), b*(c+kR)*tau for k=+-1,+-2,+-3. Not a global minimum.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');(out/'basis.json').write_text(json.dumps(dict(input_sha256=sha(out/'input.json'),basis=basis,A=A,B=B),indent=2,sort_keys=True)+'\n');print('MIXED POOL',m,'RANK',rank,'MONOMIALS',len(mons),'SECONDS',round(time.monotonic()-start,3),flush=True)
if __name__=='__main__':main()
