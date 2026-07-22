#!/usr/bin/env python3
from pathlib import Path
import pickle,re,sympy as sp
from flint import fmpq,fmpq_poly
import exact_core as ec
R=Path(__file__).resolve().parent;st=pickle.loads((R/'case1_checkpoint.pkl').read_bytes())
def dk(v):return ec.K(fmpq_poly([fmpq(n,d) for n,d in v]))
def dp(v):return ec.PP({tuple(m):dk(c) for m,c in v})
qs=[dp(q) for q in st['all_comp']]
def pp_pow(a,n):
 o=ec.PP.const(1)
 while n:
  if n&1:o=o*a
  a=a*a;n//=2
 return o
def sub_scalar(q,idx,val):
 z={}
 for m,c in q.d.items():
  mm=list(m);e=mm[idx];mm[idx]=0;mm=tuple(mm);z[mm]=z.get(mm,ec.K(0))+c*val**e
 return ec.PP(z)
def sub_pp(q,idx,val):
 o=ec.PP()
 for m,c in q.d.items():
  mm=list(m);e=mm[idx];mm[idx]=0;o+=ec.PP({tuple(mm):c})*pp_pow(val,e)
 return o
def shape(q):return (len(q.d),max((sum(m) for m in q.d),default=-1),[max((m[i] for m in q.d),default=0) for i in range(6)])
def coeff_m(q,m):return q.d.get(tuple(m),ec.K(0))
def ratio(a,b):
 if set(a.d)!=set(b.d):return None
 m=next(iter(a.d));rr=a.d[m]/b.d[m]
 return rr if all(a.d[k]==b.d[k]*rr for k in a.d) else None

def coeff_var(q,idx):
 z={}
 for m,c in q.d.items():
  if m[idx]==1:
   mm=list(m);mm[idx]=0;z[tuple(mm)]=c
 return ec.PP(z)

fl=(R/'factor_q0_exact_clean.out').read_text().splitlines();Ssym=sp.Symbol('s')
for br,line in enumerate(fl[:2],1):
 raw=line.split('=',1)[1];raw=re.sub(r'u(\d+)',r'u**\1',raw);raw=re.sub(r'(\d)u',r'\1*u',raw);expr=sp.sympify(raw,locals={'s':Ssym,'u':ec.U});sv=ec.eval_u(-sp.expand(expr-Ssym))
 red=[sub_scalar(q,1,sv) for q in qs]
 # Eliminate u1/u3 between E6 and E2; the remainder is a square.
 lam=ratio(coeff_var(red[6],3),coeff_var(red[2],3)) or ratio(coeff_var(red[6],5),coeff_var(red[2],5)); assert lam
 elim=red[6]-red[2]*lam
 print('lambda E6/E2',lam.sing())
 S=coeff_m(elim,(2,0,0,0,0,0));assert S
 # L=r+alpha*h+beta*u2+gamma using r-linear coefficients
 alpha=coeff_m(elim,(1,0,1,0,0,0))/(2*S)
 beta=coeff_m(elim,(1,0,0,0,1,0))/(2*S)
 gamma=coeff_m(elim,(1,0,0,0,0,0))/(2*S)
 L=ec.PP.var(0)+ec.PP.const(alpha)*ec.PP.var(2)+ec.PP.const(beta)*ec.PP.var(4)+ec.PP.const(gamma)
 assert not (elim-L*L*S)
 rsol=-(ec.PP.const(alpha)*ec.PP.var(2)+ec.PP.const(beta)*ec.PP.var(4)+ec.PP.const(gamma))
 print('\nBRANCH',br,'square alpha',alpha.sing(),'beta',beta.sing(),'gamma udeg',gamma.p.degree())
 rr=[sub_pp(q,0,rsol) for q in red]
 print('survivors after r',[(i,shape(q)) for i,q in enumerate(rr) if q])
 # constant K-linear relations among survivors
 inds=[i for i,q in enumerate(rr) if q]
 for i0,i in enumerate(inds):
  for j in inds[i0+1:]:
   x=ratio(rr[i],rr[j])
   if x:print(' proportional',i,j,x.sing())
 # Save reduced exact equations (dedupe proportional greedily)
 keep=[]
 for i in inds:
  if not any(ratio(rr[i],rr[j]) for j in keep):keep.append(i)
 print('keep',keep)
 (R/f'case1_branch{br}_after_r_meta.txt').write_text('s='+sv.sing()+'\nr='+rsol.sing(6)+'\nkeep='+str(keep)+'\n')
 (R/f'case1_branch{br}_after_r_eqs.txt').write_text('\n'.join(rr[i].sing(6) for i in keep))
