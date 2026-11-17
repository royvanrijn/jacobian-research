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
def coeff_var(q,idx):
 z={}
 for m,c in q.d.items():
  if m[idx]==1:
   mm=list(m);mm[idx]=0;z[tuple(mm)]=c
 return ec.PP(z)
def ratio(a,b):
 if set(a.d)!=set(b.d) or not a.d:return None
 m=next(iter(a.d));rr=a.d[m]/b.d[m]
 return rr if all(a.d[k]==b.d[k]*rr for k in a.d) else None
def cm(q,m):return q.d.get(tuple(m),ec.K(0))
def shape(q):return(len(q.d),max((sum(m) for m in q.d),default=-1),[max((m[i] for m in q.d),default=0) for i in range(6)])
fl=(R/'factor_q0_exact_clean.out').read_text().splitlines();ss=sp.Symbol('s')
for br,line in enumerate(fl[:2],1):
 raw=line.split('=',1)[1];raw=re.sub(r'u(\d+)',r'u**\1',raw);raw=re.sub(r'(\d)u',r'\1*u',raw);sv=ec.eval_u(-sp.expand(sp.sympify(raw,locals={'s':ss,'u':ec.U})-ss));red=[sub_scalar(q,1,sv) for q in qs]
 lam=ratio(coeff_var(red[6],3),coeff_var(red[2],3));elim=red[6]-red[2]*lam;S=cm(elim,(2,0,0,0,0,0));alpha=cm(elim,(1,0,1,0,0,0))/(2*S);beta=cm(elim,(1,0,0,0,1,0))/(2*S);gamma=cm(elim,(1,0,0,0,0,0))/(2*S);rsol=-(ec.PP.const(alpha)*ec.PP.var(2)+ec.PP.const(beta)*ec.PP.var(4)+ec.PP.const(gamma));rr=[sub_pp(q,0,rsol) for q in red]
 # q2 coefficient of u3 is bh*(h+eta*u2+delta); shift old h so new h=w.
 B=coeff_var(rr[2],5);bh=cm(B,(0,0,1,0,0,0));bu=cm(B,(0,0,0,0,1,0));bc=cm(B,(0,0,0,0,0,0));assert bh
 eta=bu/bh;delta=bc/bh
 hsol=ec.PP.var(2)-ec.PP.const(eta)*ec.PP.var(4)-ec.PP.const(delta)
 ww=[sub_pp(q,2,hsol) for q in rr]
 Bnew=coeff_var(ww[2],5);assert not (Bnew-ec.PP({(0,0,1,0,0,0):bh}))
 inds=[2,4,5,8,9,10,11]
 print('\nBR',br,'eta',eta.sing(),'delta udeg',delta.p.degree(),'surv shapes')
 for i in inds:print(i,shape(ww[i]))
 # proportional exact equations
 keep=[]
 for i in inds:
  found=False
  for j in keep:
   x=ratio(ww[i],ww[j])
   if x:print('prop',i,j,x.sing());found=True;break
  if not found:keep.append(i)
 print('keep',keep)
 (R/f'case1_branch{br}_after_w_meta.txt').write_text('s='+sv.sing()+'\nr='+rsol.sing(6)+'\nh_old='+hsol.sing(6)+'\neta='+eta.sing()+'\ndelta='+delta.sing()+'\n')
 (R/f'case1_branch{br}_after_w_eqs.txt').write_text('\n'.join(ww[i].sing(6) for i in keep))
