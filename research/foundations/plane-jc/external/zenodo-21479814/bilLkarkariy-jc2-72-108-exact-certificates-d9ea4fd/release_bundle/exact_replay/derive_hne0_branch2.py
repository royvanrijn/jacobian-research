from pathlib import Path
import pickle,re,sympy as sp
from flint import fmpq,fmpq_poly
import exact_core as ec
R=Path('.'); st=pickle.loads((R/'case1_checkpoint.pkl').read_bytes())
def dk(v):return ec.K(fmpq_poly([fmpq(n,d) for n,d in v]))
def dp(v):return ec.PP({tuple(m):dk(c) for m,c in v})
qs=[dp(q) for q in st['all_comp']]
def pw(a,n):
 o=ec.PP.const(1)
 while n:
  if n&1:o=o*a
  a=a*a;n//=2
 return o
def ss(q,i,v):
 z={}
 for m,c in q.d.items():
  M=list(m);e=M[i];M[i]=0;M=tuple(M);z[M]=z.get(M,ec.K(0))+c*v**e
 return ec.PP(z)
def spq(q,i,v):
 o=ec.PP()
 for m,c in q.d.items():
  M=list(m);e=M[i];M[i]=0;o+=ec.PP({tuple(M):c})*pw(v,e)
 return o
def cv(q,i):
 z={}
 for m,c in q.d.items():
  if m[i]==1:M=list(m);M[i]=0;z[tuple(M)]=c
 return ec.PP(z)
def rat(a,b):
 if set(a.d)!=set(b.d) or not a.d:return None
 m=next(iter(a.d));r=a.d[m]/b.d[m]
 return r if all(a.d[k]==b.d[k]*r for k in a.d) else None
def cm(q,m):return q.d.get(tuple(m),ec.K(0))
line=(R/'factor_q0_exact_clean.out').read_text().splitlines()[1];x=line.split('=',1)[1];x=re.sub(r'u(\d+)',r'u**\1',x);x=re.sub(r'(\d)u',r'\1*u',x);s=sp.Symbol('s');sv=ec.eval_u(-sp.expand(sp.sympify(x,locals={'s':s,'u':ec.U})-s));red=[ss(q,1,sv) for q in qs]
lam=rat(cv(red[6],3),cv(red[2],3));el=red[6]-red[2]*lam;S=cm(el,(2,0,0,0,0,0));a=cm(el,(1,0,1,0,0,0))/(2*S);b=cm(el,(1,0,0,0,1,0))/(2*S);g=cm(el,(1,0,0,0,0,0))/(2*S);rsol=-(ec.PP.const(a)*ec.PP.var(2)+ec.PP.const(b)*ec.PP.var(4)+ec.PP.const(g));rr=[spq(q,0,rsol) for q in red]
B=cv(rr[2],5);bh=cm(B,(0,0,1,0,0,0));eta=cm(B,(0,0,0,0,1,0))/bh;delta=cm(B,(0,0,0,0,0,0))/bh;hsol=ec.PP.var(2)-ec.PP.const(eta)*ec.PP.var(4)-ec.PP.const(delta);ww=[spq(q,2,hsol) for q in rr]
inds=[2,4,5,8,9,10,11];keep=[]
for i in inds:
 if not any(rat(ww[i],ww[j]) for j in keep):keep.append(i)
q0=ww[2];hp=ec.PP({(0,0,1,0,0,0):bh});assert not(cv(q0,5)-hp);N=(q0-hp*ec.PP.var(5))*(-bh.inv());hv=ec.PP.var(2)
def sub(q):
 d=max(m[5] for m in q.d);o=ec.PP()
 for m,c in q.d.items():
  k=m[5];M=list(m);M[5]=0;o+=ec.PP({tuple(M):c})*pw(N,k)*pw(hv,d-k)
 mh=min(m[2] for m in o.d)
 if mh:o=ec.PP({(m[0],m[1],m[2]-mh,m[3],m[4],m[5]):c for m,c in o.d.items()})
 return o
out=[]
for i in keep[1:]:
 q=sub(ww[i]);e={}
 for m,c in q.d.items():assert m[0]==m[1]==m[5]==0;e[(m[2],m[3],m[4])]=c
 out.append(e)
ser=[]
for e in out:
 ee={}
 for m,c in e.items():ee[m]={j:(int(c.p[j].p),int(c.p[j].q)) for j in range(len(c.p)) if c.p[j]}
 ser.append(ee)
(R/'hne0_branch2_deg35.pkl').write_bytes(pickle.dumps(ser,protocol=4));print([len(e) for e in ser],[max(map(sum,e)) for e in ser])
