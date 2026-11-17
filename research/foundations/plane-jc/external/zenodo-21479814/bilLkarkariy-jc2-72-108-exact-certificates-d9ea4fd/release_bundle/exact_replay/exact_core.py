#!/usr/bin/env python3
from __future__ import annotations
import re
from dataclasses import dataclass
from pathlib import Path
import sympy as sp
from flint import fmpq,fmpq_poly
ROOT=Path(__file__).resolve().parent
t=sp.symbols('t');a={i:sp.Symbol(f'a{i}') for i in range(1,9)};d={i:sp.Symbol(f'd{i}') for i in range(2,13)}
A0=sum(a[i]*t**i for i in a).subs({a[1]:1,a[8]:1});D0=sum(d[i]*t**i for i in d)
J4=sp.Poly(sp.expand(2*A0*sp.diff(D0,t)-3*sp.diff(A0,t)*D0-t**2),t);dsubs={};res=[]
for k in range(J4.degree()+1):
 e=sp.cancel(J4.nth(k).subs(dsubs))
 if e==0:continue
 solved=False
 for v in d.values():
  if e.has(v) and sp.Poly(e,v).degree()==1 and not any(e.has(w) for w in d.values() if w!=v):
   pv=sp.Poly(e,v);c=pv.nth(1);dsubs[v]=sp.cancel(-(e-c*v)/c);solved=True;break
 if not solved:res.append(e)
assert (len(dsubs),len(res))==(11,6)
loc={str(v):v for v in a.values()};rels=[]
for line in (ROOT/'firstblock_Q_exact.out').read_text().splitlines():
 m=re.match(r'L\[(\d+)\]=(.*)',line.strip())
 if m:rels.append(sp.sympify(m.group(2).replace('^','**'),locals=loc))
assert len(rels)==6
U=a[7];Hprim=sp.Poly(rels[0],U,domain=sp.QQ).primitive()[1].as_expr();Hexpr=sp.Poly(Hprim,U,domain=sp.QQ).monic().as_expr()
LEX={a[1]:sp.Integer(1),a[8]:sp.Integer(1),a[7]:U}
for q in rels[1:]:
 for v in [a[i] for i in range(2,7)]:
  qq=sp.expand(q);co=qq.coeff(v)
  if co!=0 and not (qq-co*v).has(v):LEX[v]=sp.cancel(-(qq-co*v)/co);break
assert len(LEX)==8

def qpoly(expr):
 p=sp.Poly(sp.cancel(expr),U,domain=sp.QQ); coeff=[fmpq(0)]*(p.degree()+1)
 for (i,),c in p.terms():c=sp.Rational(c);coeff[i]=fmpq(int(c.p),int(c.q))
 return fmpq_poly(coeff)
MOD=qpoly(Hexpr);assert MOD.degree()==35
@dataclass(frozen=True)
class K:
 p:fmpq_poly
 def __init__(self,x=0):
  if isinstance(x,K):p=x.p
  elif isinstance(x,fmpq_poly):p=x
  elif isinstance(x,fmpq):p=fmpq_poly([x])
  elif isinstance(x,sp.Rational):p=fmpq_poly([fmpq(int(x.p),int(x.q))])
  elif isinstance(x,int):p=fmpq_poly([x])
  else:p=qpoly(x)
  object.__setattr__(self,'p',p%MOD)
 def __add__(self,o):return K(self.p+K(o).p)
 __radd__=__add__
 def __neg__(self):return K(-self.p)
 def __sub__(self,o):return self+(-K(o))
 def __rsub__(self,o):return K(o)-self
 def __mul__(self,o):return K((self.p*K(o).p)%MOD)
 __rmul__=__mul__
 def inv(self):
  if not self:raise ZeroDivisionError
  gg,s,_=self.p.xgcd(MOD)
  if gg.degree()!=0:raise ZeroDivisionError
  return K(s/gg[0])
 def __truediv__(self,o):return self*K(o).inv()
 def __pow__(self,n):
  if n<0:return self.inv()**(-n)
  out=K(1);x=self
  while n:
   if n&1:out=out*x
   x=x*x;n//=2
  return out
 def __bool__(self):return not self.p.is_zero()
 def __eq__(self,o):return self.p==K(o).p
 def sing(self):
  terms=[]
  for i in range(len(self.p)):
   c=self.p[i]
   if not c:continue
   num=int(c.p);den=int(c.q);cs=str(num) if den==1 else f'({num}/{den})';mon='' if i==0 else ('u' if i==1 else f'u^{i}')
   if not mon:term=cs
   elif c==1:term=mon
   elif c==-1:term='-'+mon
   else:term=f'{cs}*{mon}'
   terms.append(term)
  return '0' if not terms else '+'.join(terms).replace('+-','-')
NP=6;PNAMES=['r','s','h','u1','u2','u3']
class PP:
 def __init__(self,data=None):self.d={m:K(c) for m,c in (data or {}).items() if K(c)}
 @staticmethod
 def const(c):return PP({(0,)*NP:K(c)}) if K(c) else PP()
 @staticmethod
 def var(i):
  m=[0]*NP;m[i]=1;return PP({tuple(m):K(1)})
 def __add__(self,o):
  o=toPP(o);z=dict(self.d)
  for m,c in o.d.items():
   q=z.get(m,K(0))+c
   if q:z[m]=q
   elif m in z:del z[m]
  return PP(z)
 __radd__=__add__
 def __neg__(self):return PP({m:-c for m,c in self.d.items()})
 def __sub__(self,o):return self+(-toPP(o))
 def __rsub__(self,o):return toPP(o)-self
 def __mul__(self,o):
  if isinstance(o,(int,K,fmpq)):return PP({m:c*o for m,c in self.d.items()})
  o=toPP(o);z={}
  for m,c in self.d.items():
   for n,e in o.d.items():
    k=tuple(i+j for i,j in zip(m,n));z[k]=z.get(k,K(0))+c*e
  return PP(z)
 __rmul__=__mul__
 def __bool__(self):return bool(self.d)
 def sing(self,nparams=NP):
  terms=[]
  for m,c in sorted(self.d.items(),key=lambda x:(sum(x[0]),x[0]),reverse=True):
   mons=[]
   for i,e in enumerate(m[:nparams]):
    if e:mons.append(PNAMES[i] if e==1 else f'{PNAMES[i]}^{e}')
   cs=c.sing();mon='*'.join(mons)
   if not mon:term=cs
   elif cs=='1':term=mon
   elif cs=='-1':term='-'+mon
   else:term=f'({cs})*{mon}'
   terms.append(term)
  return '0' if not terms else '+'.join(terms).replace('+-','-')
def toPP(x):return x if isinstance(x,PP) else PP.const(x)
def tadd(a,b):
 n=max(len(a),len(b));return [(a[i] if i<len(a) else PP())+(b[i] if i<len(b) else PP()) for i in range(n)]
def tscale(a,c):return [x*c for x in a]
def tmul(a,b):
 z=[PP() for _ in range(len(a)+len(b)-1)]
 for i,x in enumerate(a):
  for j,y in enumerate(b):z[i+j]=z[i+j]+x*y
 return z
def tder(a):return [a[i]*i for i in range(1,len(a))]
def rref(rows,rhs):
 rows=[r[:] for r in rows];rhs=rhs[:];m=len(rows);n=len(rows[0]);piv=[];rr=0
 for col in range(n):
  k=next((i for i in range(rr,m) if rows[i][col]),None)
  if k is None:continue
  rows[rr],rows[k]=rows[k],rows[rr];rhs[rr],rhs[k]=rhs[k],rhs[rr]
  inv=rows[rr][col].inv();rows[rr]=[x*inv for x in rows[rr]];rhs[rr]=rhs[rr]*inv
  for i in range(m):
   if i!=rr and rows[i][col]:
    q=rows[i][col];rows[i]=[x-q*y for x,y in zip(rows[i],rows[rr])];rhs[i]=rhs[i]-rhs[rr]*q
  piv.append(col);rr+=1
 free=[j for j in range(n) if j not in piv];compat=[rhs[i] for i in range(rr,m) if rhs[i]]
 return rows,rhs,piv,free,compat
def linear_solve(rows,rhs,param_start):
 R,Y,piv,free,compat=rref(rows,rhs);sol={}
 for j,fc in enumerate(free):sol[fc]=PP.var(param_start+j)
 for i,pc in enumerate(piv):
  q=Y[i]
  for fc in free:
   if R[i][fc]:q=q-sol[fc]*R[i][fc]
  sol[pc]=q
 return sol,free,compat,param_start+len(free)
def eval_u(expr):return K(qpoly(expr))
Avec=[K(0)]*9;Dvec=[K(0)]*13;Avec[1]=K(1);Avec[8]=K(1)
for i in range(2,8):Avec[i]=eval_u(LEX[a[i]])
for i in range(2,13):Dvec[i]=eval_u(sp.cancel(dsubs[d[i]].subs(LEX)))
for deg in range(21):
 val=K(0)
 for i,Ai in enumerate(Avec):
  for j,Dj in enumerate(Dvec):
   if i+j-1==deg:val+=Ai*Dj*(2*j-3*i)
 if deg==2:val-=1
 assert not val

def high_layers():
 rows=[[K(0) for _ in range(19)] for _ in range(20)]
 for j in range(1,9):
  for k,Dk in enumerate(Dvec):
   deg=j+k-1
   if deg<20:rows[deg][j-1]+=Dk*(k-3*j)
 for j in range(2,13):
  for i,Ai in enumerate(Avec):
   deg=i+j-1
   if deg<20:rows[deg][8+j-2]+=Ai*2*(j-i)
 rows=[r for r in rows if any(r)];sol1,f1,c1,np=linear_solve(rows,[PP() for _ in rows],0);assert not c1
 B=[PP() for _ in range(9)];E=[PP() for _ in range(13)]
 for j in range(1,9):B[j]=sol1[j-1]
 for j in range(2,13):E[j]=sol1[8+j-2]
 rows2=[[K(0) for _ in range(19)] for _ in range(20)]
 for j in range(1,8):
  for k,Dk in enumerate(Dvec):
   deg=j+k-1
   if deg<20:rows2[deg][j-1]+=Dk*(-3*j)
 for j in range(1,13):
  for i,Ai in enumerate(Avec):
   deg=i+j-1
   if deg<20:rows2[deg][7+j-1]+=Ai*(2*j-i)
 known=tadd(tmul(B,tder(E)),tscale(tmul(tder(B),E),K(-2)));keep=[];rhs=[]
 for deg,row in enumerate(rows2):
  const=PP();k=deg-7
  if 0<=k<len(Dvec):const+=PP.const(Dvec[k]*(-24))
  kn=known[deg] if deg<len(known) else PP()
  if any(row) or kn or const:keep.append(row);rhs.append(-(kn+const))
 sol2,f2,c2,np=linear_solve(keep,rhs,np)
 C=[PP() for _ in range(9)];F=[PP() for _ in range(13)];C[8]=PP.const(1)
 for j in range(1,8):C[j]=sol2[j-1]
 for j in range(1,13):F[j]=sol2[7+j-1]
 assert not c2 and (len(f1),len(f2),np)==(2,1,3)
 return B,E,C,F,np
