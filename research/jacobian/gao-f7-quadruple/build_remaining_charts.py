"""Construct finite-free unit witnesses for P1 poles and ell-leading poles.
Untrusted generator; exported inverses require literal independent replay.
"""
from pathlib import Path
import json,time,sympy as S
from flint import nmod_poly
HERE=Path(__file__).resolve().parent;START=time.monotonic();P=1000000007

def trim(v):
 v=list(v)
 while v and not v[-1]:v.pop()
 return v

def xgcd(f,g,zero,one,add,sub,mul,inv):
 def trim(v):
  v=list(v)
  while v and v[-1]==zero:v.pop()
  return v
 def padd(f,g,minus=False):
  return trim([(sub if minus else add)(f[i] if i<len(f) else zero,g[i] if i<len(g) else zero) for i in range(max(len(f),len(g)))])
 def pmul(f,g):
  if not f or not g:return []
  out=[zero]*(len(f)+len(g)-1)
  for i,c in enumerate(f):
   for j,d in enumerate(g):out[i+j]=add(out[i+j],mul(c,d))
  return trim(out)
 def remdiv(f,g):
  f=trim(f);q=[zero]*max(0,len(f)-len(g)+1);vi=inv(g[-1])
  while f and len(f)>=len(g):
   k=len(f)-len(g);c=mul(f[-1],vi);q[k]=c
   for i,d in enumerate(g):f[k+i]=sub(f[k+i],mul(c,d))
   f=trim(f)
  return trim(q),f
 u,v=[one],[];u1,v1=[],[one]
 while g:
  q,r=remdiv(f,g);f,g=g,r;u,u1=u1,padd(u,pmul(q,u1),True);v,v1=v1,padd(v,pmul(q,v1),True)
 c=inv(f[-1]);return [[mul(c,t) for t in w] for w in [f,u,v]]

class Algebra:
 def __init__(self,h,f):
  self.h=h;self.f=f;self.n=len(f)-1;self.bz=nmod_poly([],P);self.bo=nmod_poly([1],P);self.zero=tuple([self.bz]*self.n);self.one=self.base(self.bo)
 def base(self,q):return tuple([q%self.h]+[self.bz]*(self.n-1))
 def const(self,c):return self.base(nmod_poly([c%P],P))
 def add(self,u,v):return tuple((c+d)%self.h for c,d in zip(u,v))
 def sub(self,u,v):return tuple((c-d)%self.h for c,d in zip(u,v))
 def mul(self,u,v):
  co=[self.bz]*(2*self.n-1)
  for i,c in enumerate(u):
   for j,d in enumerate(v):co[i+j]=(co[i+j]+c*d)%self.h
  for k in range(len(co)-1,self.n-1,-1):
   for j in range(self.n):co[k-self.n+j]=(co[k-self.n+j]-co[k]*self.f[j])%self.h
  return tuple(co[:self.n])
 def binv(self,u):
  g,v,w=u.xgcd(self.h);assert g.degree()==0 and int(g[0])!=0
  return v*pow(int(g[0]),-1,P)%self.h
 def inv(self,u):
  g,v,w=xgcd(list(u),self.f,self.bz,self.bo,lambda a,b:(a+b)%self.h,lambda a,b:(a-b)%self.h,lambda a,b:a*b%self.h,self.binv)
  assert len(g)==1 and g[0]==self.bo
  result=tuple((v+[self.bz]*self.n)[:self.n]);assert self.mul(u,result)==self.one;return result
 def pow(self,u,n):
  out=self.one
  while n:
   if n&1:out=self.mul(out,u)
   u=self.mul(u,u);n//=2
  return out
 def encode(self,u):return [list(map(int,c.coeffs())) for c in u]
 def evaluate(self,expr,variables,values):
  terms=S.Poly(expr,*variables).terms();powers=[[self.one] for _ in values]
  for j,value in enumerate(values):
   for k in range(max(ex[j] for ex,c in terms)):powers[j].append(self.mul(powers[j][-1],value))
  out=self.zero
  for ex,c in terms:
   term=self.const(int(c))
   for j,k in enumerate(ex):term=self.mul(term,powers[j][k])
   out=self.add(out,term)
  return out
 def t_inverse(self,coeff,modulus):
  g,u,v=xgcd(coeff,modulus,self.zero,self.one,self.add,self.sub,self.mul,self.inv);assert len(g)==1 and g[0]==self.one
  assert len(u)<len(modulus)
  return u

if __name__=='__main__':
 s,a,b=S.symbols('s a b');src=json.loads((HERE/'quadruple-input.json').read_text());old=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_text());curve=S.sympify(old['equations'][3]['polynomial']);L=S.Poly(S.sympify(old['linear_equations'][0]),b);f2=S.sympify(src['equations'][0]['polynomial']);rows=[]
 for name in ['coordinate-poles','leading-minus','leading-plus']:
  if name=='coordinate-poles':
   cp=json.loads((HERE/'coordinate-poles-candidate.json').read_text());h=nmod_poly([int(c)%P for c in cp['projection_coefficients']],P);dummy=Algebra(h,[nmod_poly([],P),nmod_poly([1],P)]);sv=dummy.base(nmod_poly([0,1],P));num=nmod_poly([int(c)%P for c in cp['a_numerator']],P)%h;den=nmod_poly([int(c)%P for c in cp['a_denominator']],P)%h;av=dummy.base(num*dummy.binv(den)%h)
   fb=S.Poly(f2,b);f=[dummy.evaluate(fb.nth(j),[s,a],[sv,av])[0] for j in range(3)];lead=dummy.binv(f[-1]);f=[v*lead%h for v in f];A=Algebra(h,f);sv=A.base(nmod_poly([0,1],P));av=A.base(av[0]);bv=tuple([A.bz,A.bo])
  else:
   h=nmod_poly([-1 if name=='leading-minus' else 1,0,0,3],P);dummy=Algebra(h,[nmod_poly([],P),nmod_poly([1],P)]);sv=dummy.base(nmod_poly([0,1],P));gp=S.Poly(curve,a);f=[dummy.evaluate(gp.nth(j),[s],[sv])[0] for j in range(9)];assert not f[-1];f=f[:8];lead=dummy.binv(f[-1]);f=[v*lead%h for v in f];A=Algebra(h,f);sv=A.base(nmod_poly([0,1],P));av=tuple([A.bz,A.bo]+[A.bz]*5);pv=A.evaluate(L.nth(1),[s,a],[sv,av]);qv=A.evaluate(L.nth(0),[s,a],[sv,av]);pinv=A.inv(pv);bv=A.sub(A.zero,A.mul(qv,pinv))
  print(json.dumps({'stage':'chart','name':name,'rank':int(h.degree())*A.n,'seconds':time.monotonic()-START}),flush=True)
  quotient=[A.evaluate(S.sympify(r['polynomial']),[s,a,b],[sv,av,bv]) for r in src['quotient_coefficients']];assert quotient[-1]==A.one
  derivative=[A.mul(A.const(j+1),quotient[j+1]) for j in range(8)];dinv=A.t_inverse(derivative,quotient)
  marked=A.zero
  for c in reversed(quotient):marked=A.add(A.mul(marked,sv),c)
  minv=A.inv(marked);defining_derivative=tuple(A.f[j+1]*(j+1) for j in range(A.n));finv=A.inv(defining_derivative)
  row={'name':name,'rank':int(h.degree())*A.n,'base_modulus':list(map(int,h.coeffs())),'extension_modulus':[list(map(int,c.coeffs())) for c in A.f],'a_mod':A.encode(av),'b_mod':A.encode(bv),'quotient_coefficients_mod':[A.encode(c) for c in quotient],'derivative_inverse_mod':[A.encode(c) for c in dinv],'marked_value_mod':A.encode(marked),'marked_inverse_mod':A.encode(minv),'extension_derivative_inverse_mod':A.encode(finv)}
  if name!='coordinate-poles':row['b_denominator_inverse_mod']=A.encode(pinv)
  rows.append(row)
 rec={'status':'UNVERIFIED_FINITE_FREE_UNIT_CANDIDATE','prime':P,'charts':rows,'boundary':'Actual characteristic-zero degrees and source identities, good leading coefficients, and literal modular inverse replay in both arithmetic engines are required.'}
 (HERE/'remaining-charts-candidate.json').write_text(json.dumps(rec,separators=(',',':'))+'\n');print(json.dumps({'stage':'saved','ranks':[r['rank'] for r in rows],'seconds':time.monotonic()-START}))
