#!/usr/bin/env python3
"""Complete weighted ansatze and Euler/top-Jacobian identities."""
from pathlib import Path
from math import gcd
import sympy as sp
import json,hashlib,time
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,V):return [{'exponent':[int(k) for k in e],'coefficient':str(c)} for e,c in sorted(sp.Poly(p,*V).terms()) if c]
def main():
 a,b,c=sp.symbols('a b c');out=[];start=time.monotonic();pairs=[(p,q) for total in range(3,7) for p in range(1,total) for q in [total-p] if gcd(p,q)==1]
 for pw,qw in pairs:
  layers=2 if pw+qw==3 else 1;params=[];blocks={};blockparams={}
  for name,degree in [('Q',pw),('T',2*pw),('R',pw+1),('S',pw+1),('U',2*pw+1),('V',2*pw+1)]:
   if layers==1 and name in ['T','U','V']:blocks[name]=sp.Integer(0);blockparams[name]=[];continue
   vs=sp.symbols(name+'0:'+str(degree+1));params.extend(vs);blockparams[name]=list(vs);blocks[name]=sum(v*b**(degree-i)*c**i for i,v in enumerate(vs))
  Q,T,R,S,U,V=[blocks[k] for k in ['Q','T','R','S','U','V']];X=[a,b,c]+params;F=[a+a**(qw+1)*Q+a**(2*qw+1)*T,b+a**qw*R+a**(2*qw)*U,c+a**qw*S+a**(2*qw)*V];J=sp.expand(sp.Matrix(F).jacobian([a,b,c]).det());jac=lambda f,g:sp.expand(sp.diff(f,b)*sp.diff(g,c)-sp.diff(f,c)*sp.diff(g,b));identities=[]
  def check(name,power,zero,expected):
   sub={v:0 for k in zero for v in blockparams[k]};actual=sp.expand(sp.Poly(J.subs(sub),a).coeff_monomial(a**power));assert sp.expand(actual-expected)==0,(pw,qw,name);identities.append({'name':name,'a_power':power,'zero_blocks':zero,'coefficient':enc(actual,X)})
  check('first',qw,[],sp.diff(R,b)+sp.diff(S,c)+(qw+1)*Q)
  if layers==1:
   check('one-layer-top',3*qw,[],sp.Rational(pw+qw+1,pw+1)*Q*jac(R,S))
  else:
   check('one-layer-boundary-top',3*qw,['T','U','V'],sp.Rational(pw+qw+1,pw+1)*Q*jac(R,S))
   check('full-top',6*qw,[],sp.Rational(2*pw+2*qw+1,2*pw+1)*T*jac(U,V))
   check('quartic-first-full-top',5*qw,['T'],sp.Rational(2*pw+qw+1,2*pw+1)*Q*jac(U,V))
   check('no-high-tail-top',4*qw,['U','V'],sp.Rational(pw+2*qw+1,pw+1)*T*jac(R,S))
   check('one-high-tail-quartic-first',4*qw,['T','V'],((pw+qw+1)*Q*jac(U,S)-qw*U*jac(Q,S))/(pw+1))
   check('one-high-tail-mixed',5*qw,['V'],((pw+2*qw+1)*T*jac(U,S)-qw*U*jac(T,S))/(pw+1))
  record={'weights':[pw,-qw,-qw],'layers':layers,'parameter_names':list(map(str,params)),'blocks':{k:enc(v,X) for k,v in blocks.items()},'block_parameter_indices':{k:[params.index(v) for v in vs] for k,vs in blockparams.items()},'outputs':[enc(f,X) for f in F],'identities':identities,'degree':1+layers*(pw+qw),'nonlinear_atoms':len(params)};out.append(record);print('STRUCTURAL',pw,qw,'layers',layers,'atoms',len(params),'identities',len(identities),flush=True)
 result={'status':'complete-primary-structural-identities','source_sha256':sha(Path(__file__)),'pairs':out,'seconds':round(time.monotonic()-start,3),'scope':'All primitive positive pairs with3<=p+q<=6; larger sums are affine at degree7. Ratio1/1 is an explicit prior premise, not silently covered.'};folder=HERE/'equal-negative-smallweights-v2';folder.mkdir(exist_ok=True);(folder/'structural-models.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
