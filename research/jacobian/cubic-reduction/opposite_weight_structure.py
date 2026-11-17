#!/usr/bin/env python3
"""Small universal invariant identities; no elimination result is a premise."""
from pathlib import Path
import sympy as sp
import json,hashlib
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,V):return [{'exponent':[int(k) for k in e],'coefficient':str(c)} for e,c in sorted(sp.Poly(p,*V).terms()) if c]
def main():
 u,v,m,q=sp.symbols('u v m q');A,B,C,D,E,G,H,I,K,L,M=sp.symbols('A B C D E G H I K L M');P=[m,q,A,B,C,D,E,G,H,I,K,L,M];f=1+A*v+B*u*v+C*v*v;g=1+D*v+E*u*v+G*v*v;h=1+H*v+I*u*v+K*v*v;k=L+M*v;R=u*h+v*k;br=lambda f,g:sp.diff(f,u)*sp.diff(g,v)-sp.diff(f,v)*sp.diff(g,u);J=sp.expand(g*(f*sp.diff(R,u)+R*sp.diff(f,u)+q*v*br(R,f))+m*v*(f*br(R,g)+R*br(f,g)));poly=sp.Poly(J-1,u,v)
 expected=[(q+2)*B+(m+1)*E+2*I,(m+q+2)*B*E+(q+3)*B*I+(m+2)*E*I,(m+q+3)*B*E*I];actual=[sp.expand(poly.coeff_monomial(u**i*v**i)) for i in [1,2,3]];assert all(sp.expand(a-b)==0 for a,b in zip(actual,expected)),[sp.expand(a-b) for a,b in zip(actual,expected)]
 masks=[]
 for mask in range(1,8):
  active=[int(bool(mask&(1<<j))) for j in range(3)];degree=sum(active);e=actual[degree-1].subs({x:0 for x,on in zip([B,E,I],active) if not on});coefficient=1+(q+1)*active[0]+m*active[1]+active[2];wanted=coefficient*sp.prod(x for x,on in zip([B,E,I],active) if on);assert sp.expand(e-wanted)==0;masks.append({'mask':mask,'degree':degree,'coefficient':enc(coefficient,P)})
 pure=sp.expand(J.subs({B:0,E:0,I:0}));ff=f.subs(B,0);gg=g.subs(E,0);hh=h.subs(I,0);factor=sp.expand(hh*(gg*(ff+q*v*sp.diff(ff,v))+m*v*ff*sp.diff(gg,v)));assert pure==factor
 leading=[]
 for df in range(3):
  for dg in range(3):
   F=[1,1+A*v,1+A*v+C*v*v][df];Gg=[1,1+D*v,1+D*v+G*v*v][dg];expr=sp.expand(Gg*(F+q*v*sp.diff(F,v))+m*v*F*sp.diff(Gg,v));lc=sp.Poly(F,v).LC()*sp.Poly(Gg,v).LC();coef=sp.expand(sp.Poly(expr,v).coeff_monomial(v**(df+dg)));assert coef==sp.expand((1+q*df+m*dg)*lc);leading.append({'degrees':[df,dg],'coefficient':enc(coef,P)})
 result={'status':'complete-primary-universal-identities','source_sha256':sha(Path(__file__)),'parameter_names':list(map(str,P)),'dictionary':enc(J,[u,v]+P),'diagonal_coefficients':[enc(p,P) for p in actual],'support_masks':masks,'pure_v_factorization':enc(pure,[v]+P),'leading_cases':leading,'scope':'Formal m,q identities; application uses positive integers and characteristiczero. Field-unit/leading-degree deductions are written.'};folder=HERE/'opposite-weight-reflection';(folder/'universal-identities.json').write_text(json.dumps(result,indent=2)+'\n');print('UNIVERSAL:3 diagonal coefficients,7 support masks,9 leading-degree cases; polynomial identities exact',flush=True)
if __name__=='__main__':main()
