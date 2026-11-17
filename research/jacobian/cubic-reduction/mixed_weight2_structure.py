#!/usr/bin/env python3
"""Primary structural identities, separately exported for rational replay."""
from pathlib import Path
import sympy as s
import json,hashlib
HERE=Path(__file__).resolve().parent;OUT=HERE/'mixed-weight2'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,V):return [{'exponent':list(e),'coefficient':str(c)} for e,c in s.Poly(s.expand(p),*V).terms() if c]
def main():
 a,b,c,U,V=s.symbols('a b c U V');params=s.symbols('A B C D E G H I K L M');A,B,C,D,E,G,H,I,K,L,M=params;f=1+A*U+B*U**2+C*U*V;g=1+D*U+E*U**2+G*U*V;h=1+H*U+I*U**2+K*U*V;k=L+M*U;lift={U:a*b*b,V:a*a*c};F=[a*f.subs(lift),b*g.subs(lift),c*h.subs(lift)+b**4*k.subs(U,a*b*b)];J=s.Poly(s.expand(s.det(s.Matrix(F).jacobian([a,b,c]))),a,b,c);j=0
 for (ea,eb,ec),coef in J.terms():assert eb%2==0 and ea==eb//2+2*ec;j+=coef*U**(eb//2)*V**ec
 diag=[s.expand(j).coeff(U,i).coeff(V,i) for i in [1,2,3]];wanted=[4*C+3*G+2*K,6*C*G+5*C*K+4*G*K,7*C*G*K];assert all(s.expand(x-y)==0 for x,y in zip(diag,wanted));pure=s.expand(j.subs({C:0,G:0,K:0}));ff,gg,hh=[x.subs({C:0,G:0,K:0}) for x in [f,g,h]];rhs=hh*(gg*(ff+U*s.diff(ff,U))+2*U*ff*s.diff(gg,U));assert s.expand(pure-rhs)==0
 even={'parameters':list(map(str,params)),'forms':[enc(x,[U,V]+list(params)) for x in [f,g,h,k]],'source_outputs':[enc(x,[a,b,c]+list(params)) for x in F],'invariant_jacobian':enc(j,[U,V]+list(params)),'diagonal_coefficients':[enc(x,list(params)) for x in diag],'pure_U_factorization':enc(pure,[U]+list(params))}
 r,mu,P,A0,A1,B3,B5,N,H2,H4,M5,M7=s.symbols('r mu P A0 A1 B3 B5 N H2 H4 M5 M7');params=[r,mu,P,A0,A1,B3,B5,N,H2,H4,M5,M7];u,v=s.symbols('u v');phi=1+mu*u*u+P*u**4;A=A0+A1*u*u;B=u+B3*u**3+B5*u**5;C=N*u;D=1+H2*u*u+H4*u**4;E=M5*u**5+M7*u**7;g=A*v+B;h=C*v*v+D*v+E;bracket=lambda p,q:s.diff(p,u)*s.diff(q,v)-s.diff(p,v)*s.diff(q,u);J=s.expand(phi*bracket(g,h)+g*bracket(phi,h)/2+r*h*bracket(g,phi)/2);wanted=[s.diff(phi,u)*(B*D-r*A*E)/2+phi*(s.diff(B,u)*D-A*s.diff(E,u)),s.diff(phi,u)*((1-r)*A*D+2*B*C)/2+phi*(s.diff(A,u)*D+2*s.diff(B,u)*C-A*s.diff(D,u)),(2-r)*s.diff(phi,u)*A*C/2+phi*(2*s.diff(A,u)*C-A*s.diff(C,u))];assert all(s.expand(s.Poly(J,v).nth(i)-z)==0 for i,z in enumerate(wanted));pure=s.expand(J.subs({mu:0,P:0}));top=s.expand(s.Poly(pure,v).nth(2));assert top==s.expand(N*(3*A1*u*u-A0));odd={'parameters':list(map(str,params)),'forms':[enc(x,[u,v]+params) for x in [phi,A,B,C,D,E,g,h]],'univariate_dictionary':enc(J,[u,v]+params),'v_coefficients':[enc(x,[u]+params) for x in wanted],'constant_first_top':enc(top,[u]+params)}
 out={'status':'complete-primary-structural-identities','source_sha256':sha(Path(__file__)),'even_r4':even,'odd_univariate':odd,'written_all_weight_cutoffs':{'even':'r>=12 has only invariant U terms; all r>=8 classify as identity','odd':'r>=15 has only invariant U terms; r9,11,13 have only a constant c cross coefficient in output2','combined':'r>=14 is identity after Keller conditions'}};(OUT/'structural-identities.json').write_text(json.dumps(out,indent=2)+'\n');print('PRIMARY structural identities complete')
if __name__=='__main__':main()
