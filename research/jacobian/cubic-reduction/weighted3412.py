#!/usr/bin/env python3
"""Complete312/412 frontier models and changed boundary identities."""
from pathlib import Path
from itertools import product
import sympy as s
import json,hashlib
HERE=Path(__file__).resolve().parent;OUT=HERE/'weighted-312-412';OUT.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,V):return [{'exponent':list(e),'coefficient':str(c)} for e,c in s.Poly(s.expand(p),*V).terms() if c]
def main():
 a,b,c,u,v,U,V=s.symbols('a b c u v U V');models=[]
 for p in [3,4]:
  weights=[p,-1,-2];basis=[(i,e) for i in range(3) for e in product(range(8),repeat=3) if 2<=sum(e)<=7 and sum(x*y for x,y in zip(e,weights))==weights[i]];names=list(s.symbols('c0:'+str(len(basis))));F=[a,b,c]
  for z,(i,e) in zip(names,basis):F[i]+=z*a**e[0]*b**e[1]*c**e[2]
  f,g,h=[s.expand(z.subs({a:1,b:u,c:v})) for z in F];br=lambda x,y:s.diff(x,u)*s.diff(y,v)-s.diff(x,v)*s.diff(y,u);J=s.expand(f*br(g,h)+g*br(f,h)/p+2*h*br(g,f)/p);actual=s.expand(s.det(s.Matrix(F).jacobian([a,b,c])));lift=0
  for (eu,ev),co in s.Poly(J,u,v).terms():assert (eu+2*ev)%p==0;lift+=co*a**((eu+2*ev)//p)*b**eu*c**ev
  assert s.expand(lift-actual)==0;parent=HERE/'weight-resonance7/mixed-cases'/('p'+str(p)+'q1r2.json');d={'status':'complete-primary-model','source_sha256':sha(Path(__file__)),'positive_weight':p,'weights':weights,'basis':[{'output':i,'exponent':list(e)} for i,e in basis],'parameter_names':list(map(str,names)),'source_outputs':[enc(z,[a,b,c]+names) for z in F],'cover_outputs':[enc(z,[u,v]+names) for z in [f,g,h]],'quotient_jacobian':enc(J,[u,v]+names),'native_generators':[enc(z,names) for z in s.Poly(J-1,u,v).coeffs()],'parent_timeout_model_sha256':sha(parent)};ip=OUT/('p'+str(p)+'.json');ip.write_text(json.dumps(d,indent=2)+'\n');models.append({'positive_weight':p,'input':ip.name,'input_sha256':sha(ip),'atoms':len(basis),'generators':len(d['native_generators'])});print('MODEL',p,len(basis),len(d['native_generators']),flush=True)
 A,B,C,D,E,G,L,M=s.symbols('A B C D E G L M');params=[A,B,C,D,E,G,L,M];uu=a*b**4;vv=a*b*b*c;F=[a*(1+A*vv+B*uu),b*(1+C*vv+D*uu),c*(1+E*vv+G*uu)+b*b*(L+M*uu)];Js=s.Poly(s.expand(s.det(s.Matrix(F).jacobian([a,b,c]))),a,b,c);J=0
 for (ea,eb,ec),co in Js.terms():assert ea>=ec and eb==4*(ea-ec)+2*ec;J+=co*U**(ea-ec)*V**ec
 diag=[s.expand(J).coeff(U,0).coeff(V,j) for j in [1,2,3]];assert all(s.expand(x-y)==0 for x,y in zip(diag,[2*A+3*C+2*E,4*A*C+3*A*E+4*C*E,5*A*C*E]));pure=s.expand(J.subs({A:0,C:0,E:0}));assert s.expand(pure-(1+G*U)*(1+(2*B+5*D)*U+6*B*D*U*U))==0;p4={'parameters':list(map(str,params)),'source_outputs':[enc(z,[a,b,c]+params) for z in F],'invariant_jacobian':enc(J,[U,V]+params),'diagonal_coefficients':[enc(z,params) for z in diag],'pure_U_factorization':enc(pure,[U]+params)}
 mu,P,Q,M,N,B,S,H,J,K,X=s.symbols('mu P Q M N B S H J K X');params=[mu,P,Q,M,N,B,S,H,J,K,X];uu=a*b**3;vv=a*b*c;F=[a*(1+mu*vv+P*uu+Q*vv*vv),b*(1+M*vv+N*uu+B*vv*vv),c*(1+H*vv+J*uu+X*vv*vv)+b*b*(S+K*uu)];Js=s.Poly(s.expand(s.det(s.Matrix(F).jacobian([a,b,c]))),a,b,c);JJ=0
 for (ea,eb,ec),co in Js.terms():assert ea>=ec and eb==3*(ea-ec)+ec;JJ+=co*U**(ea-ec)*V**ec
 zero=s.expand(JJ).coeff(U,0);prod=(1+mu*V+Q*V*V)*(1+M*V+B*V*V)*(1+H*V+X*V*V);assert s.expand(zero-s.diff(V*prod,V))==0;pure=s.expand(JJ.subs({mu:0,Q:0,M:0,B:0,H:0,X:0}));assert s.expand(pure-(1+J*U)*(1+(2*P+4*N)*U+5*P*N*U*U))==0;p3={'parameters':list(map(str,params)),'source_outputs':[enc(z,[a,b,c]+params) for z in F],'invariant_jacobian':enc(JJ,[U,V]+params),'U_zero_derivative_product':enc(zero,[V]+params),'pure_U_factorization':enc(pure,[U]+params)}
 sp=OUT/'structural-identities.json';sp.write_text(json.dumps({'status':'complete-primary-structural-identities','source_sha256':sha(Path(__file__)),'p4_boundary_reduced':p4,'p3_Lzero_boundary_reduced':p3},indent=2)+'\n');(OUT/'checkpoint.json').write_text(json.dumps({'status':'complete','source_sha256':sha(Path(__file__)),'models':models,'structural_identities_sha256':sha(sp),'old_timeouts_not_rerun':True},indent=2)+'\n');print('STRUCTURAL IDENTITIES COMPLETE')
if __name__=='__main__':main()
