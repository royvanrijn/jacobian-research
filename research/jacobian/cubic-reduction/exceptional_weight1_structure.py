#!/usr/bin/env python3
"""Primary exact identities for the exceptional positive-unit weighted family."""
from pathlib import Path
import sympy as s
import json,hashlib
HERE=Path(__file__).resolve().parent;OUT=HERE/'exceptional-weight1';OUT.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,vs):return [{"exponent":list(e),"coefficient":str(c)} for e,c in s.Poly(s.expand(p),*vs).terms() if c]
def main():
 u,w,v=s.symbols('u w v');r,mu,p,z,a,b,b2,b3,b4,C,d1,d2,d3=s.symbols('r mu p z a b b2 b3 b4 C d1 d2 d3');es=s.symbols('e1:8');names=[r,mu,p,z,a,b,b2,b3,b4,C,d1,d2,d3]+list(es);phi=1+mu*u+p*u**2+z*u**3;A=a+b*u;B=u+b2*u**2+b3*u**3+b4*u**4;D=1+d1*u+d2*u**2+d3*u**3;E=sum(es[i]*u**(i+1) for i in range(7));g=A*w+B;h=C*w*w+D*w+E
 bracket=lambda x,y:s.diff(x,u)*s.diff(y,w)-s.diff(x,w)*s.diff(y,u)
 J=s.expand(g*bracket(phi,h)+r*h*bracket(g,phi)+phi*bracket(g,h));wanted=[s.diff(phi,u)*(B*D-r*A*E)+phi*(s.diff(B,u)*D-A*s.diff(E,u)),s.diff(phi,u)*((1-r)*A*D+2*B*C)+phi*(s.diff(A,u)*D+2*s.diff(B,u)*C-A*s.diff(D,u)),(2-r)*s.diff(phi,u)*A*C+2*phi*s.diff(A,u)*C]
 assert all(s.expand(s.Poly(J,w).nth(i)-x)==0 for i,x in enumerate(wanted))
 univariate={'parameters':list(map(str,names)),'forms':[enc(x,[u,w]+names) for x in [phi,g,h]],'dictionary':enc(J,[u,w]+names),'w_coefficients':[enc(x,[u]+names) for x in wanted]}
 r,A,B,alpha,beta,T=s.symbols('r A B alpha beta T');h0s=s.symbols('h0:7');h1s=s.symbols('k0:4');params=[r,A,B,alpha,beta,T]+list(h0s)+list(h1s);G0=A*u+B;G1=alpha*u+beta;H0=sum(h0s[i]*u**i for i in range(7));H1=sum(h1s[i]*u**i for i in range(4));G=G0+v*G1;H=H0+v*H1+T*v*v;graph=s.expand((G+v*s.diff(G,v))*s.diff(H,u)-s.diff(G,u)*(r*H+v*s.diff(H,v)));cofs=[G0*s.diff(H0,u)-r*A*H0,2*G1*s.diff(H0,u)+G0*s.diff(H1,u)-r*alpha*H0-(r+1)*A*H1,2*G1*s.diff(H1,u)-(r+1)*alpha*H1-(r+2)*A*T,-(r+2)*alpha*T];assert all(s.expand(s.Poly(graph,v).nth(i)-x)==0 for i,x in enumerate(cofs))
 graphs={'parameters':list(map(str,params)),'forms':[enc(x,[u,v]+params) for x in [G0,G1,H0,H1,G,H]],'dictionary':enc(graph,[u,v]+params),'v_coefficients':[enc(x,[u]+params) for x in cofs],'boundary_degree_bounds':{'4':[4,6,1,4],'5':[4,6,1,5]},'affine_vectorfield':{'q_image':'A*q','v_image':'-A*v','constant_equation':'A*(q*H_q-v*H_v-r*H)=1','source_q_v':[0,1]}}
 alpha,P,R,T=s.symbols('alpha P R T');params=[alpha,P,R,T];A=-(1+alpha)/2;mu=(3-alpha)/2;N=1+mu*u+P*u*u+R*u**3;G0=A*u-1;G1=1+alpha*u;H1=-G1**3;H0=s.series(N*G1**3,u,0,5).removeO()+T*u**5;E0=s.expand(G0*s.diff(H0,u)-5*A*H0-1);E1=s.expand(2*G1*s.diff(H0,u)+G0*s.diff(H1,u)-5*alpha*H0-6*A*H1);sol={}
 for var,k in [(P,1),(R,2),(T,4)]:sol[var]=s.factor(s.solve(s.Poly(E0,u).nth(k).subs(sol),var)[0])
 Q=s.expand(s.Poly(E1,u).nth(1).subs(sol));S=s.expand(s.Poly(E1,u).nth(2).subs(sol));L,M,gcd=s.gcdex(Q,S,alpha);assert gcd==1 and s.expand(L*Q+M*S)==1
 r5={'parameters':list(map(str,params)),'forms':[enc(x,[u]+params) for x in [N,G0,G1,H0,H1]],'E0':enc(E0,[u]+params),'E1':enc(E1,[u]+params),'triangular_elimination_indices':[1,2,4],'substitutions':{str(k):enc(z,[alpha]) for k,z in sol.items()},'residual_indices':[1,2],'residuals':[enc(x,[alpha]) for x in [Q,S]],'unit_multipliers':[enc(x,[alpha]) for x in [L,M]],'multiplier_monomials':sum(len(enc(x,[alpha])) for x in [L,M]),'primary_groebner_discovery':'Nine native reduced coefficient equations give primary basis[1]; not a proof premise. The accepted five-term unit is independently multiplied after complete triangular elimination.'}
 t,p,q=s.symbols('t p q');eq4=[p+t,2*t*p];eq5=[p+t,q+2*t*p,t*(2*q+p*p)];assert s.expand(2*t*t-(2*t*eq4[0]-eq4[1]))==0;assert s.expand(5*t**3-(eq5[2]-2*t*eq5[1]-t*(p-5*t)*eq5[0]))==0
 catalan={'parameters':['t','p','q'],'r4_equations':[enc(x,[t,p,q]) for x in eq4],'r5_equations':[enc(x,[t,p,q]) for x in eq5],'r4_target':enc(2*t*t,[t,p,q]),'r5_target':enc(5*t**3,[t,p,q])}
 out={'status':'complete-primary-identities','source_sha256':sha(Path(__file__)),'univariate':univariate,'graph':graphs,'r5':r5,'constant_first_quadratic_exclusions':catalan}
 (OUT/'structural-identities.json').write_text(json.dumps(out,indent=2)+'\n');print('PRIMARY structural identities complete; r5 unit terms',r5['multiplier_monomials'])
if __name__=='__main__':main()
