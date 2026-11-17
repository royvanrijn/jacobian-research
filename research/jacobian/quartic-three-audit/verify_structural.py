#!/usr/bin/env python3
"""Exact structural identities used in the written quartic audit, plus degree5 extension controls."""
from pathlib import Path
from math import gcd
import hashlib,json,time
import sympy as S
from verify_square_pencils import det,twice
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
x,y,z=S.symbols('x y z');X=(x,y,z)
def coefficient(p,mon):return S.Poly(S.expand(p),*X).coeff_monomial(mon)
def zcoeff(p,k):return S.Poly(S.expand(p),z).nth(k)
def main():
 start=time.monotonic();checks=[]
 for ss in [0,1]:
  degree=4-ss;cs=S.symbols('f0:'+str(3*(degree+1)));P=S.Matrix([sum(cs[i*(degree+1)+j]*x**j*y**(degree-j) for j in range(degree+1)) for i in range(3)]);rs=S.symbols('r0:15');R=S.zeros(3,1) if ss==0 else S.Matrix([sum(rs[i*5+j]*x**j*y**(4-j) for j in range(5)) for i in range(3)]);H4=z**ss*P+R;J4=H4.jacobian(X);nn=P.diff(x).cross(P.diff(y));cc=S.Matrix(S.symbols('c0:3'));bb=S.Matrix([x,y,0]);H3=z**3*cc+z*z*bb;d8=S.expand(S.trace(J4.adjugate()*H3.jacobian(X)));assert S.expand(zcoeff(d8,2*ss+2)-3*(nn.dot(cc)))==0;zero={v:0 for v in cc};assert S.expand(zcoeff(d8.subs(zero),2*ss+1)-S.Rational(8-3*ss,4-ss)*nn.dot(bb))==0;d7=twice(J4,(z*z*bb).jacobian(X));assert S.expand(zcoeff(d7,ss+3)+(8-3*ss)*P[2])==0;checks.append('binary-leading-s'+str(ss))
 gs=S.symbols('g0:5');g=z*z+gs[0]*x*z+gs[1]*y*z+gs[2]*x*x+gs[3]*x*y+gs[4]*y*y;N=S.Matrix([x*x,x*y,y*y]);nn=S.Matrix([y*y,-2*x*y,x*x]);kk=S.diff(g,z)*S.Matrix(X)-4*g*S.Matrix([0,0,1]);J4=(g*N).jacobian(X);assert (J4.adjugate()+g*kk*nn.T).applyfunc(S.expand)==S.zeros(3);us=S.symbols('u0:5');vs=S.symbols('v0:5');qm=[x*x,x*y,y*y,x*z,y*z];U=sum(c*m for c,m in zip(us,qm));V=sum(c*m for c,m in zip(vs,qm));H3=S.Matrix([2*x*U,y*U+x*V,2*y*V]);cc=S.Matrix(S.symbols('c0:3'));H2=z*z*cc;assert S.expand(zcoeff(S.trace(J4.adjugate()*H2.jacobian(X))+twice(J4,H3.jacobian(X)),5)-4*nn.dot(cc))==0;L=S.Matrix(S.symbols('l0:3'));assert S.expand(zcoeff(S.trace(J4.adjugate()*(z*L).jacobian(X)),4)-2*nn.dot(L))==0;checks.append('common-quadratic-factor-adjugate-and-top-z')
 # Cubic scalar critical-point eliminants and exact coordinate inverse.
 c,d,e,f,g0,h,i,j,b,a=S.symbols('c d e f g h i j b a');AA=d*c-e*x;BB=-x**3-h*x+d*g0;DD=-e*e/(2*d);EE=c*x-e*(x*x+h)/(2*d)+i;elim=S.expand(3*BB**2-DD*AA*BB+EE*AA**2);assert S.Poly(elim,x).LC()==3
 cusp=x*x*y+z**3+c*x*z+e*y*z+g0*x+h*y+i*z+j;zz=-(x*x+h)/e;yy=-(3*zz*zz+c*x+i)/e;v=S.expand(S.diff(cusp,x).subs({z:zz,y:yy},simultaneous=True));assert S.Poly(v,x).LC()==-6/e**3
 AA=4*b*(x+c)-a*a;BB=-a*(x*x+e)+2*b*f;DD=-a*x;EE=-x**3-e*x+b*d;elim=S.expand(b*BB**2-DD*AA*BB+EE*AA**2);assert S.Poly(elim,x).LC()==-16*b*b
 tangent=x*x*y+x*z*z+a*y*z+c*z*z+d*x+e*y+f*z+j;zz=-(x*x+e)/a;yy=-(2*(x+c)*zz+f)/a;v=S.expand(S.diff(tangent,x).subs({z:zz,y:yy},simultaneous=True));assert S.Poly(v,x).LC()==5/a**2
 ss=x*y+z*z+d;NN=S.Matrix([x,y-2*z*ss/f-x*ss*ss/f**2,z+x*ss/f]);II=S.Matrix([x,y+2*z*ss/f-x*ss*ss/f**2,z-x*ss/f]);assert S.expand(ss.subs(dict(zip(X,NN)),simultaneous=True)-ss)==0;assert NN.subs(dict(zip(X,II)),simultaneous=True).applyfunc(lambda p:S.expand(p))==S.Matrix(X);assert II.subs(dict(zip(X,NN)),simultaneous=True).applyfunc(lambda p:S.expand(p))==S.Matrix(X);assert S.expand(NN.jacobian(X).det())==1;assert max(S.Poly(v,*X).total_degree() for v in II)==5;checks.append('all-four-cubic-critical-eliminant-leading-coefficients-and-Nagata-inverse')
 # Exceptional scalar and the two new algebraic cancellation identities.
 rr=x+x*x*y;der=lambda v:S.expand(-x*x*S.diff(v,x)+(1+2*x*y)*S.diff(v,y));assert der(rr)==0 and der(1/x)==1;assert S.gcd(S.diff(rr,x),S.diff(rr,y))==1
 lam=S.symbols('lambda');qq=S.Function('q')(x,y);bb=S.Function('b')(x,y);hh=S.Function('h')(rr);cc=S.Function('c')(x,y);PP=lam*qq*z*z+(2*lam*bb+hh)*z+cc;QQ=qq*z+bb;formula=qq*der(cc)-(2*lam*bb+hh)*der(bb);bracket=S.expand(der(PP)*S.diff(QQ,z)-S.diff(PP,z)*der(QQ));assert S.expand(bracket-formula+der(qq)*(lam*qq*z*z+(2*lam*bb+hh)*z))==0
 # Identity after the invariant-leading conditions Dq=Dh=0.
 rhs=der(qq*cc-lam*bb*bb-hh*bb);assert S.expand(rhs-formula-der(qq)*cc)==0;checks.append('exceptional-rational-slice-and-proportional-invariant-leading-cancellation')
 cases=[]
 for m in range(1,5):
  for n in range(1,m+1):
   if m==n:kind='affine-pair' if m==1 else 'equal-degree-constant-shear-or-constant-leading'
   elif m>=3:
    assert m//gcd(m,n)>5-m;assert n==1 or 5-n<=3;kind='ratio-fibre-factor-or-UFD-degree-contradiction'
   else:assert (m,n)==(2,1);kind='quadratic-affine-fibre-factor-or-affine-target-reduction'
   cases.append(dict(m=m,n=n,p_degree_at_most=5-m,q_degree_at_most=5-n,closed_by=kind))
 assert len(cases)==10
 out=dict(status='pass',source_sha256=sha(Path(__file__)),paper_html_sha256=sha(HERE/'source-v1.html'),helper_sha256=sha(HERE/'verify_square_pencils.py'),sympy_version=S.__version__,identities=checks,degree5_exceptional_cases=cases,all_positive_z_degree_pairs_through_total5_covered=True,scope='Exact algebraic controls for the written quartic source audit and the separate degree5 exceptional-scalar extension. Geometric residue periods, UFD/field arguments and universal case proof remain written, not formalized. No full arbitrary coefficient search or global simplicity optimum claimed.',seconds=round(time.monotonic()-start,3));(HERE/'structural-validation.json').write_text(json.dumps(out,indent=2)+'\n');print('STRUCTURAL IDENTITIES PASS',len(checks),'groups;degree5cases',len(cases),'seconds',out['seconds'],flush=True)
if __name__=='__main__':main()
