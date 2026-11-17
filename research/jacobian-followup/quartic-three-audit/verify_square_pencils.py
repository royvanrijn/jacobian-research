#!/usr/bin/env python3
"""Independent exact reconstruction of Wang-Yan3.4/3.9 coefficient exclusions.
No author code or numerical parameter samples; monomial-ideal reductions only.
"""
from pathlib import Path
from itertools import permutations
import sympy as S
import hashlib,json,time
HERE=Path(__file__).resolve().parent
x,y,z=S.symbols('x y z');X=(x,y,z)
a,b,al,be,ga,de,ep,rho,A,B,C,D,E,G=S.symbols('a b alpha beta gamma delta epsilon rho A B C D E G')
M=S.symbols('M1:4');N=S.symbols('N1:4');L=S.symbols('L1:4');ALL=X+(a,b,al,be,ga,de,ep,rho,A,B,C,D,E,G)+M+N+L
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def det(T):
 return S.expand(sum((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))*S.prod(T[i,p[i]] for i in range(3)) for p in permutations(range(3))))
def mixed(U,V,W):
 return S.expand(sum(det(S.Matrix.hstack(U[:,p[0]],V[:,p[1]],W[:,p[2]]))*(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3)) for p in permutations(range(3))))
def twice(U,V):
 return S.expand(sum(det(S.Matrix.hstack(*[U[:,i] if i==j else V[:,i] for i in range(3)])) for j in range(3)))
def red(p,gens):
 f=S.Poly(S.expand(p),*ALL);out=0
 for es,c in f.terms():
  if any(all(es[ALL.index(v)]>=k for v,k in g.items()) for g in gens):continue
  out+=c*S.prod(v**e for v,e in zip(ALL,es) if e)
 return S.expand(out)
def coeff(p,e):return S.Poly(S.expand(p),*X).coeff_monomial(e)
def mons(d):return [x**i*y**j*z**(d-i-j) for i in range(d,-1,-1) for j in range(d-i,-1,-1)]
def vector(p,d):
 rows=[]
 for f in p:
  poly=S.Poly(S.expand(f),*X);rows.extend(poly.coeff_monomial(m) for m in mons(d))
 return S.Matrix(rows)
def operator_matrix(q,degree):
 nn=S.Matrix([q*q,-2*x*x*q,x**4]);op=(lambda f:z*S.diff(f,z)-y*S.diff(f,y)) if q==y*z else (lambda f:-x*S.diff(f,y)+2*y*S.diff(f,z));cols=[]
 for i in range(3):
  for m in mons(degree):
   candidate=4*x*nn[i]*op(m);j4=S.Matrix([x**4,x*x*q,q*q]).jacobian(X);direct=sum(j4.adjugate()[j,i]*S.diff(m,X[j]) for j in range(3));assert S.expand(candidate-direct)==0;cols.append(vector([candidate],degree+5))
 return S.Matrix.hstack(*cols),op,nn

def run(tag):
 q=y*z if tag=='yz' else y*y+x*z;H4=S.Matrix([x**4,x*x*q,q*q]);J4=H4.jacobian(X);adj=J4.adjugate();ell=a*y+b*z;I0=[{a:1,al:1},{b:1,al:1}]
 # Complete unnormalized cubic kernel:30 slots,10-dimensional family.
 uu=S.symbols('u1:4');vv=S.symbols('v1:4');u=sum(t*s for t,s in zip(uu,X));v=sum(t*s for t,s in zip(vv,X));fam=S.Matrix([2*x*x*u+al*x*q,q*u+x*x*v+be*x**3+ga*x*q,2*q*v+de*x**3]);params=uu+vv+(al,be,ga,de);cubic_op,op,nn=operator_matrix(q,3);assert cubic_op.rank()==20;basis=S.Matrix.hstack(*(vector(fam.diff(p),3) for p in params));assert basis.rank()==10 and cubic_op*basis==S.zeros(cubic_op.rows,10)
 # Quadratic kernel is exactly the six constant combinations of x^2,q.
 quadratic_op,_,_=operator_matrix(q,2);assert quadratic_op.rank()==12;kb=[]
 for i in range(3):
  for m in [x*x,q]:v0=[0]*3;v0[i]=m;kb.append(vector(v0,2))
 kernel=S.Matrix.hstack(*kb);assert kernel.rank()==6 and quadratic_op*kernel==S.zeros(quadratic_op.rows,6)
 if tag=='yz':
  H3=S.Matrix([2*x*x*ell+al*x*q,be*x**3+ga*x*q+q*ell,de*x**3+ep*x*q]);H2=S.Matrix([ell**2-ga*x*ell+A*x*x+B*q,(6*be-ep)*x*ell/4+C*x*x+D*q,3*de*x*ell/2+E*x*x+G*q]);U=S.Matrix([A+ga**2/2,C-ga*(6*be+ep)/8,E-(3*de*ga+ep**2)/4])
 else:
  H3=S.Matrix([2*x*x*ell+al*x*q,q*ell+rho*x*x*z+be*x**3+ga*x*q,2*rho*q*z+de*x**3]);H2=S.Matrix([ell**2+A*x*x+B*q-ga*x*ell+al*rho*x*z,rho*z*ell+C*x*x+D*q+3*be*x*ell/2+ga*rho*x*z,rho**2*z*z+E*x*x+G*q+3*de*x*ell/2]);U=S.Matrix([A+ga**2/2,C-3*be*ga/4,E-3*de*ga/4])
 J3=H3.jacobian(X);J2=H2.jacobian(X);H1=S.Matrix([x*M[i]+y*N[i]+z*L[i] for i in range(3)]);J1=H1.jacobian(X);rhs=vector([-twice(J4,J3)],7);assert quadratic_op.rows==rhs.rows
 residuals=[S.expand((v.T*rhs)[0]) for v in quadratic_op.T.nullspace()];res=[]
 for v in residuals:
  if v:assert red(v,I0)==0;res.append([S.Poly(v,a,b,al,be,ga,de,ep,rho).coeff_monomial(a*al),S.Poly(v,a,b,al,be,ga,de,ep,rho).coeff_monomial(b*al)])
 assert all(S.expand(v-c[0]*a*al-c[1]*b*al)==0 for v,c in zip([v for v in residuals if v],res));assert S.Matrix(res).rank()==2;assert S.expand(S.trace(adj*J3))==0;d7=S.expand(S.trace(adj*J2)+twice(J4,J3));assert red(d7,I0)==0;print('PASS',tag,'complete cubic/quadratic kernels and exact compatibility',flush=True)
 d6=S.expand(S.trace(adj*J1)+mixed(J4,J3,J2)+det(J3));checks=0
 if tag=='yz':
  for mon,j,mult,which in [(x*y**3*z**2,0,4,'a'),(x*y*y*z**3,0,-4,'b'),(x**3*y*y*z,1,-8,'a'),(x**3*y*z*z,1,8,'b'),(x**5*y,2,4,'a'),(x**5*z,2,-4,'b')]:
   expected=mult*((a*U[j]-N[j]) if which=='a' else (b*U[j]-L[j]));assert red(coeff(d6,mon)-expected,I0)==0;checks+=1
  assert det(S.Matrix.hstack(S.Matrix(M),a*U,b*U))==0
 else:
  assert red(coeff(d6,x*x*y**4)-coeff(d6,x**4*z*z)-18*a*b*be,I0)==0;assert red(coeff(d6,x**5*z)-coeff(d6,x**4*y*y)-9*a*b*de,I0)==0;checks+=2;I=I0+[{a:1,b:1,be:1},{a:1,b:1,de:1}]
  for mon,j,mult in [(x*x*y**4,0,4),(x**4*y*y,1,-8),(x**6,2,4)]:assert red(coeff(d6,mon)-mult*(a*U[j]-N[j]),I)==0;checks+=1
  sub={N[i]:a*U[i] for i in range(3)};d6sub=S.expand(d6.subs(sub,simultaneous=True));d5=S.expand(twice(J4,J2)+mixed(J4,J3,J1)+twice(J2,J3));d5=S.expand(d5.subs(sub,simultaneous=True));relations=[(coeff(d5,x*y**4)-a*coeff(d6sub,x**3*y**3),0,2),(coeff(d5,x**3*y*y)-2*a*coeff(d6sub,x**5*y),1,-4),(coeff(d5,x**5),2,2)]
  for value,j,mult in relations:assert red(value-mult*a*(M[j]-ga*U[j]),I)==0;checks+=1
  assert S.expand(det(S.Matrix.hstack(S.Matrix(M),a*U,S.Matrix(L)))-det(S.Matrix.hstack(a*S.Matrix(M),U,S.Matrix(L))))==0;assert det(S.Matrix.hstack(a*ga*U,U,S.Matrix(L)))==0
 print('PASS',tag,'all selected coefficient identities',checks,flush=True)
 return dict(pencil=tag,cubic_operator_rank=20,cubic_kernel_dimension=10,quadratic_operator_rank=12,quadratic_kernel_dimension=6,compatibility_ideal=['a*alpha','b*alpha'],compatibility_left_residuals=[str(v) for v in residuals if v],compatibility_relation_span_rank=2,complete_quadratic_family_checked=True,selected_linear_identities_checked=checks,no_parameter_division=True,linear_jacobian_forced_singular=True)
def main():
 start=time.monotonic();assert sha(HERE/'source-v1.html')=='04dfa22c70f97783543eb1a011d90e4d4b9f4fb8f7af394ef6969c1aa491ca29';results=[]
 for tag in ['yz','parabolic']:
  results.append(run(tag));(HERE/'square-pencils-checkpoint.json').write_text(json.dumps(dict(completed=results,source_sha256=sha(Path(__file__))),indent=2)+'\n')
 out=dict(status='pass',source_sha256=sha(Path(__file__)),paper_html_sha256=sha(HERE/'source-v1.html'),pencils=results,scope='Exact coefficient/operator replays for the two nonbinary square conic pencils only. Mathematical leading normal forms and other sections are separately written/audited; not a complete quartic theorem proof.',seconds=round(time.monotonic()-start,3));(HERE/'square-pencils-validation.json').write_text(json.dumps(out,indent=2)+'\n');print('ALL SQUARE PENCIL IDENTITIES PASS',out['seconds'],flush=True)
if __name__=='__main__':main()
