"""Necessary rational-univariate reduction at fivefold roots; candidate only."""
from pathlib import Path
import json,time,resource
import sympy as S
from flint import fmpz_poly
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
HERE=Path(__file__).resolve().parent;START=time.monotonic()
a,s,b=S.symbols('a s b')
cofactor=json.loads((HERE/'cofactor-candidate.json').read_text())['records'][0]
num=fmpz_poly(cofactor['a_cofactor']);den=fmpz_poly(cofactor['constant_cofactor']);common=fmpz_poly(cofactor['removed_common_factor'])
old=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-certificate.json').read_text());p=fmpz_poly(old['resultants'][1]['primitive_reduced_coefficients'])
squarefree=p//p.gcd(p.derivative());bad=squarefree.gcd(common);good=squarefree//bad
assert good.gcd(den).degree()==0
source=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_text())
def homogeneous(expression):
 poly=S.Poly(expression,a,s,domain=S.ZZ);D=poly.degree(a)
 terms=[fmpz_poly([]) for _ in range(D+1)]
 for i in range(D+1):
  q=S.Poly(poly.as_expr(),a).nth(i);q=S.Poly(q,s)
  terms[i]=fmpz_poly([int(q.nth(j)) for j in range(q.degree()+1)]) if not q.is_zero else fmpz_poly([])
 N=[fmpz_poly([1])];Dens=[fmpz_poly([1])]
 for i in range(D):N.append(N[-1]*num);Dens.append(Dens[-1]*den)
 result=fmpz_poly([])
 for i,q in enumerate(terms):result+=q*N[i]*Dens[D-i]
 return result
expression=S.sympify(source['equations'][4]['polynomial']);residual=homogeneous(expression);h=good.gcd(residual)
print(json.dumps({'stage':'finite projection','squarefree_degree':squarefree.degree(),'cofactor_exception_degree':bad.degree(),'good_degree':good.degree(),'extra_equation_degree':residual.degree(),'retained_degree':h.degree(),'seconds':time.monotonic()-START}),flush=True)
L=S.Poly(S.sympify(source['linear_equations'][0]),b);maxdegree=max(S.Poly(L.nth(i),a).degree() for i in [0,1])
P=homogeneous(L.nth(1))*den**(maxdegree-S.Poly(L.nth(1),a).degree());Q=homogeneous(L.nth(0))*den**(maxdegree-S.Poly(L.nth(0),a).degree())
bad_b=h.gcd(P);print(json.dumps({'stage':'b recovery','common_b_exception_degree':bad_b.degree(),'P_degree':P.degree(),'Q_degree':Q.degree(),'seconds':time.monotonic()-START}),flush=True)
data={'scope':'Necessary generic finite projection only; no theorem or residual squarefreeness yet.','squarefree_resultant':list(map(int,squarefree.coeffs())),'cofactor_exception':list(map(int,bad.coeffs())),'generic_initial':list(map(int,good.coeffs())),'extra_equation_substitution':list(map(int,residual.coeffs())),'generic_retained':list(map(int,h.coeffs())),'a_numerator':list(map(int,num.coeffs())),'a_denominator':list(map(int,den.coeffs())),'b_numerator':list(map(int,(-Q).coeffs())),'b_denominator':list(map(int,P.coeffs())),'b_exception':list(map(int,bad_b.coeffs()))}
(HERE/'rur-candidate.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'status':'CANDIDATE_ONLY','bytes':(HERE/'rur-candidate.json').stat().st_size,'seconds':time.monotonic()-START}),flush=True)
