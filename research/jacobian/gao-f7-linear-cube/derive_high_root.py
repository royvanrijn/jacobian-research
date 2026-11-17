"""Reproducible quadratic-linear elimination for the sixfold-root gate.

Generation only. The separate verifier rederives the mathematical source.
"""
from pathlib import Path
import hashlib,json,time
import sympy as S
HERE=Path(__file__).resolve().parent
START=time.monotonic()
s,a,b=S.symbols('s a b')
source=HERE.parent/'gao-f7-fibres/sixfold-input.json'
f=[S.Poly(S.sympify(x),b) for x in json.loads(source.read_text())['polynomials']]
assert len(f)==4 and all(q.degree()==2 for q in f)
A,B,C=[f[0].nth(i) for i in [2,1,0]]
assert S.expand(A-4*s*s*(3*s**3-2)**2)==0
contents=[4*s*(3*s**3-2),S.Integer(4),4*s*s]
linear=[]
for j,content in zip(range(1,4),contents):
 expression=S.cancel((A*f[j].as_expr()-f[j].nth(2)*f[0].as_expr())/content)
 assert S.denom(expression)==1
 L=S.Poly(expression,b);assert L.degree()==1
 linear.append(L)
def strip_exception_factors(expression):
 q=S.Poly(expression,a,s,domain=S.QQ).clear_denoms()[1].primitive()[1].as_expr()
 factors=[]
 for divisor in [s,3*s**3-2]:
  while True:
   quotient,remainder=S.div(q,divisor,a,s)
   if remainder!=0:break
   q=quotient;factors.append(str(divisor))
 return S.Poly(q,a,s),factors
rows=[]
for i in range(3):
 for j in range(i+1,3):
  expression=linear[i].nth(1)*linear[j].nth(0)-linear[j].nth(1)*linear[i].nth(0)
  q,factors=strip_exception_factors(S.expand(expression))
  rows.append({'type':'linear_cross','indices':[i,j],'polynomial':str(q.as_expr()),'removed_exception_factors':factors,'degrees':q.degree_list(),'terms':len(q.terms())})
for i,L in enumerate(linear):
 P,Q=L.nth(1),L.nth(0)
 q,factors=strip_exception_factors(S.expand(A*Q*Q-B*P*Q+C*P*P))
 rows.append({'type':'quadratic_resultant','indices':[i],'polynomial':str(q.as_expr()),'removed_exception_factors':factors,'degrees':q.degree_list(),'terms':len(q.terms())})
data={'source':'../gao-f7-fibres/sixfold-input.json','variables':['a','s'],'scope':'Necessary two-variable equations for an F7 root of multiplicity at least6 on s*(3s^3-2)!=0. The exceptional root values are not discarded and require a separate proof. No exclusion yet.','source_quadratic_coeff':str(A),'linear_contents':list(map(str,contents)),'linear_equations':[str(L.as_expr()) for L in linear],'equations':rows}
(HERE/'high-root-bivariate-input.json').write_text(json.dumps(data,indent=2)+'\n')
for name,V,expressions in [('high-root-szero-integral',[a,b],[q.as_expr().subs(s,0) for q in f]),('high-root-cubic-exception-integral',[a,b,s],[S.rem(S.Poly(q.as_expr(),s),S.Poly(3*s**3-2,s)).as_expr() for q in f]+[3*s**3-2])]:
 polys=[S.Poly(q,*V,domain=S.QQ).clear_denoms()[1].primitive()[1].as_expr() for q in expressions]
 case={'variables':list(map(str,V)),'equations':[{'polynomial':str(q)} for q in polys],'scope':'Necessary equations at '+('s=0' if len(V)==2 else '3s^3=2')+'; no exclusion without replay.'}
 (HERE/(name+'-input.json')).write_text(json.dumps(case,indent=2)+'\n')
print(json.dumps({'status':'NECESSARY_SYSTEMS_DERIVED','seconds':time.monotonic()-START,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'theorem_exclusion':False}))
