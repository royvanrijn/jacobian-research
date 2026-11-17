"""Exact algebra checks for the standalone paper.

The proof for every order is in main.tex. The two finite moment values
below check normalization, not infinitude or the reductive classification.
Run with the repository's Python environment (SymPy is required).
"""
import json
import sympy as S

x,y,z,w,t,W=S.symbols('x y z w t W')
checks={}
for n in (8,12):
    f=x**8+14*x**4*y**4+y**8 if n==8 else x*y*(x**10+11*x**5*y**5-y**10)
    g=x**2*y**2*(x**4-y**4)**2 if n==8 else x**20+y**20-228*x**15*y**5+494*x**10*y**10+228*x**5*y**15
    v=S.Poly(g,x,y).total_degree()
    j=S.expand((S.diff(f,x)*S.diff(g,y)-S.diff(f,y)*S.diff(g,x))/(n*v))
    T=x**12-33*x**8*y**4-33*x**4*y**8+y**12 if n==8 else x**30+y**30+522*(x**25*y**5-x**5*y**25)-10005*(x**20*y**10+x**10*y**20)
    stated=x*y*(x**4-y**4)*T/6 if n==8 else -T/12
    assert S.Poly(S.expand(j-stated),x,y).is_zero
    relation=j**2-f**3*g/36+3*g**3 if n==8 else j**2-12*f**5-g**3/144
    assert S.Poly(S.expand(relation),x,y).is_zero
    checks[f'degree_{n}_Jacobian_and_Klein_identity']='PASS'

f=z**2*w**2
R=(z+f*x)**2*(w+f*y)**2/f
assert S.cancel(R-(1+x*z*w**2)**2*(1+y*z**2*w)**2)==0
checks['normalizer_polynomial_factorization']='PASS'
lhs=1/(x*(x+1)*(W-x)*(W+1-x))
rhs=(1/x+1/(W-x))/(W*(W+1))-(1/(x+1)+1/(W+1-x))/((W+1)*(W+2))
assert S.cancel(lhs-rhs)==0
checks['symbolic_partial_fraction_identity']='PASS'

moments={}
for m in (1,2):
    L=2*m-2
    beta=S.factorial(4*m)**2/S.factorial(6*m-1)/S.factorial(L)*sum(
        (-1)**r*S.binomial(L,r)/S.Integer((4*m-r)*(4*m-r-1)*(2*m+r+1)*(2*m+r+2))
        for r in range(L+1))
    closed=S.Rational(2*m,3*m+1)*S.factorial(2*m)*S.factorial(4*m)/S.factorial(6*m)
    assert beta==closed>0
    moments[str(m)]=str(beta)
checks['radial_normalizations']=moments
print(json.dumps({'exact_checks':checks,'scope':'Polynomial identities and two normalizations only; main.tex contains the all-order proofs.','formal_verification':False,'external_review':False},indent=2))
