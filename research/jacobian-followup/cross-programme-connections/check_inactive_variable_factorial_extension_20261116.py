from pathlib import Path
import json
from math import comb, factorial
import sympy as s
x,y,z=s.symbols("x y z")
a=s.I
b=2+s.I
phase=x*x*y+a*x*x+b*x*z
def moment(poly, shapes):
    out=0
    for mon, coeff in s.Poly(s.expand(poly), x,y,z).terms():
        out+=coeff*s.prod(s.rf(alpha,k) for alpha,k in zip(shapes,mon))
    return s.expand(out)
def separate_sum(m,q,shapes):
    ax,ay,az=shapes
    out=0
    for k in range(m+1):
        for ell in range(m-k+1):
            j=m-k-ell
            out+=s.Integer(factorial(m)//(factorial(j)*factorial(k)*factorial(ell)))*a**k*b**ell*s.rf(ax,2*m-ell)*s.rf(ay,j)*s.rf(az,ell+q)
    return s.expand(out)
count=0
for shapes in [(s.Integer(1),)*3,(s.Rational(1,2),s.Rational(2,3),s.Rational(5,4))]:
    for m in range(7):
        for q in range(3):
            assert s.simplify(moment(z**q*phase**m,shapes)-separate_sum(m,q,shapes))==0
            count+=1
for m in range(7):
    assert moment((z-1)*(x*x*y+a*x*x)**m,(1,1,1))==0
    count+=1
assert s.simplify(moment((z-1)*phase,(1,1,1))-b)==0
count+=1
import mpmath as mp
mp.mp.dps=60
ratios=[]
for m in [10,30,80]:
    value=mp.mpc(0)
    for ell in range(m+1):
        ratio=mp.factorial(2*m-ell)/mp.factorial(2*m)
        for k in range(m-ell+1):
            value+=mp.j**k/mp.factorial(k)*(2+mp.j)**ell*ratio
    ratios.append({"m":m,"normalized_real":str(mp.re(value)),"normalized_imag":str(mp.im(value)),"distance_to_exp_i":str(abs(value-mp.exp(mp.j)))})
receipt={"scope":"finite coefficient identities and illustrative numerical consistency only","exact_checks":count,"all_exact_checks_pass":True,"numerical_consistency":ratios,"universal_limit_proof":"INACTIVE_VARIABLE_FACTORIAL_EXTENSION_2026_11_16.tex","independent_review":False,"registered_theorem":False}
target=Path(__file__).with_suffix(".json")
target.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(receipt,indent=2))
