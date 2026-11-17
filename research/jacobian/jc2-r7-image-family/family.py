"""Exact coefficient charts for all normalized degree-(6,10) gap-21 pairs."""
from pathlib import Path
import sys
import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'foundations/scripts'))
from verify_f2_degree_6_10_gap import construct_obstructions

a,b,c,d,B,D,E = S.symbols('a b c d B D E')
t,x,X,Y,u = S.symbols('t x X Y u')
VARIABLES = (a,c,d,B,D,E)
_, rho, _, _, powers = construct_obstructions()
p = t**6+a*t**4+b*t**3+c*t**2+d*t
q = S.expand(sum(co*sum(powers[j][n].as_expr()*t**(2*j-n)
                        for n in range(2*j+1))
                 for j,co in [(5,1),(4,B),(2,D),(1,E)]))

g0 = 18*B*d+45*d**2-6*B*a-30*d*a+7*a**2-6*c
g1 = (108*B**2*a+405*d**2*a-108*B*a**2-378*d*a**2+a**3
      +216*B*c-216*c*d+252*c*a+324*D-72)
g2 = (810*d**3-810*d**2*a-36*B*a**2+189*d*a**2+7*a**3
      +108*B*c+243*c*d-126*c*a+81*D-18)
targets = (g0,g1,g2)

def primitive(expr, variables=VARIABLES):
    return S.Poly(S.expand(expr), *variables).clear_denoms()[1].primitive()[1].as_expr()

generators = [primitive(rho[k].as_expr().subs(b,1)) for k in [1,3,5,7,9]]
c_value = S.expand(S.solve(g0,c)[0])
D_value = S.expand(S.solve(g2,D)[0].subs(c,c_value))
E_value = S.expand(S.solve(rho[1].as_expr().subs(b,1),E)[0].subs({c:c_value,D:D_value}))
graph = {c:c_value,D:D_value,E:E_value}
reduced = [primitive(f.subs(graph),(a,d,B)) for f in
           [rho[k].as_expr().subs(b,1) for k in [3,7,9]]]

Q = 36*B**2-108*B*a+288*B*d+95*a**2-510*a*d+675*d**2
Z = 5*X-3*Y
N = 5*X**2-2*X*Y+3*Y**2
R = 5*X**4*Y**2+2*Z**3
ds = N/(6*Z)
conic = {d:ds,a:2*ds+(X-Y)/6,B:(X+7*(2*ds+(X-Y)/6)-21*ds)/6}

def complete(seed):
    result = dict(seed)
    result.update({v:S.factor(f.subs(seed)) for v,f in graph.items()})
    result[b] = S.Integer(1)
    return result

line = complete({a:2*d,d:d,B:-S.Rational(7,6)*d})
laurent = complete({a:-(15*x**6+8*x**3+4)/(4*x**2),
                    d:-(15*x**6+4*x**3+4)/(12*x**2),
                    B:-(5*x**3+1)/(3*x**2)})
