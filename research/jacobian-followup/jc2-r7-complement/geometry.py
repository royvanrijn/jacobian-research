"""Exact infinity jets and closest degree-(3,5) target equations."""
from pathlib import Path
import sys
import sympy as S
from sympy.polys.rings import ring

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'jc2-r7-image-family'))
import family as F
sys.path.insert(0,str(HERE.parent/'jc2-r7-laurent-collisions'))
import laurent as K

U,V=S.symbols('U V')
MONOMIALS=tuple((i,j) for j in range(4) for i in range(6) if 3*i+5*j<=15)
PIVOTS=((5,0),(3,1),(1,2),(4,0),(2,1),(0,2),(3,0),(1,1),(2,0),(0,1),(1,0),(0,0))

def approximate_quintic(p,q,t):
    residual=S.Poly(S.expand(q**3-p**5),t)
    polynomial=V**3-U**5
    steps=[]
    while residual and residual.degree()>=10:
        degree=residual.degree()
        matches=[(i,j) for i,j in MONOMIALS if 6*i+10*j==degree]
        assert len(matches)==1
        i,j=matches[0];coefficient=residual.LC()
        polynomial-=coefficient*U**i*V**j
        residual=S.Poly(residual.as_expr()-coefficient*p**i*q**j,t)
        steps.append((degree,i,j,coefficient))
    # The remaining two lower-pole freedoms are separated from the pencil.
    c6,c0=residual.nth(6),residual.nth(0)
    polynomial-=c6*U+c0
    residual=S.Poly(residual.as_expr()-c6*p-c0,t)
    return S.expand(polynomial),S.expand(residual.as_expr()),steps

def infinity_power(alpha,order):
    """Coefficient of tau^order in (1+sum e_2j*tau^2j+A*tau^21)^alpha."""
    R,*g=ring('e2,e4,e6,e8,e10,e12,e14,e16,e18,e20,A',S.QQ)
    terms={2*(i+1):g[i] for i in range(10)};terms[21]=g[-1]
    values=[R.one]+[R.zero]*order
    for n in range(1,order+1):
        values[n]=sum((((alpha+1)*i-n)*coefficient*values[n-i]
                       for i,coefficient in terms.items() if i<=n),R.zero)/n
    return R,values[order]
