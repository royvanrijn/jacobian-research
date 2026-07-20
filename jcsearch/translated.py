"""Translated two-divisor chart with divisors x=0 and xy+1=0."""
from dataclasses import dataclass
import sympy as sp
from .laurent import linear_certificate

x,y,u,v=sp.symbols("x y u v")
h=1+x*y

@dataclass(frozen=True)
class TranslatedTwoDivisorChart:
    @property
    def expressions(self): return (h/x,x**2/h)
    @property
    def inverse(self): return (u*v,u-1/(u*v))
    @property
    def jacobian(self): return sp.factor(sp.Matrix(self.expressions).jacobian((x,y)).det())
    @property
    def reciprocal_outer_jacobian(self): return -u
    @property
    def signature(self):
        return {"divisors":("x=0","xy+1=0"),"jacobian":"-1/u",
                "inverse_poles":("u=0","v=0"),"kind":"translated_two_divisor"}
    def pullback(self,f): return sp.cancel(f.subs({u:self.expressions[0],v:self.expressions[1]},simultaneous=True))

def generic_laurent(support,prefix):
    coefficients=sp.symbols(f"{prefix}0:{len(support)}")
    return sum(c*u**i*v**j for c,(i,j) in zip(coefficients,support)),coefficients

def cancellation_system(support,prefix="c"):
    """Exact common-denominator divisibility conditions for a Laurent support."""
    expression,coefficients=generic_laurent(support,prefix)
    ab=[(i-j,-i+2*j) for i,j in support] # h^alpha x^beta after pullback
    H=max(0,-min(alpha for alpha,_ in ab));X=max(0,-min(beta for _,beta in ab))
    denominator=sp.expand(x**X*h**H)
    remainders=[]
    for alpha,beta in ab:
        numerator=sp.expand(h**(alpha+H)*x**(beta+X))
        remainder=sp.reduced(numerator,[denominator],x,y)[1]
        remainders.append(sp.expand(remainder))
    combined=sp.expand(sum(c*r for c,r in zip(coefficients,remainders)))
    equations=sp.Poly(combined,x,y).coeffs() if combined else []
    certificate,matrix,rhs=linear_certificate(equations,coefficients)
    monomials=[m for m,_ in sp.Poly(combined,x,y).terms()] if combined else []
    return {"expression":expression,"coefficients":coefficients,"denominator":denominator,
      "denominator_orders":{"x":X,"xy+1":H},"remainders":remainders,"equations":equations,
      "matrix":matrix,"rhs":rhs,"certificate":certificate,
      "earliest_remainder_monomial":min(monomials,key=lambda m:(sum(m),m)) if monomials else None}

def kernel_expressions(system):
    support_expr=sp.Matrix([system["expression"]]).jacobian(system["coefficients"])
    return tuple(sp.expand((support_expr*vector)[0]) for vector in system["matrix"].nullspace())

def support_families():
    planted={(0,0),(1,1),(1,0),(-1,-1)}
    diamond2={(i,j) for i in range(-2,3) for j in range(-2,3) if abs(i)+abs(j)<=2}|planted
    box2={(i,j) for i in range(-2,3) for j in range(-2,3)}|planted
    asymmetric={(i,j) for i in range(-2,4) for j in range(-2,3)}|planted
    box3={(i,j) for i in range(-3,4) for j in range(-3,4)}|planted
    shifted42={(i,j) for i in range(-3,4) for j in range(-2,4)}|planted
    return {"diamond_r2":tuple(sorted(diamond2)),"box_r2":tuple(sorted(box2)),
            "asymmetric_30":tuple(sorted(asymmetric)),"shifted_42":tuple(sorted(shifted42)),
            "box_r3":tuple(sorted(box3))}
