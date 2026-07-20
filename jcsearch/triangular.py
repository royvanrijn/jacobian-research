"""The one-exceptional-divisor charts Phi_p and their Laurent kernels."""
from __future__ import annotations
from dataclasses import dataclass
import sympy as sp
from .laurent import laurent_dict, linear_certificate

x,y,u,v=sp.symbols("x y u v")

@dataclass(frozen=True)
class TriangularChart:
    p: sp.Expr

    @property
    def expressions(self): return (y+self.p.subs(v,1/x),1/x)
    @property
    def jacobian(self): return sp.factor(sp.Matrix(self.expressions).jacobian((x,y)).det())
    @property
    def inverse(self): return (1/v,u-self.p)
    @property
    def signature(self):
        # p is removable by a triangular target automorphism. The geometric
        # one-divisor signature is therefore independent of its coefficients.
        return {"valuation_xy":(-1,0),"jacobian_divisor":{"v=0":2},
                "pole_semigroup_generators":((-1,0),(0,1)),
                "exponent_lattice_snf":(1,1),"exceptional_divisors":1}

    def pullback(self,f): return sp.cancel(f.subs({u:self.expressions[0],v:self.expressions[1]},simultaneous=True))
    def to_chart(self,f): return sp.cancel(f.subs({x:1/v,y:u-self.p},simultaneous=True))

def rectangular_support(imax=6,jmin=-8,jmax=8):
    return tuple((i,j) for i in range(imax+1) for j in range(jmin,jmax+1))

def sparse_weighted_support(imax, jmin, jmax, weight=1):
    """A valuation staircase, preferred over a dense rectangle."""
    return tuple((i,j) for i in range(imax+1) for j in range(jmin,jmax+1)
                 if j+weight*i <= jmax)

def generic_laurent(support,prefix="c"):
    coefficients=sp.symbols(f"{prefix}0:{len(support)}")
    return sum(c*u**i*v**j for c,(i,j) in zip(coefficients,support)),coefficients

def polynomialization_system(chart:TriangularChart,support,prefix="c"):
    f,coefficients=generic_laurent(support,prefix)
    # Work in (y,v): v=1/x, so positive powers of v are forbidden x-poles.
    expanded=sp.expand(f.subs(u,y+chart.p))
    equations=[]; earliest=None
    for (ey,ev),coefficient in laurent_dict(expanded,(y,v)).items():
        if ev>0:
            equations.append(coefficient)
            key=(ev,ey)
            earliest=key if earliest is None or key<earliest else earliest
    cert,matrix,rhs=linear_certificate(equations,coefficients)
    return {"expression":f,"coefficients":coefficients,"equations":equations,
            "certificate":cert,"matrix":matrix,"rhs":rhs,"earliest_forbidden":earliest}

def kernel_basis(system):
    return system["matrix"].nullspace()

def p_catalog(max_degree=4):
    """Canonical sparse strata; leading coefficient normalized to one."""
    out=[sp.Integer(0)]
    for d in range(1,max_degree+1):
        out.extend([v**d, sum(v**k for k in range(1,d+1))])
    unique=[]
    for p in out:
        if p not in unique: unique.append(p)
    return tuple(unique)

