"""General 3D weighted-lift construction and seed geometry."""
from __future__ import annotations

from dataclasses import dataclass
import sympy as sp

x,y,z,w=sp.symbols("x y z w")


@dataclass(frozen=True)
class WeightedSeedModel:
    """Exact one-variable data controlling a weighted lift.

    ``seed`` is p(w), ``primitive`` is H(w)=integral_0^w p(s) ds, and the
    inverse pencil over a target (A,B,C) is

        H(w) - B*C*w + c*A*C**2.
    """

    seed: sp.Expr
    c: sp.Expr = sp.Integer(1)
    b: sp.Expr = sp.Integer(1)

    def __post_init__(self):
        seed = sp.expand(sp.sympify(self.seed))
        c = sp.sympify(self.c)
        b = sp.sympify(self.b)
        assert c != 0 and b != 0
        assert seed.subs(w, 0) == 0
        assert seed.subs(w, 1) == -c
        primitive = sp.integrate(seed, (w, 0, w))
        assert sp.expand(primitive.subs(w, 1)) == 0
        kappa = sp.cancel(sp.diff(seed, w).subs(w, 1)/c)
        assert kappa != -2
        object.__setattr__(self, "seed", seed)
        object.__setattr__(self, "c", c)
        object.__setattr__(self, "b", b)
        object.__setattr__(self, "primitive", sp.expand(primitive))
        object.__setattr__(self, "kappa", kappa)
        object.__setattr__(self, "a", sp.cancel(-(1+kappa)/(2+kappa)))

    @property
    def seed_degree(self) -> int:
        return sp.Poly(self.seed, w).degree()

    @property
    def fiber_degree(self) -> int:
        return sp.Poly(self.primitive, w).degree()

    def inverse_polynomial(self, A, B, C):
        return sp.expand(self.primitive - B*C*w + self.c*A*C**2)

    def discriminant(self, s, t):
        return sp.factor(sp.discriminant(self.primitive - s*w + t, w))

    def branch_parameterization(self):
        """Return (s(w),t(w)) for the repeated-root discriminant."""
        return self.seed, sp.expand(w*self.seed - self.primitive)

    def zero_profile(self):
        """Return endpoint multiplicities and the remaining H(w) factor."""
        polynomial = sp.Poly(self.primitive, w)
        multiplicity_zero = 0
        while polynomial.eval(0) == 0:
            polynomial = polynomial.exquo(sp.Poly(w, w))
            multiplicity_zero += 1
        multiplicity_one = 0
        factor_one = sp.Poly(w-1, w)
        while polynomial.eval(1) == 0:
            polynomial = polynomial.exquo(factor_one)
            multiplicity_one += 1
        return multiplicity_zero, multiplicity_one, sp.factor(polynomial.as_expr())

    def mapping(self):
        vv=x*y;ss=x**2*z;uu=1+vv;gamma=1+self.a*vv+self.b*ss;W=uu*gamma
        q=sp.integrate(w*sp.diff(self.seed,w)/self.c,w)
        alpha=sp.cancel(uu+q.subs(w,W)/gamma**2)
        beta=sp.cancel(self.c+self.seed.subs(w,W)/gamma)
        F=tuple(sp.cancel(f) for f in (alpha/x**2,beta/x,x*gamma))
        assert all(sp.denom(f)==1 for f in F)
        return F

    def boundary_map(self):
        return tuple(sp.expand(component.subs(x, 0)) for component in self.mapping())


def canonical_seed(degree: int, c=sp.Integer(1)):
    """Seed p=H' for H=c*w**degree*(1-w), with degree >= 2."""
    if degree < 2:
        raise ValueError("seed degree must be at least two")
    return sp.expand(sp.diff(sp.sympify(c)*w**degree*(1-w), w))


def deformation_basis(index: int):
    """Endpoint- and integral-preserving seed deformation of degree index+2."""
    if index < 1:
        raise ValueError("deformation index must be positive")
    correction = sp.Rational(6, (index+2)*(index+3))
    return sp.expand(w*(1-w)*(w**index-correction))

def construct(seed,c=sp.Integer(1),b=sp.Integer(1)):
    model = WeightedSeedModel(seed, c, b)
    F = model.mapping()
    assert sp.factor(sp.Matrix(F).jacobian((x,y,z)).det())==model.b*model.c
    return F
