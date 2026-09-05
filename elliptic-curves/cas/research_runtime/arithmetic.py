"""Immutable, labelled arithmetic inputs. Expensive discovery lives in adapters."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .store import digest


def rationals(values):
    return tuple(str(Fraction(str(value))) for value in values)


@dataclass(frozen=True)
class CurveModel:
    coefficients: tuple[str, ...]
    base: str = "QQ"

    def __post_init__(self):
        if self.base != "QQ" or len(self.coefficients) != 5:
            raise ValueError("CurveModel requires five rational Weierstrass coefficients")
        object.__setattr__(self, "coefficients", rationals(self.coefficients))
        if not self.discriminant:
            raise ValueError("singular Weierstrass model")

    @property
    def b_invariants(self):
        a1, a2, a3, a4, a6 = map(Fraction, self.coefficients)
        return (a1*a1 + 4*a2, a1*a3 + 2*a4, a3*a3 + 4*a6,
                a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4)

    @property
    def discriminant(self):
        b2, b4, b6, b8 = self.b_invariants
        return -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6

    @property
    def two_division_polynomial(self):
        # Generator z=4*x, so the monic polynomial is integral on an
        # integral Weierstrass model and disc(f)=256*disc(E).
        b2, b4, b6, _ = self.b_invariants
        return rationals((16*b6, 8*b4, b2, 1))

    @property
    def key(self):
        return digest(self)

    def contains(self, point):
        if point is None:
            return True
        x, y = map(Fraction, point)
        a1, a2, a3, a4, a6 = map(Fraction, self.coefficients)
        return y*y + a1*x*y + a3*y == x**3 + a2*x*x + a4*x + a6


@dataclass(frozen=True)
class TwoTorsionContext:
    """Identity of a *labelled* etale cubic, shared across 2-congruent curves.

    Root permutations are not silently identified by a reduced polynomial.
    A varying curve supplies its generator map separately in ArithmeticContext.
    All optional arithmetic is accessed through this identity, never a curve ID.
    """

    polynomial: tuple[str, ...]
    labels: tuple[str, ...] = ("theta",)

    def __post_init__(self):
        object.__setattr__(self, "polynomial", rationals(self.polynomial))
        object.__setattr__(self, "labels", tuple(self.labels))
        if len(self.polynomial) != 4 or self.polynomial[-1] != "1" or not self.labels:
            raise ValueError("a labelled monic cubic is required")
        c, b, a, _ = map(Fraction, self.polynomial)
        if a*a*b*b - 4*b**3 - 4*a**3*c - 27*c*c + 18*a*b*c == 0:
            raise ValueError("the cubic algebra is not etale")

    @property
    def key(self):
        return digest(self)


@dataclass(frozen=True)
class ArithmeticContext:
    model: CurveModel
    minimal_model: CurveModel | None
    # Sage/PARI convention: input = [u^2*x+r, u^3*y+s*u^2*x+t].
    minimal_to_input: tuple[str, ...] | None
    discriminant_factorization: tuple[tuple[int, int], ...] | None
    two_torsion: TwoTorsionContext
    curve_generator_in_algebra: tuple[str, ...]

    def __post_init__(self):
        object.__setattr__(self,"curve_generator_in_algebra",rationals(self.curve_generator_in_algebra))
        if self.minimal_model is None:
            if self.minimal_to_input is not None or self.discriminant_factorization is not None:
                raise ValueError("unprepared arithmetic must not claim minimal data or a factorization")
            factors = None
        else:
            object.__setattr__(self, "minimal_to_input", rationals(self.minimal_to_input))
            factors = tuple(tuple(row) for row in self.discriminant_factorization)
            object.__setattr__(self, "discriminant_factorization", factors)
            if len(self.minimal_to_input) != 4 or Fraction(self.minimal_to_input[0]) == 0:
                raise ValueError("missing invertible model transport")
            u,r,s,t = map(Fraction,self.minimal_to_input)
            a1,a2,a3,a4,a6 = map(Fraction,self.model.coefficients)
            transformed = ((a1+2*s)/u, (a2-s*a1+3*r-s*s)/u**2,
                (a3+r*a1+2*t)/u**3,
                (a4-s*a3+2*r*a2-(t+r*s)*a1+3*r*r-2*s*t)/u**4,
                (a6+r*a4+r*r*a2+r**3-t*a3-r*t*a1-t*t)/u**6)
            if rationals(transformed) != self.minimal_model.coefficients:
                raise ValueError("model isomorphism does not preserve the exact equation")
            if any(Fraction(c).denominator != 1 for c in self.minimal_model.coefficients):
                raise ValueError("global minimal model must be integral")
        if len(self.curve_generator_in_algebra) != 3:
            raise ValueError("curve generator needs three labelled power coordinates")
        c,b,a,_ = map(Fraction,self.two_torsion.polynomial)
        def multiply(left,right):
            raw=[Fraction(0)]*5
            for i,x in enumerate(left):
                for j,y in enumerate(right):raw[i+j]+=x*y
            for i in (4,3):
                for j,coefficient in enumerate((c,b,a)):raw[i-3+j]-=raw[i]*coefficient
            return tuple(raw[:3])
        generator=tuple(map(Fraction,self.curve_generator_in_algebra))
        powers=[(Fraction(1),Fraction(0),Fraction(0)),generator]
        powers.extend([multiply(generator,generator),multiply(multiply(generator,generator),generator)])
        polynomial=tuple(map(Fraction,self.kummer_model.two_division_polynomial))
        if any(sum(polynomial[i]*powers[i][j] for i in range(4)) for j in range(3)):
            raise ValueError("curve generator does not satisfy its 2-division polynomial")
        determinant=(powers[1][1]*powers[2][2]-powers[1][2]*powers[2][1])
        if determinant == 0:
            raise ValueError("curve generator does not generate the labelled cubic algebra")
        if factors is None:
            return
        if len({p for p, _ in factors}) != len(factors):
            raise ValueError("duplicate discriminant factors")
        product = 1
        for p, e in factors:
            if type(p) is not int or type(e) is not int or p < 2 or e <= 0:
                raise ValueError("invalid discriminant factorization")
            product *= p**e
        if product != abs(self.minimal_model.discriminant):
            raise ValueError("incomplete or wrong discriminant factorization")

    @classmethod
    def for_search(cls, model):
        """Exact raw model identity; no minimalization, factorization or field setup."""
        if not isinstance(model,CurveModel):model=CurveModel(tuple(model))
        return cls(model,None,None,None,TwoTorsionContext(model.two_division_polynomial),(0,1,0))

    @property
    def prepared(self):
        return self.minimal_model is not None

    def require_prepared(self):
        if not self.prepared:
            raise ValueError("descent requires an explicitly prepared ArithmeticContext")
        return self

    @property
    def kummer_model(self):
        return self.minimal_model or self.model

    @property
    def kummer_to_input(self):
        return self.minimal_to_input or ("1","0","0","0")

    @property
    def bad_primes(self):
        self.require_prepared()
        return tuple(p for p, _ in self.discriminant_factorization)

    @property
    def key(self):
        return digest(self)

    def record(self):
        return {"model": self.model.coefficients,
                "minimal_model": self.minimal_model.coefficients if self.minimal_model else None,
                "minimal_to_input": self.minimal_to_input,
                "discriminant_factorization": self.discriminant_factorization,
                "two_torsion": {"polynomial": self.two_torsion.polynomial,
                                "labels": self.two_torsion.labels},
                "curve_generator_in_algebra": self.curve_generator_in_algebra}

    @classmethod
    def from_record(cls, row):
        return cls(CurveModel(tuple(row["model"])), CurveModel(tuple(row["minimal_model"])) if row["minimal_model"] is not None else None,
                   tuple(row["minimal_to_input"]) if row["minimal_to_input"] is not None else None,
                   tuple(map(tuple,row["discriminant_factorization"])) if row["discriminant_factorization"] is not None else None,
                   TwoTorsionContext(tuple(row["two_torsion"]["polynomial"]), tuple(row["two_torsion"]["labels"])),
                   tuple(row["curve_generator_in_algebra"]))
