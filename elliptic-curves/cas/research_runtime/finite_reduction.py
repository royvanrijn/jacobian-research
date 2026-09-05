"""Stable finite quotients, per-point signatures and incremental column bases."""

from dataclasses import dataclass
from fractions import Fraction

from .arithmetic import CurveModel, rationals
from .binary import BinaryBasis
from .store import FiniteFieldFacts, default_store, digest


def short_presentation(model, point=None):
    """Exact completion of the square/cube; no minimization or factorization."""
    a1, a2, a3, a4, a6 = map(Fraction, model.coefficients)
    shift2, shift4, shift6 = a2+a1*a1/4, a4+a1*a3/2, a6+a3*a3/4
    coefficients = (0, 0, 0, shift4-shift2**2/3,
                    shift6-shift2*shift4/3+2*shift2**3/27)
    if point is None:
        return coefficients
    x, y = map(Fraction, point)
    return (x+shift2/3, y+(a1*x+a3)/2)


class ReductionCache:
    def __init__(self, store=None):
        self.store=store or default_store()
        self.facts=FiniteFieldFacts(self.store)
        self._quotients={}
        self._points={}
        self.quotient_builds=0
        self.point_evaluations=0

    def quotient(self, coefficients, prime):
        from mod2_reduction_independence import (_is_prime,_reduce_rational,finite_curve_points,
            finite_multiply,finite_add,finite_subtract)
        model=CurveModel(tuple(coefficients))
        if prime<=2 or not _is_prime(prime):
            raise ValueError("certificate primes must be odd")
        a,b=(_reduce_rational(Fraction(c),prime) for c in short_presentation(model)[3:])
        if (-16*(4*a**3+27*b**2))%prime==0:
            raise ValueError(f"bad reduction at {prime}")
        identity=(model.key,prime)
        if identity in self._quotients:return self._quotients[identity]

        def build():
            self.quotient_builds+=1
            points=finite_curve_points(a,b,prime)
            doubled={finite_multiply(point,2,a,prime) for point in points}
            span=[None]
            for point in points:
                if any(finite_subtract(point,rep,a,prime) in doubled for rep in span):continue
                span.extend(finite_add(rep,point,a,prime) for rep in tuple(span))
            table=[]
            for point in points:
                mask=next(mask for mask,rep in enumerate(span) if finite_subtract(point,rep,a,prime) in doubled)
                table.append([point,mask])
            return {"group_order":len(points),"doubled_order":len(doubled),
                    "dimension":len(span).bit_length()-1,"span":span,"table":table,
                    "frobenius_trace":prime+1-len(points)}

        row=self.facts.query(model.coefficients,prime,"mod2-quotient",build=build,version="exhaustive-labelled-1")
        # Replay retained quotient witnesses once per process. This checks the
        # whole finite group and all cosets, without finding a new quotient basis.
        points=finite_curve_points(a,b,prime)
        doubled={finite_multiply(point,2,a,prime) for point in points}
        span=tuple(None if point is None else tuple(point) for point in row["span"])
        table={None if point is None else tuple(point):mask for point,mask in row["table"]}
        dimension=row["dimension"]
        if (type(dimension) is not int or not 0<=dimension<=2 or len(span)!=(1<<dimension)
            or row["group_order"]!=len(points) or row["doubled_order"]!=len(doubled)
            or len(span)*len(doubled)!=len(points) or len(table)!=len(row["table"])
            or set(table)!=set(points) or span[0] is not None
            or row["frobenius_trace"]!=prime+1-len(points)):
            raise ArithmeticError("invalid retained finite quotient witness")
        for point,mask in table.items():
            if type(mask) is not int or not 0<=mask<len(span) or finite_subtract(point,span[mask],a,prime) not in doubled:
                raise ArithmeticError("invalid quotient coset coordinate")
        for i,left in enumerate(span):
            for j,right in enumerate(span):
                if table[finite_add(left,right,a,prime)] != i^j:
                    raise ArithmeticError("quotient labels are not binary coordinates")
        result=(row,table,digest(row))
        self._quotients[identity]=result
        return result

    def point_signature(self, coefficients, point, prime):
        from mod2_reduction_independence import _reduce_rational
        model=CurveModel(tuple(coefficients));point=rationals(point)
        identity=(model.key,point,prime)
        if identity in self._points:return self._points[identity]
        if not model.contains(point):raise ValueError("point does not lie on the exact curve")
        row,table,quotient_hash=self.quotient(model.coefficients,prime)
        affine = short_presentation(model,point)
        # On a good integral Weierstrass model, every nonintegral affine
        # point reduces to the identity. Its denominator does not make the
        # prime unusable for all the other cached columns.
        reduced = (None if any(c.denominator % prime == 0 for c in affine) else
                   tuple(_reduce_rational(c,prime) for c in affine))
        expected=table[reduced]
        cached=self.facts.query(model.coefficients,prime,"point-mod2",labels=(point,quotient_hash),
            build=lambda:{"mask":expected},version="exhaustive-labelled-1")
        if cached["mask"]!=expected:raise ArithmeticError("tampered per-point reduction signature")
        result=(expected,row["dimension"])
        self._points[identity]=result
        self.point_evaluations+=1
        return result

    def signature(self, coefficients, points, prime):
        from mod2_reduction_independence import Mod2ReductionSignature
        row,_,_=self.quotient(coefficients,prime)
        values=[self.point_signature(coefficients,point,prime)[0] for point in points]
        return Mod2ReductionSignature(prime,row["group_order"],row["doubled_order"],row["dimension"],
            tuple(tuple((mask>>i)&1 for mask in values) for i in range(row["dimension"])))


@dataclass(frozen=True)
class IncrementalReductions:
    model: CurveModel
    primes: tuple[int, ...]
    dimensions: tuple[int, ...]
    points: tuple[tuple[str, ...], ...]
    columns: tuple[int, ...]
    basis: BinaryBasis

    def __post_init__(self):
        for name in ("primes", "dimensions", "columns"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        object.__setattr__(self, "points", tuple(rationals(p) for p in self.points))
        if (len(self.primes) != len(self.dimensions) or len(set(self.primes)) != len(self.primes)
            or len(self.points) != len(self.columns) or any(d not in (0, 1, 2) for d in self.dimensions)):
            raise ValueError("inconsistent finite-reduction state")
        expected = BinaryBasis(sum(self.dimensions))
        for point, column in zip(self.points, self.columns):
            if not self.model.contains(point):
                raise ValueError("finite-reduction state contains an off-curve point")
            expected, _ = expected.append(column)
        if expected != self.basis:
            raise ValueError("binary column basis is inconsistent with retained columns")

    @classmethod
    def empty(cls, coefficients, primes, cache):
        model=CurveModel(tuple(coefficients));primes=tuple(primes)
        if len(set(primes))!=len(primes):raise ValueError("duplicate certificate primes")
        dimensions=tuple(cache.quotient(model.coefficients,p)[0]["dimension"] for p in primes)
        return cls(model,primes,dimensions,(),(),BinaryBasis(sum(dimensions)))

    @property
    def independent_images(self):
        return self.basis.rank==len(self.points)

    def append(self, point, cache):
        """Evaluate only the new point; an ambiguous column stays ambiguous."""
        point=rationals(point);column=0;offset=0
        if not self.model.contains(point):raise ValueError("point is off the exact model")
        for p,dimension in zip(self.primes,self.dimensions):
            mask,actual=cache.point_signature(self.model.coefficients,point,p)
            if actual!=dimension:raise ArithmeticError("stable local quotient changed")
            column|=mask<<offset;offset+=dimension
        basis,relation=self.basis.append(column)
        return IncrementalReductions(self.model,self.primes,self.dimensions,
            (*self.points,point),(*self.columns,column),basis),relation

    def escalate(self, prime, cache):
        """Add one ambiguous-case prime; old signatures are never reconstructed."""
        if prime in self.primes:raise ValueError("prime already in stable reduction set")
        dimension=cache.quotient(self.model.coefficients,prime)[0]["dimension"]
        columns=tuple(column|(cache.point_signature(self.model.coefficients,point,prime)[0]<<self.basis.width)
                      for point,column in zip(self.points,self.columns))
        basis=BinaryBasis(self.basis.width+dimension)
        for column in columns:basis,_=basis.append(column)
        return IncrementalReductions(self.model,(*self.primes,prime),(*self.dimensions,dimension),
            self.points,columns,basis)


_default=None

def default_reduction_cache():
    global _default
    if _default is None:_default=ReductionCache()
    return _default
