"""Finite-field rank bounds and regulator compatibility for elliptic surfaces.

The source-specific arithmetic replayer establishes good reduction, the full
Frobenius polynomial and fibre data. This layer handles exact rank extraction
and determinant compatibility without assuming a particular construction.
"""

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt

from .store import digest


def rational_square(value):
    value = Fraction(str(value))
    return value >= 0 and isqrt(value.numerator)**2 == value.numerator and isqrt(value.denominator)**2 == value.denominator


def squareclass(value):
    from sympy import factorint
    value = Fraction(str(value))
    if value <= 0:
        raise ValueError("positive height/regulator required")
    result = 1
    # A rational and numerator*denominator have the same squareclass.
    for p, exponent in factorint(value.numerator * value.denominator).items():
        if exponent % 2:
            result *= int(p)
    return result


def frobenius_invariants(coefficients, prime, *, expected_degree, moments=None):
    """Extract arithmetic and geometric MW bounds from det(Z-Frob|H1).

    Geometric rank counts ALL root-of-unity factors of F(pZ); it does not
    truncate a cyclotomic search. A rank bound is only useful after the input
    polynomial and its removal of the trivial lattice have been verified.
    """
    from sympy import Poly, Rational, Symbol, factor_list, isprime
    if type(prime) is not int or not 5 <= prime < 2**64 or not isprime(prime):
        raise ValueError("good finite-field prime must be a certified prime >=5")
    coeffs = tuple(Fraction(str(c)) for c in coefficients)
    if type(expected_degree) is not int or expected_degree < 0 or len(coeffs) != expected_degree+1 or coeffs[-1] != 1:
        raise ValueError("incomplete or nonmonic Frobenius polynomial")
    if any(c.denominator != 1 for c in coeffs):
        raise ValueError("Frobenius polynomial must be integral")
    T = Symbol("T")
    L = Poly(sum(Rational(c.numerator,c.denominator)*T**i for i,c in enumerate(reversed(coeffs))),T,domain="QQ")
    if moments is not None:
        first = -Fraction(L.nth(1))
        second = Fraction(L.nth(1))**2-2*Fraction(L.nth(2))
        if tuple(map(lambda v: Fraction(str(v)),moments)) != (first,second):
            raise ValueError("Frobenius orientation/moment mismatch")
    quotient, arithmetic_rank = L, 0
    linear = Poly(1-prime*T,T,domain="QQ")
    while quotient.eval(Rational(1,prime)) == 0:
        quotient,remainder=divmod(quotient,linear)
        if not remainder.is_zero:
            raise ArithmeticError("nonexact Frobenius rank division")
        arithmetic_rank += 1
    leading = Fraction(quotient.eval(Rational(1,prime)))
    if leading <= 0:
        raise ValueError("nonpositive normalized L-value")
    scaled = Poly(sum(Rational(c.numerator,c.denominator)*prime**i*T**i for i,c in enumerate(coeffs)),T,domain="QQ")
    geometric_rank = sum(factor.degree()*multiplicity for factor,multiplicity in factor_list(scaled)[1] if factor.is_cyclotomic)
    if arithmetic_rank > geometric_rank:
        raise ArithmeticError("arithmetic rank exceeds geometric rank")
    return {"arithmetic_rank_upper": arithmetic_rank, "geometric_rank_upper": int(geometric_rank),
            "L_coefficients_low_to_high": [str(c) for c in reversed(coeffs)],
            "L_star": str(leading)}


@dataclass(frozen=True)
class Surface:
    """Exact d(t)y²=x³+A(t)x+B(t), including twist constant and base labels."""
    A: tuple[str, ...]
    B: tuple[str, ...]
    d: tuple[str, ...] = ("1",)
    base: str = "QQ(t)"

    def __post_init__(self):
        for key in ("A","B","d"):
            row=list(map(lambda c:str(Fraction(str(c))),getattr(self,key)))
            while len(row)>1 and row[-1]=="0":row.pop()
            if not row:raise ValueError("missing surface polynomial")
            object.__setattr__(self,key,tuple(row))
        if self.d == ("0",) or self.base != "QQ(t)":
            raise ValueError("unsupported surface base or zero twist")

    @property
    def key(self):
        return digest(self)


@dataclass(frozen=True)
class VerifiedReduction:
    surface_key: str
    prime: int
    arithmetic_rank_upper: int
    geometric_rank_upper: int
    regulator_if_rank_equality: str
    witness_hash: str
    height_specialization_isometry: bool


def replay_reduction(surface, witness, *, replay):
    """Bind a cheap arithmetic witness to an exact model before extracting facts.

    ``replay`` must check the full polynomial, all closed-place Tamagawa
    factors (including infinity), chi, and good smooth surface reduction.
    A hash or two point-count moments alone do not prove a full polynomial.
    """
    if witness.get("surface_key") != surface.key:
        raise ValueError("reduction witness belongs to a different surface")
    if replay(surface,witness) is not True:
        raise ArithmeticError("finite-field proof replay failed")
    chi, tamagawa = witness["chi"], witness["tamagawa_product"]
    if type(chi) is not int or chi < 1 or type(tamagawa) is not int or tamagawa < 1:
        raise ValueError("missing exact fibre/BSD normalization")
    if witness.get("good_reduction") is not True:
        raise ValueError("rank specialization requires verified good surface reduction")
    if witness.get("nonisotrivial") is not True:
        raise ValueError("height/regulator gate requires a nonisotrivial surface")
    data=frobenius_invariants(witness["frobenius_coefficients"],witness["prime"],
        expected_degree=witness["expected_L_degree"],moments=witness.get("moments"))
    regulator=Fraction(witness["prime"])**(chi-1)*Fraction(data["L_star"])/tamagawa
    return VerifiedReduction(surface.key,witness["prime"],data["arithmetic_rank_upper"],
        data["geometric_rank_upper"],str(regulator),digest(witness),
        witness.get("height_specialization_isometry") is True)


def pre_search_gate(surface, reductions, *, candidate_rank, candidate_regulator=None):
    """Rank bound, then two-prime regulator test BEFORE generating sections.

    Regulator comparison is valid only when a hypothetical rank-r subgroup
    forces reduction rank=analytic rank=r at BOTH distinct primes. In higher
    rank a mismatch rules out rank r, hence gives upper bound r-1. It does not
    prove rank zero unless r=1. Compatible determinants remain UNKNOWN.
    """
    if type(candidate_rank) is not int or candidate_rank < 1:
        raise ValueError("positive candidate rank required")
    reductions=tuple(reductions)
    if any(not isinstance(row,VerifiedReduction) or row.surface_key != surface.key for row in reductions):
        raise ValueError("misbound or unverified finite reduction")
    if len({row.prime for row in reductions}) != len(reductions):
        raise ValueError("repeated primes cannot supply independent reductions")
    arithmetic=min((r.arithmetic_rank_upper for r in reductions),default=None)
    geometric=min((r.geometric_rank_upper for r in reductions),default=None)
    # A reduction may already improve the requested target. Test equality at
    # that smaller bound too: two rank-one reductions can prove rank zero even
    # when the original search asked for two or more sections.
    rank_hypothesis=min(candidate_rank,arithmetic) if arithmetic is not None else candidate_rank
    eligible=[r for r in reductions if rank_hypothesis>0 and
              r.arithmetic_rank_upper==rank_hypothesis and r.height_specialization_isometry]
    comparisons=[]
    if eligible:
        first=eligible[0]
        for second in eligible[1:]:
            ratio=Fraction(first.regulator_if_rank_equality)/Fraction(second.regulator_if_rank_equality)
            comparisons.append({"primes":[first.prime,second.prime],"ratio":str(ratio),
                                "compatible":rational_square(ratio)})
    incompatible=any(not row["compatible"] for row in comparisons)
    if incompatible:
        arithmetic=min(arithmetic,rank_hypothesis-1)
    height_tests=[]
    if candidate_regulator is not None:
        candidate=Fraction(str(candidate_regulator))
        if candidate<=0:raise ValueError("positive candidate determinant required")
        for row in reductions:
            if row.arithmetic_rank_upper!=candidate_rank or not row.height_specialization_isometry:continue
            ratio=candidate/Fraction(row.regulator_if_rank_equality)
            height_tests.append({"prime":row.prime,"ratio":str(ratio),"compatible":rational_square(ratio)})
    excluded=(arithmetic is not None and arithmetic<candidate_rank)
    candidate_excluded=any(not row["compatible"] for row in height_tests)
    return {"schema":"elliptic-surface.pre-search-regulator.v1","surface_key":surface.key,
            "candidate_rank":candidate_rank,"arithmetic_rank_upper":arithmetic,
            "geometric_rank_upper":geometric,"rank_target_excluded":excluded,
            "candidate_regulator_excluded":candidate_excluded,
            "section_search_eligible":not(excluded or candidate_excluded),
            "two_prime_regulator_test":"INCOMPATIBLE" if incompatible else "COMPATIBLE" if comparisons else "NOT_APPLICABLE",
            "comparisons":comparisons,"height_tests":height_tests,
            "witness_hashes":[row.witness_hash for row in reductions],
            "status":"PROVED_EXCLUSION" if excluded or candidate_excluded else "UNKNOWN"}
