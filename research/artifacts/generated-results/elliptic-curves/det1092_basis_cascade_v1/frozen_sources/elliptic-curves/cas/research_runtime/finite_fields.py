"""Reusable exact point-count tables on rational elliptic families."""
from fractions import Fraction

from .store import FiniteFieldFacts, default_store


def family_traces(A, B, prime, *, a_degree=8, b_degree=12, facts=None, discover=True):
    """Cache all good affine fibres plus infinity; singular traces stay unknown.

    Exact rational coefficients and homogeneous degrees bind the table. Scores
    are intentionally outside this store: different scheduling policies reuse
    the same exact counts. No rank assertion follows from these finite traces.
    """
    from mod2_reduction_independence import _is_prime
    if type(prime) is not int or prime < 5 or not _is_prime(prime):
        raise ValueError("prime >=5 required for a short family table")
    A, B = tuple(map(lambda c: Fraction(str(c)), A)), tuple(map(lambda c: Fraction(str(c)), B))
    if len(A) > a_degree+1 or len(B) > b_degree+1:
        raise ValueError("coefficient degree exceeds homogeneous presentation")
    if a_degree < 0 or b_degree < 0 or 3*a_degree != 2*b_degree:
        raise ValueError("incompatible homogeneous Weierstrass degrees")
    identity = {"A": list(map(str, A)), "B": list(map(str, B)),
                "a_degree": a_degree, "b_degree": b_degree, "base": "P1"}
    def reduce(row):
        return tuple(c.numerator*pow(c.denominator, -1, prime) % prime for c in row)
    a, b = reduce(A), reduce(B)
    def evaluate(row, t):
        value = 0
        for c in reversed(row):
            value = (value*t+c) % prime
        return value
    def build():
        characters = [-1]*prime
        characters[0] = 0
        for x in range(1, prime):
            characters[x*x % prime] = 1
        rows = []
        for t in range(prime+1):
            av = (a[a_degree] if len(a) > a_degree else 0) if t == prime else evaluate(a, t)
            bv = (b[b_degree] if len(b) > b_degree else 0) if t == prime else evaluate(b, t)
            singular = (4*av**3+27*bv**2) % prime == 0
            trace = None if singular else -sum(characters[(x*x*x+av*x+bv) % prime] for x in range(prime))
            rows.append({"parameter": "infinity" if t == prime else t,
                         "a": av, "b": bv, "singular": singular, "trace": trace,
                         "group_order": None if singular else prime+1-trace})
        return {"prime": prime, "fibres": rows, "quadratic_characters": characters,
                "rank_claim": None}
    return (facts or FiniteFieldFacts(default_store())).query(identity, prime, "projective-family-traces",
        build=build if discover else None, version="short-weierstrass-counts-1")
