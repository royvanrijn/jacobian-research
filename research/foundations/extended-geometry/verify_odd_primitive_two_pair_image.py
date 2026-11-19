"""Exact finite controls for odd-composite primitive two-pair Image witnesses.

The unrestricted parameter/order theorem uses the accompanying written proof.
Original contractions below use four-variable monomials, not inverse shears.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, factorial, prod
from pathlib import Path
import platform

if not __debug__:
    raise RuntimeError("This exact certificate requires assertions enabled")


@dataclass
class Balanced:
    degree: int
    # (i,j) denotes xi1^i xi2^(degree-i) z1^j z2^(degree-j).
    terms: dict[tuple[int, int], int]

    def __post_init__(self):
        assert self.degree >= 0
        assert all(0 <= i <= self.degree and 0 <= j <= self.degree
                   and type(c) is int for (i, j), c in self.terms.items())
        self.terms = {p: c for p, c in self.terms.items() if c}


def add(*polynomials):
    degree = polynomials[0].degree
    assert all(p.degree == degree for p in polynomials)
    terms = {}
    for polynomial in polynomials:
        for p, c in polynomial.terms.items():
            terms[p] = terms.get(p, 0) + c
    return Balanced(degree, terms)


def scale(polynomial, coefficient):
    return Balanced(polynomial.degree, {p: c*coefficient for p, c in polynomial.terms.items()})


def multiply(left, right):
    terms = {}
    for (i, j), a in left.terms.items():
        for (k, l), b in right.terms.items():
            p = i+k, j+l
            terms[p] = terms.get(p, 0) + a*b
    return Balanced(left.degree+right.degree, terms)


def power(polynomial, exponent):
    answer = Balanced(0, {(0, 0): 1})
    for _ in range(exponent):
        answer = multiply(answer, polynomial)
    return answer


R = Balanced(1, {(1, 1): 1, (0, 0): 1})
T = Balanced(1, {(1, 1): 1, (0, 0): -1})
Z = Balanced(1, {(1, 0): 1})
W = Balanced(1, {(0, 1): 2})


def construct(r, s):
    assert r >= 1 and s >= 0
    a, b = 2*r+1, 2*r+1+2*s
    A = add(R, Z)
    series = add(*(scale(multiply(power(R, 4*r-j), power(Z, j-1)), comb(4*r, j))
                   for j in range(1, 4*r+1)))
    B = add(scale(multiply(power(R, 4*r), W), 2),
            scale(multiply(power(T, 2), series), -1))
    quotient = multiply(power(A, 3*r), power(B, r))
    if s:
        extra = add(power(R, 2*a*s), multiply(power(T, 2*s), power(A, 4*r*s)))
        quotient = multiply(quotient, extra)
    F = multiply(T, quotient)
    assert F.degree == a*b
    return F, quotient


def contraction(polynomial, marker_power=0):
    """E2(Z^marker_power * polynomial), directly in the original variables."""
    n = polynomial.degree
    return sum(c*factorial(i+marker_power)*factorial(n-i)
               for (i, j), c in polynomial.terms.items() if j == i+marker_power)


def evaluate(polynomial, xi1, xi2, z1, z2):
    n = polynomial.degree
    return sum(c*xi1**i*xi2**(n-i)*z1**j*z2**(n-j)
               for (i, j), c in polynomial.terms.items())


def univariate_multiply(a, b):
    result = {}
    for i, av in a.items():
        for j, bv in b.items():
            result[i+j] = result.get(i+j, 0) + av*bv
    return {p: c for p, c in result.items() if c}


@lru_cache(maxsize=None)
def profile_power(r, s, m):
    H = {1+2*j: (-1)**j*comb(r, j) for j in range(r+1)}
    if s:
        H = univariate_multiply(H, {0: 1, 2*s: 1})
    answer = {0: 1}
    for _ in range(m):
        answer = univariate_multiply(answer, H)
    return answer


def profile_mean(r, s, m):
    return sum((Fraction(c, degree+1) for degree, c in profile_power(r, s, m).items()
                if degree % 2 == 0), Fraction(0))


def profile_mean_binomial(r, s, m):
    if m % 2:
        return Fraction(0)
    extra = range(m+1) if s else range(1)
    return sum((Fraction((-1)**j*comb(r*m, j)*(comb(m, h) if s else 1),
                         m+2*j+2*s*h+1)
                for j in range(r*m+1) for h in extra), Fraction(0))


def angular_coefficient(r, s, m, marker_power):
    if m % 2:
        return Fraction(0)
    index = r*m-marker_power
    if index < 0:
        return Fraction(0)
    return sum((Fraction(c, degree+1)*comb(r*m+2*r*degree, index)
                for degree, c in profile_power(r, s, m).items()), Fraction(0))


def require_equal(actual, expected):
    assert actual == expected


def rejects(check):
    try:
        check()
    except AssertionError:
        return True
    return False


def run():
    assert add(power(T, 2), scale(power(R, 2), -1),
               scale(multiply(Z, W), 2)).terms == {}
    angular_checks = 0
    endpoint_checks = 0
    for r in range(1, 5):
        for s in range(5):
            for m in range(1, 11):
                C = profile_mean(r, s, m)
                assert C == profile_mean_binomial(r, s, m)
                assert (C == 0) if m % 2 else (C > 0)
                assert angular_coefficient(r, s, m, 0) == 0
                assert angular_coefficient(r, s, m, 2*r) == C
                angular_checks += 1
                if m % 2 == 0:
                    # Literal low Taylor coefficients of J_m((1+x)^(2r)).
                    for j in range(1, r*m+1):
                        value = sum((Fraction(c, degree+1)*comb(2*r*(degree+1), j)
                                     for degree, c in profile_power(r, s, m).items()),
                                    Fraction(0))
                        assert value == 0
                        endpoint_checks += 1
    original_records = []
    for r, s, maximum in [(1,0,6),(1,1,4),(1,2,3),(2,0,3),(2,1,2),(3,0,2)]:
        F, quotient = construct(r, s)
        a, b = 2*r+1, 2*r+1+2*s
        d = a*b
        assert evaluate(F, 1, 1, 1, -1) == (-1)**r*2**b != 0
        at_T_zero = 3**(3*r)*2**(r*(4*r+2)+(2*a*s if s else 0))
        assert evaluate(quotient, 1, 1, 1, 1) == at_T_zero != 0
        current = Balanced(0, {(0, 0): 1})
        moments = []
        for m in range(1, maximum+1):
            current = multiply(current, F)
            pure = contraction(current)
            marked = contraction(current, 2*r)
            C = profile_mean(r, s, m)
            expected = factorial(d*m+2*r+1)*C
            assert expected.denominator == 1
            require_equal(pure, 0)
            require_equal(marked, expected)
            if r == 1 and s == 0 and m % 2 == 0:
                closed = Fraction(2**m*factorial(m), prod(m+2*j+1 for j in range(m+1)))
                assert C == closed
                assert contraction(current, 1) == 0  # The old quartic marker is too small here.
            moments.append({"m": m, "expanded_support": len(current.terms),
                            "pure": str(pure), "marked": str(marked),
                            "marked_angular_value": str(C)})
        encoded = [[i, j, c] for (i, j), c in sorted(F.terms.items())]
        original_records.append({"r": r, "s": s, "balanced_degree": d,
                                 "marker_balanced_degree": 2*r,
                                 "F_expanded_support": len(F.terms),
                                 "F_matrix_sha256": sha256(json.dumps(encoded).encode()).hexdigest(),
                                 "R_zero_value": str((-1)**r*2**b),
                                 "F_over_T_at_T_zero": str(at_T_zero), "moments": moments})
    # Reject a changed ORIGINAL coefficient and a false marked normalization.
    F9, _ = construct(1, 0)
    tampered = Balanced(F9.degree, {**F9.terms})
    tampered.terms[(0, 0)] = tampered.terms.get((0, 0), 0)+1
    controls = {
        "altered_original_coefficient_rejected": rejects(
            lambda: require_equal(contraction(tampered), 0)),
        "incorrect_radial_factor_rejected": rejects(
            lambda: require_equal(contraction(power(F9, 2), 2),
                                  factorial(20)*profile_mean(1, 0, 2))),
        "quartic_marker_does_not_claim_new_survival": contraction(power(F9, 2), 1) == 0,
        "new_marker_does_survive": contraction(power(F9, 2), 2) != 0,
    }
    assert all(controls.values())
    return {"status": "passed", "python": platform.python_version(),
            "arithmetic": "integers in original four-variable balanced monomials; Fraction over Q",
            "finite_angular_parameter_cases": angular_checks,
            "literal_endpoint_taylor_coefficients_checked": endpoint_checks,
            "original_contraction_records": original_records,
            "negative_controls": controls,
            "scope": "Original contractions in the six listed finite parameter/order ranges; angular identities at r=1..4,s=0..4,m=1..10. The all-r,s,m theorem uses the written endpoint, radial, primitivity and multiplicity proofs.",
            "independent_external_review": False, "formal_verification": False,
            "worldwide_priority_established": False,
            "odd_prime_degree_result": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = here.parents[2]
    result = run()
    sources = [Path(__file__).resolve(), here / "ODD_PRIMITIVE_TWO_PAIR_IMAGE_WITNESSES.md",
               here / "odd_primitive_two_pair_image_source_ledger.json",
               here / "TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md"]
    result["source_sha256"] = {
        source.relative_to(root).as_posix(): sha256(source.read_bytes()).hexdigest()
        for source in sources
    }
    output = args.output or here / "odd_primitive_two_pair_image_verification.json"
    output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "angular_cases": result["finite_angular_parameter_cases"],
                      "original_parameter_cases": len(result["original_contraction_records"]),
                      "output": str(output)}))


if __name__ == "__main__":
    main()
