"""Finite exact QQ controls for the literal regular-sector proof.

The unrestricted global-image and filtered statements are written proofs.
This checker contracts original markers, verifies the spatial lift and
radial recurrences, and computes finite filtered intersections exactly.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
from math import factorial, prod
from pathlib import Path
import platform

import sympy as sp

if not __debug__:
    raise RuntimeError("Run without -O so all exact assertions execute")

Z, W, Y, t = sp.symbols("Z W Y t")
S = Z*W
P0, P = S**2, S**2+Z


def zero(value):
    assert sp.expand(value) == 0, value


def even_moment(power):
    return 0 if power % 2 else factorial(power)//(2**(power//2)*factorial(power//2))


@lru_cache(maxsize=None)
def raw_monomial(a, b, c, n, added_z):
    """Original binomial expansion, followed only by Wick balance."""
    if not added_z:
        return sp.Integer(0) if a != b else sp.Integer(factorial(a+2*n)*even_moment(c))
    answer = sp.Integer(0)
    for j in range(n+1):
        iz, iw = a+2*j+n-j, b+2*j
        if iz == iw:
            answer += sp.binomial(n, j)*factorial(iw)*even_moment(c)
    return answer


def raw_marker(marker, n, added_z=True):
    return sp.expand(sum(coefficient*raw_monomial(a, b, c, n, added_z)
                         for (a, b, c), coefficient
                         in sp.Poly(sp.expand(marker), Z, W, Y).terms()))


def direct_expanded(marker, n):
    result = sp.Integer(0)
    for (a, b, c), coefficient in sp.Poly(sp.expand(marker*P**n), Z, W, Y).terms():
        if a == b:
            result += coefficient*factorial(a)*even_moment(c)
    return sp.expand(result)


def T(marker):
    return sp.expand((W-4*W**2)*marker-(1-2*Z*W**2)*sp.diff(marker, Z)
                     -2*W**3*sp.diff(marker, W))


def radial_series(j, order):
    return sum(sp.Rational(factorial(2*n+j), factorial(n))*t**n for n in range(order+1))


def laurent_coefficients(value):
    result = {}
    for term in sp.Add.make_args(sp.expand(value)):
        power = int(term.as_powers_dict().get(t, 0))
        coefficient = sp.cancel(term/t**power)
        assert not coefficient.has(t)
        result[power] = result.get(power, 0)+coefficient
    return {power: value for power, value in result.items() if value}


def filtered_control(degree, hs, ks):
    pairs = []
    for a in range(degree+1):
        for b in range(a, degree-a+1):
            scalar = sp.Rational(1, 2**b*factorial(b-a))
            pairs.append((sp.expand(scalar*t**(-a)*hs[b]),
                          sp.expand(-scalar*t**(-a)*ks[b])))
    constraints = sorted({("f", e) for a, _ in pairs for e in laurent_coefficients(a)}
                         | {("b", e) for _, b in pairs for e in laurent_coefficients(b) if e < 0})
    rows = [[laurent_coefficients(a if kind == "f" else b).get(e, 0)
             for a, b in pairs] for kind, e in constraints]
    matrix = sp.Matrix(rows) if rows else sp.zeros(0, len(pairs))
    values = []
    for vector in matrix.nullspace():
        value = sp.expand(sum(vector[i]*pairs[i][1] for i in range(len(pairs))))
        assert all(e >= 0 for e in laurent_coefficients(value))
        values.append(value)
    bound = (degree-1)//2 if degree >= 1 else -1
    assert all(not value or sp.Poly(value, t).degree() <= bound for value in values)
    columns = [sp.Matrix([value.coeff(t, j) for j in range(max(0, bound+1))])
               for value in values]
    rank = sp.Matrix.hstack(*columns).rank() if columns and bound >= 0 else 0
    assert rank == max(0, bound+1)
    return {"marker_degree": degree, "polynomial_sector_dimension": rank,
            "maximum_polynomial_degree": bound}


def checks():
    zero((1-2*Z*W**2)*sp.diff(P, Z)+2*W**3*sp.diff(P, W)-1)
    markers = [1, Z, W, S, W**2, Z*W**2, Y, Y**2, W*Y**2, Z**2*W**3]
    for marker in markers:
        for n in range(7):
            zero(raw_marker(T(marker), n)-n*raw_marker(marker, n-1) if n else raw_marker(T(marker), 0))
    for marker in markers:
        for n in range(5):
            zero(raw_marker(marker, n)-direct_expanded(marker, n))
    order, max_j = 10, 10
    fs = [radial_series(j, order) for j in range(max_j+2)]
    zero(sp.expand(2*t*fs[1]-fs[0]+1).series(t, 0, order+1).removeO())
    for j in range(1, max_j+1):
        zero(sp.expand(j*fs[j-1]+2*t*fs[j+1]-fs[j]).series(t, 0, order+1).removeO())
    hs, ks = [sp.Integer(1), sp.Integer(1)], [sp.Integer(0), sp.Integer(1)]
    for j in range(1, max_j):
        hs.append(sp.expand(hs[-1]-2*j*t*hs[-2]))
        ks.append(sp.expand(ks[-1]-2*j*t*ks[-2]))
    for j in range(1, max_j+1):
        assert sp.Poly(hs[j], t).degree() <= j//2
        assert sp.Poly(ks[j], t).degree() <= (j-1)//2
        value = sp.expand(hs[j]*radial_series(0, order+j)-ks[j])
        value = sp.series(value/(2*t)**j, t, 0, order+1).removeO()
        zero(value-fs[j])
    # Exact original angular formula, with the exponential factorials retained.
    for a in range(5):
        for b in range(5):
            for c in range(5):
                for n in range(7):
                    expected = 0
                    shift = b-a
                    if shift >= 0 and n >= shift:
                        expected = (even_moment(c)*sp.Rational(factorial(n)*factorial(2*(n-shift)+b),
                                                             factorial(shift)*factorial(n-shift)))
                    zero(raw_monomial(a, b, c, n, True)-expected)
    # The Laurent regularizations have the predicted core multiplier polynomials.
    for j in range(9):
        for n in range(9):
            ratio = 2**j*prod(2*n+2*i-1 for i in range(1, j+1))
            assert sp.Rational(factorial(2*(n+j))*factorial(n),
                               factorial(n+j)*factorial(2*n)) == ratio
    lifts = []
    marker = 1-2*W
    for n in range(7):
        assert sp.Poly(marker, Z, W, Y).total_degree() == 2*n+1
        for m in range(11):
            zero(raw_marker(marker, m)-factorial(m)*(1 if m == n else 0))
        lifts.append({"target_degree": n, "least_marker_degree": 2*n+1,
                      "marker": str(marker)})
        marker = T(marker)
    filtered = [filtered_control(degree, hs, ks) for degree in range(11)]
    return {"status": "passed", "coefficient_field": "QQ",
            "original_lift_markers": len(markers), "lift_moment_orders": [0, 6],
            "independent_expanded_Wick_orders": [0, 4],
            "radial_recurrence_indices": [0, max_j], "series_order": order,
            "angular_formula_exponents_each": [0, 4], "angular_formula_moments": [0, 6],
            "Laurent_regularization_indices": [0, 8],
            "polynomial_lifts": lifts, "filtered_intersections": filtered,
            "scope": "Finite exact controls. Written proofs establish the unrestricted global module and all-degree literal filtration.",
            "independent_replay": False, "formal_verification": False, "external_review": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    result.update(python=platform.python_version(), sympy=sp.__version__,
                  proof_sha256=digest(directory/"GAUSSIAN_ANGULAR_CORE_LITERAL_REGULAR_SECTOR.md"),
                  checker_sha256=digest(Path(__file__)),
                  angular_core_import_sha256=digest(directory/"GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md"))
    report = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(report, encoding="utf-8", newline="\n")
    print(report, end="")


if __name__ == "__main__":
    main()
