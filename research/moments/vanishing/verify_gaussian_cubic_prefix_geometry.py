"""Exact algebra supporting the cubic Gaussian prefix geometry deduction.

The dimension/existence conclusions use the written proof and its classification
input. This checker proves only the explicit finite polynomial identities.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial, prod
from pathlib import Path
import platform


if not __debug__:
    raise RuntimeError("The exact checker requires assertions enabled")


ZERO = (0, 0, 0)


def monomial(powers, coefficient=1):
    return {tuple(powers): Fraction(coefficient)} if coefficient else {}


def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for powers, coefficient in polynomial.items():
            value = result.get(powers, Fraction(0)) + coefficient
            if value:
                result[powers] = value
            else:
                result.pop(powers, None)
    return result


def scale(polynomial, coefficient):
    return {p: c * coefficient for p, c in polynomial.items() if c * coefficient}


def multiply(left, right):
    result = {}
    for a, av in left.items():
        for b, bv in right.items():
            powers = tuple(x + y for x, y in zip(a, b))
            result[powers] = result.get(powers, Fraction(0)) + av * bv
    return {p: c for p, c in result.items() if c}


def power(polynomial, exponent):
    answer = monomial(ZERO)
    for _ in range(exponent):
        answer = multiply(answer, polynomial)
    return answer


def derivative(polynomial, coordinate):
    result = {}
    for powers, coefficient in polynomial.items():
        if powers[coordinate]:
            new_powers = list(powers)
            new_powers[coordinate] -= 1
            result[tuple(new_powers)] = coefficient * powers[coordinate]
    return result


def homogeneous(polynomial, degree):
    return {p: c for p, c in polynomial.items() if sum(p) == degree}


def mean(polynomial):
    answer = Fraction(0)
    for (a, b, c), coefficient in polynomial.items():
        if a == b and c % 2 == 0:
            answer += coefficient * factorial(a) * prod(range(1, c, 2))
    return answer


def invariant(polynomial):
    return mean(multiply(power(homogeneous(polynomial, 1), 3),
                         homogeneous(polynomial, 3)))


def rank(matrix):
    rows = [[Fraction(value) for value in row] for row in matrix]
    if not rows:
        return 0
    pivot_row = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(pivot_row, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        divisor = rows[pivot_row][column]
        rows[pivot_row] = [value / divisor for value in rows[pivot_row]]
        for i in range(len(rows)):
            if i != pivot_row and rows[i][column]:
                coefficient = rows[i][column]
                rows[i] = [a - coefficient * b
                           for a, b in zip(rows[i], rows[pivot_row])]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def encoded(polynomial):
    return [{"powers_Z_W_Y": list(p), "coefficient": str(c)}
            for p, c in sorted(polynomial.items())]


def rejects(check):
    try:
        check()
    except AssertionError:
        return True
    return False


def require_rank(matrix, expected):
    assert rank(matrix) == expected


def require_zero(polynomial):
    assert polynomial == {}


def run():
    z, w, y = (monomial(p) for p in ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    q = add(scale(multiply(z, w), 2), power(y, 2))
    pstar = add(w, multiply(z, y), scale(power(z, 3), Fraction(-1, 2)))
    operators = (
        lambda p: add(multiply(z, derivative(p, 0)),
                      scale(multiply(w, derivative(p, 1)), -1)),
        lambda p: add(multiply(z, derivative(p, 2)),
                      scale(multiply(y, derivative(p, 1)), -1)),
        lambda p: add(multiply(w, derivative(p, 2)),
                      scale(multiply(y, derivative(p, 0)), -1)),
    )
    assert all(op(q) == {} for op in operators)
    tangents = [pstar] + [op(pstar) for op in operators]
    expected = [
        pstar,
        add(scale(w, -1), multiply(z, y), scale(power(z, 3), Fraction(-3, 2))),
        add(power(z, 2), scale(y, -1)),
        add(multiply(z, w), scale(power(y, 2), -1),
            scale(multiply(power(z, 2), y), Fraction(3, 2))),
    ]
    assert tangents == expected
    monomials = [(a, b, d - a - b)
                 for d in range(4) for a in range(d + 1) for b in range(d - a + 1)]
    assert len(monomials) == 20
    tangent_matrix = [[t.get(p, Fraction(0)) for p in monomials] for t in tangents]
    require_rank(tangent_matrix, 4)
    assert invariant(pstar) == -3
    assert invariant(add(w, power(z, 3))) == 6
    assert invariant(add(pstar, scale(power(z, 3), Fraction(1, 2)))) == 0

    weight_controls = []
    for sign in (1, -1):
        linear = [p for p in monomials if sum(p) == 1 and sign * (p[0] - p[1]) > 0]
        cubic = [p for p in monomials if sum(p) == 3 and sign * (p[0] - p[1]) > 0]
        assert len(linear) == 1 and len(cubic) == 4
        count = 0
        for a in linear:
            for b in linear:
                for c in linear:
                    for d in cubic:
                        powers = tuple(sum(values) for values in zip(a, b, c, d))
                        assert sign * (powers[0] - powers[1]) > 0
                        assert mean(monomial(powers)) == 0
                        count += 1
        one_sided = add(*(monomial(p) for p in monomials if sign * (p[0] - p[1]) > 0))
        assert invariant(one_sided) == 0
        assert invariant(homogeneous(one_sided, 1)) == 0
        assert invariant(homogeneous(one_sided, 3)) == 0
        weight_controls.append({"weight_sign": sign, "linear_monomials": len(linear),
                                "cubic_monomials": len(cubic), "products_checked": count})

    # Derivative controls use direct polynomial multiplication, not the shear map.
    moment_jacobian = []
    p_power = monomial(ZERO)
    for n in range(1, 11):
        moment_jacobian.append([n * mean(multiply(monomial(p), p_power)) for p in monomials])
        p_power = multiply(p_power, pstar)
    require_rank(moment_jacobian, 10)

    bad_tangents = tangent_matrix[:]
    bad_tangents[1] = tangent_matrix[0]
    wrong_operator_value = add(multiply(z, derivative(q, 2)),
                               multiply(y, derivative(q, 1)))
    assert wrong_operator_value != {}
    controls = {
        "duplicated_orbit_tangent_rejected": rejects(lambda: require_rank(bad_tangents, 4)),
        "wrong_rotation_sign_rejected": rejects(lambda: require_zero(wrong_operator_value)),
        "mixed_weight_phase_has_nonzero_invariant": invariant(add(w, power(z, 3))) == 6,
    }
    assert all(controls.values())
    return {
        "status": "passed",
        "arithmetic": "fractions.Fraction over Q; sparse original spatial polynomials",
        "python": platform.python_version(),
        "coefficient_space_dimension": len(monomials),
        "separating_invariant_at_Pstar": str(invariant(pstar)),
        "orbit_tangents": [encoded(t) for t in tangents],
        "orbit_tangent_rank": rank(tangent_matrix),
        "first_ten_moment_jacobian_rank": rank(moment_jacobian),
        "one_sided_weight_controls": weight_controls,
        "rejecting_and_boundary_controls": controls,
        "written_conclusions_not_computationally_proved": [
            "The retained classification of every pure ternary cubic",
            "Exceptional orbit dimension and smoothness beyond the displayed tangents",
            "Krull dimension lower bound for moment zero sets",
            "Existence of a nonpure algebraic phase with its first fifteen moments zero",
            "The unspecified later nonzero moment and algebraic coefficient values",
            "The all-order vanishing of later moment differentials",
            "Nonreducedness of the full moment scheme, which uses the written inputs",
        ],
        "not_claimed": ["external review", "formal verification", "universal upper prefix",
                        "sixteenth moment nonzero", "historical novelty"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = here.parents[2]
    result = run()
    sources = [Path(__file__).resolve(), here / "GAUSSIAN_CUBIC_PREFIX_GEOMETRY.md",
               here / "MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md",
               root / "results/gaussian-boundary/manuscript.tex",
               root / "results/gaussian-boundary/affine_classification_proof.tex"]
    result["source_sha256"] = {
        path.relative_to(root).as_posix(): sha256(path.read_bytes()).hexdigest()
        for path in sources
    }
    destination = args.output or here / "gaussian_cubic_prefix_geometry_verification.json"
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "invariant": -3,
                      "orbit_tangent_rank": 4, "moment_jacobian_rank": 10,
                      "receipt": str(destination)}))


if __name__ == "__main__":
    main()
