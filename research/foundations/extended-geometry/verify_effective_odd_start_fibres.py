#!/usr/bin/env python3
"""Exact finite controls for the effective fixed-coordinate fibre theorem.

This checks the displayed bounded factor searches, one complete finite box,
ramified pole twists, and multiplication matrices retaining all jets.  It is
not an implementation of the general Orbit decision algorithm or of the
whole classification over arbitrary number fields.
"""

from __future__ import annotations

if not __debug__:
    raise RuntimeError("This verifier requires assertions; do not use -O or -OO.")

import argparse
from hashlib import sha256
from itertools import product
import json
from math import factorial
from pathlib import Path
import platform

import sympy as sp


x, y, z = sp.symbols("x y z")


def same(left, right):
    return sp.cancel(left - right) == 0


def fraction(value, variable=x):
    numerator, denominator = sp.cancel(value).as_numer_denom()
    return (sp.Poly(numerator, variable, domain=sp.QQ),
            sp.Poly(denominator, variable, domain=sp.QQ))


def degree(value):
    numerator, denominator = fraction(value)
    return max(numerator.degree(), denominator.degree())


def at_infinity(value):
    numerator, denominator = fraction(value)
    if numerator.degree() > denominator.degree():
        return sp.oo
    if numerator.degree() < denominator.degree():
        return sp.Integer(0)
    return numerator.LC() / denominator.LC()


def marked_component(u):
    assert same(u.subs(x, 0), 0)
    assert same(sp.diff(u, x).subs(x, 0), 1)
    equation = sp.cancel(u + u.subs(x, y)).as_numer_denom()[0]
    factors = sp.factor_list(equation, x, y)[1]
    marked = [factor for factor, _ in factors
              if factor.subs({x: 0, y: 0}) == 0]
    assert len(marked) == 1
    component = sp.Poly(marked[0], x, y, domain=sp.QQ).monic().as_expr()
    assert sp.diff(component, x).subs({x: 0, y: 0}) != 0
    assert sp.diff(component, y).subs({x: 0, y: 0}) != 0
    return component


def vanishes_on(numerator, component):
    remainder = sp.reduced(sp.expand(numerator), [component],
                           x, y, domain=sp.QQ)[1]
    return remainder == 0


def reciprocal(eta, component):
    assert same(eta.subs(x, 0), 1)
    numerator, denominator = sp.cancel(
        eta * eta.subs(x, y) - 1).as_numer_denom()
    # Denominator clearing must not turn a function undefined identically on
    # the selected component into an accepted norm identity.
    assert not vanishes_on(denominator, component)
    return vanishes_on(numerator, component)


def canonical_search(u):
    """The proof's finite search over Q, for the explicit rational examples."""
    component = marked_component(u)
    map_degree = degree(u)
    numerator, denominator = fraction(u)
    c = at_infinity(u)
    fibre = (denominator.as_expr() if c is sp.oo else
             numerator.as_expr() + c * denominator.as_expr())
    factors = []
    for factor, _ in sp.factor_list(fibre, x)[1]:
        if factor.subs(x, 0) != 0:
            factors.append(sp.Poly(factor / factor.subs(x, 0),
                                   x, domain=sp.QQ))
    accepted, tested = [], 0
    ranges = [range(map_degree // factor.degree() + 1)
              for factor in factors]
    for powers in product(*ranges):
        size = sum(power * factor.degree()
                   for power, factor in zip(powers, factors))
        if size == 0 or map_degree % size != 0:
            continue
        candidate = sp.prod(factor.as_expr()**power
                            for factor, power in zip(factors, powers))
        candidate = sp.expand(candidate)
        tested += 1
        if reciprocal(candidate, component):
            accepted.append(candidate)
    accepted.sort(key=lambda value: sp.degree(value, x))
    generator = accepted[0] if accepted else None
    if generator is not None:
        smallest = sp.degree(generator, x)
        assert sum(sp.degree(value, x) == smallest
                   for value in accepted) == 1
        for value in accepted:
            exponent, remainder = divmod(sp.degree(value, x), smallest)
            assert remainder == 0 and same(value, generator**exponent)
    return generator, tested, accepted, component


def factor_orders(value, variable=x):
    numerator, denominator = fraction(value, variable)
    orders = {}
    for sign, polynomial in ((1, numerator), (-1, denominator)):
        for factor, multiplicity in sp.factor_list(
                polynomial.as_expr(), variable)[1]:
            key = sp.Poly(factor, variable, domain=sp.QQ).monic().as_expr()
            orders[key] = orders.get(key, 0) + sign * multiplicity
    return {factor: order for factor, order in orders.items() if order}


def pole_admissible(eta, u, a):
    u_orders = factor_orders(u)
    for factor, order in factor_orders(eta).items():
        if order >= 0:
            continue
        pole_order = -u_orders.get(factor, 0)
        if pole_order <= 0 or -order > 1 + pole_order // (a + 2):
            return False
    return True


def jet_equals(eta, target, modulus):
    numerator, denominator = fraction(eta)
    jj = sp.Poly(modulus, x, domain=sp.QQ)
    assert sp.gcd(denominator, jj).degree() == 0
    difference = numerator - sp.Poly(target, x, domain=sp.QQ) * denominator
    return difference.rem(jj).is_zero


def polynomial_square(value):
    coefficient, factors = sp.factor_list(value, x)
    # Only the parity is needed in the control with nonsquare discriminant.
    return all(exponent % 2 == 0 for _, exponent in factors) and (
        sp.sqrt(coefficient).is_Rational is True)


def multiplication_matrix(value, modulus, variable=z):
    jj = sp.Poly(modulus, variable, domain=sp.QQ)
    size = jj.degree()
    numerator, denominator = fraction(value, variable)
    assert sp.gcd(denominator, jj).degree() == 0
    reduced = (numerator * sp.invert(denominator, jj)).rem(jj)
    columns = []
    for index in range(size):
        remainder = (reduced * sp.Poly(variable**index, variable)).rem(jj)
        columns.append(sp.Matrix([remainder.nth(i) for i in range(size)]))
    return sp.Matrix.hstack(*columns)


def coordinates(value, modulus, variable=z):
    numerator, denominator = fraction(value, variable)
    jj = sp.Poly(modulus, variable, domain=sp.QQ)
    assert sp.gcd(denominator, jj).degree() == 0
    remainder = (numerator * sp.invert(denominator, jj)).rem(jj)
    return sp.Matrix([remainder.nth(i) for i in range(jj.degree())])


def matrix_json(matrix):
    return [[str(entry) for entry in row] for row in matrix.tolist()]


def verify():
    counts = {
        "canonical_search_examples": 0,
        "canonical_factor_candidates": 0,
        "finite_box_vectors": 0,
        "finite_box_reciprocal_vectors": 0,
        "ramified_twist_vectors": 0,
        "ramified_twist_accepted_vectors": 0,
        "orbit_algebras": 0,
        "orbit_exponents_checked": 0,
        "restriction_of_scalars_examples": 0,
        "negative_controls": 0,
    }
    details = {}

    quadratic = 1 + x + x**2 / 8
    smaller = 1 + x
    examples = [
        ("linear_polynomial_coordinate", x, None),
        ("nonlinear_polynomial_coordinate", x + x**2, None),
        ("mobius_coordinate", x / (1 + x), 1 + 2*x),
        ("quadratic_nonrational_return",
         2 * (quadratic - 1) / (quadratic + 1), quadratic),
        ("nonminimal_initial_generator",
         (smaller**2 - 1) / (smaller**2 + 1), smaller),
    ]
    search_records = []
    for name, u, expected in examples:
        found, tested, accepted, component = canonical_search(u)
        assert (found is None) == (expected is None)
        if expected is not None:
            assert same(found, expected)
        if name == "quadratic_nonrational_return":
            assert sp.degree(component, y) == 2
            discriminant = sp.discriminant(component, y)
            assert not polynomial_square(discriminant)
            counts["negative_controls"] += 1
        if name == "nonminimal_initial_generator":
            assert any(same(value, smaller**2) for value in accepted)
            assert sp.degree(found, x) < sp.degree(smaller**2, x)
            counts["negative_controls"] += 1
        counts["canonical_search_examples"] += 1
        counts["canonical_factor_candidates"] += tested
        search_records.append({
            "name": name, "coordinate": str(sp.cancel(u)),
            "map_degree": int(degree(u)), "candidates": tested,
            "canonical_generator": None if found is None else str(found),
            "marked_component": str(component),
        })
    details["canonical_searches"] = search_records

    # Enumerate the entire proved box for this no-generator example.
    u = x / (1 - x**2)
    found, tested, _, component = canonical_search(u)
    assert found is None and tested == 0 and same(component, x + y)
    a, support_size, map_degree = 1, 2, int(degree(u))
    bb = 1 + map_degree // (a + 2)
    hh = factorial(support_size) * bb * (2*map_degree)**(support_size - 1)
    assert bb == 1 and hh == 8
    assert factor_orders(u) == {x: 1, x - 1: -1, x + 1: -1}
    assert same(sp.diff(u, x), (1 + x**2)/(1 - x**2)**2)
    modulus = 1 + x**2
    norm_vectors, final_vectors = [], []
    for n_minus, n_plus in product(range(-1, hh + 1), repeat=2):
        eta = sp.cancel((1 - x)**n_minus * (1 + x)**n_plus)
        counts["finite_box_vectors"] += 1
        assert pole_admissible(eta, u, a)
        if reciprocal(eta, component):
            norm_vectors.append([n_minus, n_plus])
            assert not jet_equals(eta, -1, modulus)
            if jet_equals(eta, -1, modulus):
                final_vectors.append([n_minus, n_plus])
    assert norm_vectors == [[-1, 1], [0, 0], [1, -1]]
    assert final_vectors == [] and counts["finite_box_vectors"] == 100
    counts["finite_box_reciprocal_vectors"] = len(norm_vectors)
    valid_eta = (1 + x)/(1 - x)
    assert reciprocal(valid_eta, component)
    assert not reciprocal(valid_eta, x*y - 1)
    assert not reciprocal(valid_eta, (x + y)*(x*y - 1))
    counts["negative_controls"] += 2
    details["finite_box"] = {
        "N": map_degree, "s": support_size, "a": a, "B": bb, "H": hh,
        "lower_bounds": [-1, -1], "tested_vectors": 100,
        "reciprocal_vectors": norm_vectors,
        "jet_modulus": str(modulus), "final_vectors": final_vectors,
        "other_correspondence_component": "x*y - 1",
    }

    # The maximum ramification, rather than an arbitrary point in a fibre,
    # creates the asymmetric exponent interval in this example.
    rr = 1 + x + x**2
    alpha, inverse_alpha = sp.Rational(3, 4), sp.Rational(4, 3)
    uu = -(z**2 - 1)/(24*(z - alpha)*(z - inverse_alpha))
    u = sp.cancel(uu.subs(z, rr))
    assert same(uu.subs(z, 1/z), -uu)
    found, tested, _, component = canonical_search(u)
    assert same(found, rr) and degree(u) == 4
    counts["canonical_search_examples"] += 1
    counts["canonical_factor_candidates"] += tested
    pole_values = factor_orders(uu, z)
    assert pole_values[z - alpha] == pole_values[z - inverse_alpha] == -1
    ramification = [
        max(exponent for _, exponent in sp.factor_list(rr - value, x)[1])
        for value in (alpha, inverse_alpha)
    ]
    assert ramification == [2, 1]
    q_values = [int(sp.floor(sp.Rational(1, 3) + sp.Rational(1, ee)))
                for ee in ramification]
    assert q_values == [0, 1]
    assert fraction(uu, z)[1].eval(0) != 0  # k_min = 0.
    twist = sp.cancel((rr - alpha)/(1 - alpha*rr))
    assert reciprocal(twist, component) and reciprocal(rr, component)
    accepted = []
    for n, j in product(range(-2, 4), range(-2, 5)):
        eta = sp.cancel(twist**n * rr**j)
        actual = pole_admissible(eta, u, 1)
        expected = n in (0, 1) and j >= 0
        assert actual == expected, (n, j, factor_orders(eta))
        counts["ramified_twist_vectors"] += 1
        if actual:
            accepted.append([n, j])
    assert len(accepted) == 10
    counts["ramified_twist_accepted_vectors"] = len(accepted)
    assert not pole_admissible(1/twist, u, 1)
    assert not pole_admissible(1/rr, u, 1)
    counts["negative_controls"] += 2
    details["ramified_twists"] = {
        "canonical_generator": str(rr), "generator_search_candidates": tested,
        "outer_coordinate": str(uu), "pole_pair": ["3/4", "4/3"],
        "maximum_ramification": ramification, "pole_orders": [1, 1],
        "q_values": q_values, "k_min": 0,
        "tested_n_range": [-2, 3], "tested_j_range": [-2, 4],
        "accepted_vectors": accepted,
    }

    # Full quotient-algebra vectors, including nilpotent directions.
    orbit_records = []
    modulus = z**2 + 1
    mm = multiplication_matrix(z, modulus)
    start, target = coordinates(1, modulus), coordinates(-1, modulus)
    assert mm**2 == -sp.eye(2) and mm**4 == sp.eye(2)
    assert mm != sp.eye(2) and mm**2 != sp.eye(2)
    for j in range(25):
        assert (mm**j * start == target) == (j % 4 == 2)
        counts["orbit_exponents_checked"] += 1
    orbit_records.append({"modulus": str(modulus), "matrix": matrix_json(mm),
                          "solution": "j = 2 modulo 4"})
    counts["orbit_algebras"] += 1

    modulus = (z - 1)**2
    mm = multiplication_matrix(z, modulus)
    start = coordinates(1, modulus)
    target = coordinates(1 + 7*(z - 1), modulus)
    nilpotent = mm - sp.eye(2)
    assert nilpotent**2 == sp.zeros(2)
    assert nilpotent*start != sp.zeros(2, 1)
    assert target == start + 7*nilpotent*start
    # (I+N)^j = I+jN when N^2=0.  The nonzero N*start makes j=7
    # the sole integer solution; the loop checks the concrete realization.
    for j in range(25):
        assert mm**j == sp.eye(2) + j*nilpotent
        assert (mm**j*start == target) == (j == 7)
        counts["orbit_exponents_checked"] += 1
    orbit_records.append({"modulus": str(modulus), "matrix": matrix_json(mm),
                          "solution": "j = 7 only"})
    counts["orbit_algebras"] += 1

    modulus = (z + 1)**2
    mm = multiplication_matrix(z, modulus)
    start, target = coordinates(1, modulus), coordinates(-1, modulus)
    nilpotent = mm + sp.eye(2)
    assert nilpotent**2 == sp.zeros(2)
    assert sp.Matrix.hstack(start, nilpotent*start).det() != 0
    # (-I+N)^j=(-1)^j(I-jN).  Equality with -start needs odd j
    # from the first basis coefficient and j=0 from the second: impossible.
    for j in range(25):
        assert mm**j == (-1)**j*(sp.eye(2) - j*nilpotent)
        assert mm**j*start != target
        if j % 2:
            assert (-1)**j == -1  # Point values alone falsely accept.
        counts["orbit_exponents_checked"] += 1
    orbit_records.append({"modulus": str(modulus), "matrix": matrix_json(mm),
                          "solution": "empty; odd point-value matches fail jets"})
    counts["orbit_algebras"] += 1
    counts["negative_controls"] += 1

    modulus = z**2
    mm = multiplication_matrix(z, modulus)
    start = coordinates(1, modulus)
    assert mm.det() == 0 and mm**2 == sp.zeros(2)
    for j in range(10):
        assert (mm**j*start == start) == (j == 0)
        assert mm**j*coordinates(z, modulus) != -start
        counts["orbit_exponents_checked"] += 1
    counts["orbit_algebras"] += 1
    counts["negative_controls"] += 2
    orbit_records.append({"modulus": str(modulus), "matrix": matrix_json(mm),
                          "solution": "nonunit R: target 1 only at j=0; "
                                      "nonunit base can never yield -1"})
    details["orbit_controls"] = orbit_records

    # Restriction of scalars must retain both number-field and jet directions.
    aa, epsilon = sp.symbols("a epsilon")
    relations = sp.groebner([aa**2 - 2, epsilon**2],
                           aa, epsilon, domain=sp.QQ)
    basis = [sp.Integer(1), aa, epsilon, aa*epsilon]

    def vector(value):
        remainder = relations.reduce(sp.expand(value))[1]
        polynomial = sp.Poly(remainder, aa, epsilon, domain=sp.QQ)
        return sp.Matrix([polynomial.coeff_monomial(term) for term in basis])

    algebra_element = aa + epsilon
    mm = sp.Matrix.hstack(*(vector(algebra_element*term) for term in basis))
    expected_matrix = sp.Matrix([[0, 2, 0, 0], [1, 0, 0, 0],
                                 [1, 0, 0, 2], [0, 1, 1, 0]])
    assert mm == expected_matrix
    start = vector(1)
    target = vector(2*aa + 6*epsilon)
    assert vector(algebra_element**3) == target == sp.Matrix([0, 2, 6, 0])
    assert mm**3*start == target
    assert (mm**2 - 2*sp.eye(4))**2 == sp.zeros(4)
    assert mm**2 - 2*sp.eye(4) != sp.zeros(4)  # The algebra is not reduced.
    augmented = sp.diag(mm, sp.Integer(1))
    augmented_start = start.col_join(sp.Matrix([1]))
    augmented_target = target.col_join(sp.Matrix([1]))
    assert augmented**3*augmented_start == augmented_target
    assert sp.Matrix.hstack(2*target, target).rank() == 1
    assert sp.Matrix.hstack((2*target).col_join(sp.Matrix([1])),
                           augmented_target).rank() == 2
    for j in range(10):
        image = mm**j*start
        assert image == vector(algebra_element**j)
        augmented_image = augmented**j*augmented_start
        line_hit = sp.Matrix.hstack(augmented_image, augmented_target).rank() == 1
        assert line_hit == (image == target)
        counts["orbit_exponents_checked"] += 1
    counts["restriction_of_scalars_examples"] = 1
    counts["negative_controls"] += 1
    details["restriction_of_scalars"] = {
        "field": "Q(sqrt(2))", "jet_relation": "epsilon^2 = 0",
        "rational_basis": ["1", "sqrt(2)", "epsilon", "sqrt(2)*epsilon"],
        "multiplier": "sqrt(2) + epsilon",
        "matrix": matrix_json(mm), "target": ["0", "2", "6", "0"],
        "witness_exponent": 3, "augmented_matrix_dimension": 5,
        "unaugmented_line_accepts_wrong_scalar": True,
        "augmented_line_rejects_wrong_scalar": True,
    }
    return counts, details


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path,
                        default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    directory = Path("research/foundations/extended-geometry")
    ledger_path = directory / "effective_odd_start_fibres_source_ledger.json"
    ledger = json.loads((repo / ledger_path).read_text())
    for source in ledger["sources"]:
        assert digest(repo / source["path"]) == source["sha256"], source["path"]
    counts, details = verify()
    proof = directory / "EFFECTIVE_ODD_START_FIBRE_CLASSIFICATION.md"
    checker = directory / Path(__file__).name
    receipt = {
        "schema": "effective_odd_start_fibres_verification_v1",
        "status": "PASS",
        "arithmetic": "exact SymPy rational polynomial and matrix arithmetic",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "counts": counts,
        "details": details,
        "source_pins_checked": len(ledger["sources"]),
        "packet_pins": {str(path): digest(repo / path)
                        for path in (proof, checker, ledger_path)},
        "scope": "Finite exact controls for the written terminating "
                 "classification. The canonical-unit and pole-bound inputs "
                 "come from the retained theorem; the rational Orbit "
                 "decision subroutine is imported from its primary source.",
        "not_claimed": [
            "Implementation of the general Orbit decision algorithm",
            "End-to-end implementation over arbitrary number fields",
            "A finite-control proof of the universal classification",
            "A polynomial-time bound for the complete classification",
            "Classification of all endpoint coordinates simultaneously",
            "External human review, formal verification or historical priority",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "counts": counts,
                      "receipt": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
