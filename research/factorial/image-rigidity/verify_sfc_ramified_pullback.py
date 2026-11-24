#!/usr/bin/env python3
"""Exact polynomial and actual-moment certificates for constrained pullbacks.

Standard library only. The unbounded statements have written residue proofs;
this checker records its complete finite inputs and exact arithmetic outputs.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import comb, factorial, gcd, prod
from pathlib import Path

if not __debug__:
    raise SystemExit("Refusing optimized Python: assertions must remain enabled.")

HERE = Path(__file__).resolve().parent
NOTE_NAME = "SFC_RAMIFIED_PULLBACK_CLOSURE.md"
NOTE_PIN = "eb444a8f02feb92834336589024e0d920b4eba23970582291af406c11bd7f519"
CONTACT_NAME = "research/factorial/image-rigidity/SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md"
CONTACT_PIN = "409ca1ed8a1ef3218fb07673d14830d04cf10562f5714c88bd8de8e9a681b417"
SEED_PINS = {
    "POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md": "40278040a10d4813c473e8256b868e4fb6d42f20912972b801136398436609d6",
    "WEIGHTED_RADIAL_LINEAR_SFC_TRANSFER.md": "f9ee0dd06d5e6bc85a5b8d35aacb6d55f5c8ec18e483a59c00d3275177495645",
}


def trim(a):
    a = [Q(x) for x in a]
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a or [Q(0)]


def add(a, b):
    out = [Q(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return trim(out)


def scale(a, k):
    return trim([Q(k) * x for x in a])


def mul(a, b, cap=None):
    size = len(a) + len(b) - 1
    if cap is not None:
        size = min(size, cap + 1)
    out = [Q(0)] * size
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j < size:
                out[i + j] += x * y
    return trim(out)


def power(a, n, cap=None):
    assert n >= 0
    out = [Q(1)]
    while n:
        if n & 1:
            out = mul(out, a, cap)
        n //= 2
        if n:
            a = mul(a, a, cap)
    return out


def compose(a, h):
    out = [Q(0)]
    for x in reversed(a):
        out = add(mul(out, h), [x])
    return out


def derivative(a):
    return trim([i * a[i] for i in range(1, len(a))])


def divide(a, b):
    a, b = trim(a), trim(b)
    assert b != [0]
    quotient = [Q(0)] * max(1, len(a) - len(b) + 1)
    while a != [0] and len(a) >= len(b):
        shift = len(a) - len(b)
        value = a[-1] / b[-1]
        quotient[shift] += value
        for i, x in enumerate(b):
            a[shift + i] -= value * x
        a = trim(a)
    return trim(quotient), a


def exact_divide(a, b):
    quotient, remainder = divide(a, b)
    if remainder != [0]:
        raise ValueError(f"Nonzero polynomial remainder: {record_poly(remainder)}")
    return quotient


def gcd_poly(a, b):
    a, b = trim(a), trim(b)
    while b != [0]:
        _, remainder = divide(a, b)
        a, b = b, remainder
    return scale(a, 1 / a[-1])


def degree(a):
    assert trim(a) != [0]
    return len(trim(a)) - 1


def coefficient(a, i):
    return a[i] if i < len(a) else Q(0)


def inverse_series(a, maximum):
    assert a[0] != 0
    out = [1 / a[0]]
    for n in range(1, maximum + 1):
        out.append(-sum(coefficient(a, i) * out[n - i] for i in range(1, n + 1)) / a[0])
    return out


def normalized_sequence(A, P, maximum):
    assert A[0] == P[0] == 1
    return [inverse_series(mul(power(A, m, m), P, m), m)[m]
            for m in range(maximum + 1)]


def record_poly(a):
    return [str(x) for x in trim(a)]


def seed(values, multiplicities, aggregate_exponents):
    assert len(values) == len(multiplicities) == len(aggregate_exponents)
    assert len(set(values)) == len(values) and all(values)
    assert all(p >= 1 for p in multiplicities)
    assert all(r >= 0 for r in aggregate_exponents)
    P, A = [Q(1)], [Q(1)]
    for c, p, r in zip(values, multiplicities, aggregate_exponents):
        P = mul(P, power([1, -Q(c)], p))
        A = mul(A, power([1, -Q(c)], r))
    return A, P


def pullback(A0, P0, h):
    h = trim(h)
    if h[0] != 0 or coefficient(h, 1) != 1:
        raise ValueError("The cover must have h(0)=0 and h'(0)=1.")
    hp = derivative(h)
    Q0h = mul(h, compose(P0, h))
    Qnew = exact_divide(Q0h, hp)
    P = exact_divide(Qnew, [0, 1])
    A = exact_divide(mul(h, compose(A0, h)), [0, 1])
    assert P[0] == A[0] == 1
    assert mul(P, mul([0, 1], hp)) == Q0h
    assert mul(A, compose(P0, h)) == mul(P, mul(compose(A0, h), hp))
    d = degree(h)
    assert degree(P) == d * degree(P0)
    assert degree(A) + 1 == d * (degree(A0) + 1)
    assert gcd(degree(P), degree(A) + 1) % d == 0
    old_distinct = degree(P0) - degree(gcd_poly(P0, derivative(P0)))
    new_distinct = degree(P) - degree(gcd_poly(P, derivative(P)))
    assert new_distinct == d * old_distinct
    return A, P


def compositions(m, s):
    if s == 1:
        yield (m,)
    else:
        for first in range(m + 1):
            for rest in compositions(m - first, s - 1):
                yield (first,) + rest


def qadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def qmul(a, b, discriminant):
    return (a[0] * b[0] + discriminant * a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def qpower(a, n, discriminant):
    out = (Q(1), Q(0))
    while n:
        if n & 1:
            out = qmul(out, a, discriminant)
        n //= 2
        if n:
            a = qmul(a, a, discriminant)
    return out


def actual_moment(exponents, coefficients, m, discriminant):
    out = (Q(0), Q(0))
    for k in compositions(m, len(exponents)):
        scalar = factorial(m) // prod(factorial(x) for x in k)
        scalar *= prod(factorial(m * r + j) for r, j in zip(exponents, k))
        term = (Q(scalar), Q(0))
        for c, j in zip(coefficients, k):
            term = qmul(term, qpower(c, j, discriminant), discriminant)
        out = qadd(out, term)
    return out


def direct_moment_control(name, r, c, discriminant, A, P, maximum):
    g = normalized_sequence(A, P, maximum)
    moments = []
    for m in range(maximum + 1):
        value = actual_moment(r, c, m, discriminant)
        normalizer = factorial(m) * prod(factorial(m * a) for a in r)
        assert value == (normalizer * g[m], 0)
        moments.append({"order": m, "actual_moment": [str(x) for x in value],
                        "normalizer": str(normalizer), "normalized_moment": str(g[m])})
    return {"name": name, "common_exponents": r,
            "coefficient_field": f"Q(sqrt({discriminant}))",
            "coefficients_in_basis_1_sqrt": [[str(x) for x in a] for a in c],
            "moments": moments}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=HERE.parents[2])
    parser.add_argument("--note", type=Path, default=HERE / NOTE_NAME)
    parser.add_argument("--output", type=Path, default=HERE / "sfc_ramified_pullback_verification.json")
    args = parser.parse_args()
    pins = []
    sources = [(args.note, NOTE_PIN), (args.repo_root / CONTACT_NAME, CONTACT_PIN)]
    sources += [(args.repo_root / "research/factorial/image-rigidity" / name, expected)
                for name, expected in SEED_PINS.items()]
    for path, expected in sources:
        try:
            actual = sha256(path.read_bytes()).hexdigest()
        except OSError as error:
            raise SystemExit(f"Missing or unreadable pinned source {path}: {error}") from error
        if actual != expected:
            raise SystemExit(f"Source pin mismatch for {path}: expected {expected}, got {actual}")
        try:
            display_path = path.resolve().relative_to(args.repo_root.resolve())
        except ValueError:
            display_path = path.resolve()
        pins.append({"path": str(display_path), "sha256": actual})
    receipt = {
        "status": "PASS",
        "scope": "Exact finite certificates for written all-parameter constrained-cover, cyclic-family and descent proofs, with exact failure controls; no arbitrary-coefficient row-translation or general sharp-contact theorem.",
        "sources": pins,
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "polynomial_coefficient_order": "ascending",
    }
    maximum = 12
    cases = []

    def retain(name, A0, P0, h, expected=None):
        A, P = pullback(A0, P0, h)
        before = normalized_sequence(A0, P0, maximum)
        after = normalized_sequence(A, P, maximum)
        assert before == after
        if expected is not None:
            assert after == expected
        cases.append({"name": name, "A0": record_poly(A0), "P0": record_poly(P0),
                      "h": record_poly(h), "A": record_poly(A), "P": record_poly(P),
                      "normalized_moments_orders_0_through_12": [str(x) for x in after],
                      "actual_support": degree(P),
                      "common_exponent_total": degree(A),
                      "support_degree_gcd": gcd(degree(P), degree(A) + 1),
                      "distinct_coefficients": degree(P) - degree(gcd_poly(P, derivative(P)))})
        return A, P

    # The one-term seed also checks d=1, where h is exactly the identity.
    c = Q(3)
    for r in (0, 1, 3):
        P0 = [Q(1), -c]
        A0 = power(P0, r)
        expected = [Q(comb((r + 1) * m, m)) * c**m for m in range(maximum + 1)]
        for d in (1, 2, 4):
            h = scale(add([1], scale(power([1, -c / d], d), -1)), 1 / c)
            retain(f"one_term_seed_r{r}_degree{d}", A0, P0, h, expected)

    # The all-parameter cyclic formula is proved in the note. These finite
    # witnesses retain every polynomial, not merely a pass flag.
    for q in (2, 3, 4, 5):
        P0 = [Q(1)] + [Q(0)] * (q - 1) + [Q(-1)]
        for b in (0, 1, 2):
            A0 = power(P0, b)
            expected = [Q(comb((b * q + 1) * (m // q), m // q)) if m % q == 0 else Q(0)
                        for m in range(maximum + 1)]
            for d in (2, 3):
                h = add([1], scale(power([1, -Q(1, d)], d), -1))
                A, P = retain(f"cyclic_q{q}_b{b}_degree{d}", A0, P0, h, expected)
                # The unique ramified active fiber is z=d. Its P multiplicity
                # is one and its A multiplicity is db, including b=0.
                linear = [Q(1), -Q(1, d)]
                assert divide(P, linear)[1] == [0]
                assert divide(P, power(linear, 2))[1] != [0]
                A_quotient = exact_divide(A, power(linear, d * b))
                assert divide(A_quotient, linear)[1] != [0]

    # Integer quadratic covers, with cyclic seed rescaled to their critical value.
    for q in (2, 3):
        P0 = [Q(1)] + [Q(0)] * (q - 1) + [Q(-(4**q))]
        for b in (1, 2):
            expected = [Q(4**m * comb((b * q + 1) * (m // q), m // q))
                        if m % q == 0 else Q(0) for m in range(maximum + 1)]
            retain(f"integer_quadratic_q{q}_b{b}", power(P0, b), P0, [0, 1, -1], expected)

    # Complete integer-polynomial factorization with two critical values.
    h3 = [Q(0), Q(1), Q(-3), Q(2)]
    D = [Q(1), Q(-6), Q(6)]
    C = [Q(1), Q(12), Q(-12)]
    assert derivative(h3) == D
    factor_left = add([1], scale(power(h3, 2), -108))
    factor_right = mul(power(D, 2), C)
    assert factor_left == factor_right
    P0_integer = [Q(1), Q(0), Q(-108)]
    integer_cubic_cases = {}
    for b in (0, 1, 2, 3):
        expected = [Q(108 ** (m // 2) * comb((2 * b + 1) * (m // 2), m // 2))
                    if m % 2 == 0 else Q(0) for m in range(maximum + 1)]
        A, P = retain(f"integer_cubic_b{b}", power(P0_integer, b), P0_integer, h3, expected)
        assert P == mul(mul([1, -1], [1, -2]), mul(D, C))
        assert A == mul(mul([1, -1], [1, -2]), mul(power(D, 2 * b), power(C, b)))
        integer_cubic_cases[b] = (A, P)
    receipt["integer_critical_value_certificate"] = {
        "h": record_poly(h3), "h_prime": record_poly(D), "C": record_poly(C),
        "one_minus_108_h_squared": record_poly(factor_left),
        "h_prime_squared_times_C": record_poly(factor_right),
    }

    # An iterated cover has a critical value over zero as well as over one.
    h1 = [Q(0), Q(1), -Q(1, 4)]
    h2 = [Q(0), Q(1), -Q(1, 16)]
    cyclic2 = [Q(1), Q(0), Q(-1)]
    A_composite, P_composite = retain(
        "iterated_degree_four_with_critical_value_zero", cyclic2, cyclic2, compose(h1, h2))
    A_stage1, P_stage1 = pullback(cyclic2, cyclic2, h1)
    A_stage2, P_stage2 = pullback(A_stage1, P_stage1, h2)
    assert (A_stage2, P_stage2) == (A_composite, P_composite)
    receipt["composition_control"] = {
        "first_cover": record_poly(h1), "second_cover": record_poly(h2),
        "composite": record_poly(compose(h1, h2)),
        "two_stage_equals_one_stage": True,
    }

    # This exact descent explains the permanent t=2 residue-rank factor.
    h_t2 = [Q(0), Q(1), Q(-1)]
    seed_t2 = [Q(1), Q(-4)]
    expected_t2 = [Q(4**m * comb(2*m, m)) for m in range(maximum + 1)]
    A_t2, P_t2 = retain("retained_t2_observability_descent", seed_t2, seed_t2, h_t2, expected_t2)
    assert A_t2 == mul([1, -1], power([1, -2], 2))
    assert P_t2 == mul([1, -1], [1, -2])
    a_t2, Q_t2 = mul([0, 1], A_t2), mul([0, 1], P_t2)
    residue_controls = []
    for m in range(7):
        denominator = mul(power(a_t2, m), Q_t2)
        residues = []
        for beta in (Q(0), Q(1), Q(1, 2)):
            shifted = compose(denominator, [beta, 1])
            order = next(i for i, value in enumerate(shifted) if value)
            residue = inverse_series(shifted[order:], order - 1)[order - 1]
            residues.append(residue)
        assert residues == [expected_t2[m], expected_t2[m], -2*expected_t2[m]]
        residue_controls.append({"order": m, "residues_at_0_1_half": list(map(str, residues))})
    receipt["retained_t2_residue_controls"] = residue_controls

    # Zero common exponents are retained as active points of P.
    h_boundary = [Q(0), Q(1), -Q(1, 2)]
    A0_zero, P0_zero = seed([Q(2), Q(-3)], [1, 1], [0, 2])
    A_zero, P_zero = retain("zero_common_exponent", A0_zero, P0_zero, h_boundary)
    fhalf, fone, fquad = [1, -Q(1, 2)], [1, -1], [1, 3, -Q(3, 2)]
    assert P_zero == mul(mul(fhalf, fone), fquad)
    assert A_zero == mul(fhalf, power(fquad, 2))
    assert normalized_sequence(A_zero, exact_divide(P_zero, fone), 1)[1] != normalized_sequence(A_zero, P_zero, 1)[1]

    # Repeated coefficient fibers have multiplicity e(p-1)+1, not one.
    A0_repeat, P0_repeat = seed([Q(2), Q(-3)], [2, 1], [1, 3])
    A_repeat, P_repeat = retain("repeated_coefficient_fiber", A0_repeat, P0_repeat, h_boundary)
    assert P_repeat == mul(mul(fhalf, power(fone, 3)), fquad)
    assert A_repeat == mul(mul(fhalf, power(fone, 2)), power(fquad, 3))
    bad_simple = mul(mul(fhalf, fone), fquad)
    assert normalized_sequence(A_repeat, bad_simple, 1)[1] != normalized_sequence(A_repeat, P_repeat, 1)[1]

    receipt["polynomial_pullback_cases"] = cases
    receipt["coefficient_comparison_maximum_order"] = maximum

    # Independent actual factorial expansions in explicitly represented fields.
    actual_controls = []
    A_four, P_four = pullback(cyclic2, cyclic2, h1)
    c_four = [(Q(1, 4), Q(0)), (Q(1, 2), Q(0)),
              (-Q(1, 2), Q(1, 2)), (-Q(1, 2), -Q(1, 2))]
    control = direct_moment_control("four_term_cancellation", [1, 2, 1, 1], c_four, 2, A_four, P_four, 8)
    assert control["moments"][1]["actual_moment"] == ["0", "0"]
    assert control["moments"][2]["actual_moment"] == ["1152", "0"]
    actual_controls.append(control)

    # The q=3 cover has four nonreal coefficient points. The exact sum-of-
    # squares identity rules out real roots of their quartic, without using
    # an ordering on the coefficient field or on the complex moments.
    cyclic3 = [Q(1), Q(0), Q(0), Q(-1)]
    A_complex, P_complex = pullback(cyclic3, cyclic3, h1)
    quartic = add(add([1], h1), power(h1, 2))
    positive_identity = add(power(add(scale(h1, 2), [1]), 2), [3])
    assert scale(quartic, 4) == positive_identity
    assert P_complex == mul(mul([1, -Q(1, 4)], [1, -Q(1, 2)]), quartic)
    assert A_complex == mul(P_complex, [1, -Q(1, 2)])
    complex_g = normalized_sequence(A_complex, P_complex, 3)
    assert complex_g == [1, 0, 0, 4]
    assert factorial(3)**6 * factorial(6) * complex_g[3] == 134369280
    receipt["nonreal_six_term_cancellation_control"] = {
        "h": record_poly(h1), "quartic_with_four_nonreal_roots": record_poly(quartic),
        "four_times_quartic_equals_two_h_plus_one_squared_plus_three": record_poly(positive_identity),
        "P": record_poly(P_complex), "A": record_poly(A_complex),
        "normalized_moments_0_through_3": [str(x) for x in complex_g],
        "ordinary_moments_1_through_3": ["0", "0", "134369280"],
        "scope": "Exact polynomial factorization and residue identity; the four nonreal coefficients are reciprocals of the displayed quartic roots.",
    }

    c_six = [(Q(1), Q(0)), (Q(2), Q(0)), (Q(3), Q(1)), (Q(3), Q(-1)),
             (Q(-6), Q(4)), (Q(-6), Q(-4))]
    for b in (1, 2):
        A, P = integer_cubic_cases[b]
        actual_controls.append(direct_moment_control(
            f"integer_cubic_six_terms_b{b}", [1, 1, 2 * b, 2 * b, b, b], c_six, 3, A, P, 6))

    c_boundary = [(Q(1, 2), Q(0)), (Q(1), Q(0)),
                  (-Q(3, 2), Q(1, 2)), (-Q(3, 2), -Q(1, 2))]
    actual_controls.append(direct_moment_control(
        "zero_exponent_actual_phase", [1, 0, 2, 2], c_boundary, 15, A_zero, P_zero, 6))
    c_repeat = [(Q(1, 2), Q(0)), (Q(1), Q(0)), (Q(1), Q(0)), (Q(1), Q(0)),
                (-Q(3, 2), Q(1, 2)), (-Q(3, 2), -Q(1, 2))]
    for allocation in ([0, 0, 2], [1, 1, 0]):
        actual_controls.append(direct_moment_control(
            f"repeated_fiber_allocation_{allocation}", [1] + allocation + [3, 3],
            c_repeat, 15, A_repeat, P_repeat, 6))
    actual_controls.append(direct_moment_control(
        "repeated_seed_actual_phase", [1, 0, 3], [(Q(2), Q(0)), (Q(2), Q(0)), (Q(-3), Q(0))],
        15, A0_repeat, P0_repeat, 6))
    receipt["actual_factorial_expansion_controls"] = actual_controls

    # Negative controls must reject before a purported pullback is retained.
    bad_h = [Q(0), Q(1), Q(1)]
    bad_P0 = [Q(1), Q(-1)]
    numerator = mul(bad_h, compose(bad_P0, bad_h))
    _, remainder = divide(numerator, derivative(bad_h))
    assert remainder == [-Q(5, 16)]
    rejected = []
    for name, A0, P0, h in (
            ("critical_value_outside_seed_fibers", bad_P0, bad_P0, bad_h),
            ("unnormalized_linear_coefficient", cyclic2, cyclic2, [0, 2]),
            ("nonzero_constant_term", cyclic2, cyclic2, [1, 1])):
        try:
            pullback(A0, P0, h)
        except ValueError as error:
            rejected.append({"name": name, "reason": str(error)})
        else:
            raise AssertionError(f"Invalid pullback accepted: {name}")
    receipt["negative_controls"] = {
        "inadmissible_cover_remainder": record_poly(remainder),
        "rejected_inputs": rejected,
        "zero_exponent_fiber_cannot_be_deleted": True,
        "repeated_fiber_multiplicity_cannot_be_replaced_by_one": True,
    }

    # Complete exact certificate refuting the derivative-quotient-only bound.
    H_contact = [Q(0), Q(1), -Q(1, 2)]
    D_contact = [Q(1), Q(3, 5), Q(6, 25)]
    a_contact = [Q(0), Q(1), -Q(1, 5), -Q(3, 25), -Q(3, 50)]
    b_contact = [Q(0), Q(1), -Q(3, 10)]
    assert derivative(a_contact) == mul(derivative(H_contact), D_contact)
    error = add(H_contact, scale(compose(b_contact, a_contact), -1))
    order = next(i for i, value in enumerate(error) if value)
    assert order == 5 and error[5] == -Q(27, 1250)
    assert order > degree(b_contact) + degree(D_contact)
    critical_value = compose(a_contact, [Q(1)])[0]
    assert compose(derivative(H_contact), [Q(1)]) == [0]
    assert critical_value == Q(31, 50)
    receipt["weakened_contact_counterexample"] = {
        "H": record_poly(H_contact), "D": record_poly(D_contact),
        "a": record_poly(a_contact), "b": record_poly(b_contact),
        "H_minus_b_of_a": record_poly(error), "contact_order": order,
        "degree_bound_refuted": degree(b_contact) + degree(D_contact),
        "a_at_shared_critical_point": str(critical_value),
        "scope": "Refutes a derivative-quotient-only contact assertion; violates the required zero-fiber condition and is not an SFC counterexample.",
    }
    degree_controls = 0
    for S in range(1, 9):
        for R in range(S, S + 12):
            assert gcd(R + 1, R - S + 1) == gcd(R + 1, S)
            degree_controls += 1
    assert gcd(3, 5) == gcd(5, 2) == 1
    grouped_degree_controls = []
    for p, r, d in (([2, 1], [1, 3], 2), ([1, 1], [0, 2], 2),
                    ([3, 2], [0, 0], 3), ([2, 4], [3, 1], 5)):
        S0, R0 = sum(p), sum(r)
        S, R = d * S0, d * (R0 + 1) - 1
        assert d * (S0 + 1) - (d - 1) == S + 1
        assert gcd(S, R + 1) % d == 0
        assert Q(R + 1, S) == Q(R0 + 1, S0)
        grouped_degree_controls.append({"seed_multiplicities": p, "seed_totals": r,
                                        "cover_degree": d, "target_S": S, "target_R": R})
    receipt["degree_obstruction_controls"] = {
        "positive_gcd_identities": degree_controls,
        "grouped_degree_arithmetic": grouped_degree_controls,
        "primitive_112_gcd": 1,
        "scope": "Finite arithmetic controls only; the degree obstruction and descent converse have separate general written proofs.",
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    actual_orders = sum(len(control["moments"]) for control in actual_controls)
    print(f"PASS: {len(cases)} exact polynomial pullbacks through order {maximum}; {actual_orders} direct quadratic-field factorial moments; multiplicity, descent and invalid-cover controls.")
    print(f"Receipt: {args.output}")


if __name__ == "__main__":
    main()
