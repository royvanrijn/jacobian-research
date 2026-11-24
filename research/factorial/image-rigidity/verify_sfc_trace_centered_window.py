#!/usr/bin/env python3
"""Exact finite controls for trace-centered uniform algebraic windows.

The written proof supplies universal quantifiers. These four cases authenticate
their critical-value radicals, actual polynomial connections, minimal-order
cofactor equations, degree bounds and coefficient contacts. Standard library.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import permutations
import json
from math import comb, lcm
from pathlib import Path

if not __debug__:
    raise SystemExit("REJECT: run the verification without -O")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / "research/factorial/image-rigidity"
NOTE = HERE / "SFC_ALGEBRAIC_UNIFORM_WINDOWS.md"
NOTE_SHA = "5066e9fef20d387d80bf865cb07f19563bb5907ace2c8e5c0032b29028f3b2f4"
PINS = {
    "SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md":
        "409ca1ed8a1ef3218fb07673d14830d04cf10562f5714c88bd8de8e9a681b417",
    "verify_sfc_monomial_linear_contact.py":
        "c80ea0aa4ade828a89d18e028072c65fc0525c6700afc661183b5cb14c2e1055",
    "POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md":
        "40278040a10d4813c473e8256b868e4fb6d42f20912972b801136398436609d6",
}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def load_helper():
    for name, expected in PINS.items():
        path = SOURCE / name
        need(path.is_file(), "Missing pinned source: " + name)
        need(sha256(path.read_bytes()).hexdigest() == expected,
             "Source pin mismatch: " + name)
    spec = importlib.util.spec_from_file_location(
        "pinned_contact_polynomial_helper", SOURCE / "verify_sfc_monomial_linear_contact.py")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


def division(a, b, h):
    a, b = h.trim(list(map(F, a))), h.trim(list(map(F, b)))
    need(b != [0], "Polynomial division by zero")
    out = [F(0)] * max(1, len(a) - len(b) + 1)
    while a != [0] and len(a) >= len(b):
        k, c = len(a) - len(b), a[-1] / b[-1]
        out[k] += c
        for j, value in enumerate(b):
            a[j + k] -= c * value
        a = h.trim(a)
    return h.trim(out), a


def polynomial_gcd(a, b, h):
    while b != [0]:
        _, remainder = division(a, b, h)
        a, b = b, remainder
    return h.scale(a, 1 / F(a[-1]))


def degree(a):
    return -1 if a == [0] else len(a) - 1


def order(a):
    return next((j for j, c in enumerate(a) if c), None)


def value(a, x):
    out = F(0)
    for c in reversed(a):
        out = out * x + c
    return out


def bareiss(matrix):
    a = [list(row) for row in matrix]
    n, old, sign = len(a), 1, 1
    if n == 0:
        return 1
    for k in range(n - 1):
        row = next((i for i in range(k, n) if a[i][k]), None)
        if row is None:
            return 0
        if row != k:
            a[k], a[row] = a[row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = pivot * a[i][j] - a[i][k] * a[k][j]
                need(numerator % old == 0, "Nonexact Bareiss division")
                a[i][j] = numerator // old
            a[i][k] = 0
        old = pivot
    return sign * a[-1][-1]


def determinant(matrix):
    integral, den = [], 1
    for row in matrix:
        factor = lcm(*(F(x).denominator for x in row))
        integral.append([int(F(x) * factor) for x in row])
        den *= factor
    return F(bareiss(integral), den)


def resultant(f, g):
    m, n = len(f) - 1, len(g) - 1
    matrix = [[F(0)] * (m + n) for _ in range(m + n)]
    for j in range(n):
        matrix[j][j:j + m + 1] = list(reversed(f))
    for j in range(m):
        matrix[n + j][j:j + n + 1] = list(reversed(g))
    return determinant(matrix)


def reduce_z(polynomial, a, h):
    """Return coefficients in Q[t] of f(z) modulo a(z)-t, basis 1,z,... ."""
    d = degree(a)
    rem = [[F(x)] for x in polynomial]
    rem += [[F(0)] for _ in range(max(0, d - len(rem)))]
    for k in range(len(rem) - 1, d - 1, -1):
        lead = h.scale(rem[k], 1 / F(a[-1]))
        for j in range(d):
            rem[k - d + j] = h.sub(rem[k - d + j], h.scale(lead, a[j]))
        rem[k - d] = h.add(rem[k - d], [F(0)] + lead)
        rem[k] = [F(0)]
    return [h.trim(x) for x in rem[:d]]


def polynomial_matrix_det(matrix, h):
    n, ans = len(matrix), [F(0)]
    for perm in permutations(range(n)):
        sign = (-1) ** sum(perm[i] > perm[j]
                          for i in range(n) for j in range(i + 1, n))
        term = [F(sign)]
        for i, j in enumerate(perm):
            term = h.mul(term, matrix[i][j])
        ans = h.add(ans, term)
    return ans


def matrix_vector(matrix, vector, h):
    out = []
    for row in matrix:
        total = [F(0)]
        for x, y in zip(row, vector):
            total = h.add(total, h.mul(x, y))
        out.append(total)
    return out


def trace_basis(a, h):
    d, traces = degree(a), []
    for j in range(d):
        total = [F(0)]
        for k in range(d):
            reduced = reduce_z([F(0)] * (j + k) + [F(1)], a, h)
            total = h.add(total, reduced[k])
        traces.append(total)
    return traces


def trace(vector, traces, h):
    out = [F(0)]
    for x, y in zip(vector, traces):
        out = h.add(out, h.mul(x, y))
    return out


def on_branch(vector, inverse, limit, h):
    result = [F(0)]
    for j, coefficient in enumerate(vector):
        result = h.add(result, h.mul(coefficient, h.power(inverse, j, limit), limit))
    return h.trim(result[:limit + 1])


def apply_operator(coefficients, polynomial, h, limit=None):
    out, derivative = [F(0)], list(polynomial)
    for coefficient in coefficients:
        out = h.add(out, h.mul(coefficient, derivative, limit))
        derivative = h.derivative(derivative)
    if limit is not None:
        out = out[:limit + 1]
    return h.trim(out)


def check_case(case, h):
    a, observable, sigma = (list(map(F, case[key])) for key in ("a", "H", "sigma"))
    d, hdegree, c = degree(a), degree(observable), degree(sigma)
    need(a[:2] == [0, 1] and d >= 2 and sigma[-1] == 1, "Normalization")
    need(polynomial_gcd(sigma, h.derivative(sigma), h) == [1], "Sigma not squarefree")
    aprime = h.derivative(a)
    connection_multiplier, remainder = division(h.compose(sigma, a), aprime, h)
    need(remainder == [0], "Sigma(a) is not divisible by a'")
    # This determinant has t-degree <= deg(a')=d-1. Exactly d distinct
    # evaluations prove the complete stated critical-value resultant.
    cv = h.scale(sigma, F(case["critical_resultant_scale"]))
    need(degree(cv) <= d - 1, "Critical resultant degree certificate")
    for t in range(d):
        shifted = a[:]
        shifted[0] -= t
        need(resultant(aprime, shifted) == value(cv, t),
             "Critical-value resultant polynomial identity")

    columns = [[F(0)] for _ in range(d)]
    columns[0] = [[F(0)] for _ in range(d)]
    for j in range(1, d):
        polynomial = [F(0)] * (j - 1) + h.scale(connection_multiplier, j)
        columns[j] = reduce_z(polynomial, a, h)
    connection = [[columns[j][i] for j in range(d)] for i in range(d)]
    need(all(degree(entry) <= c - 1 for row in connection for entry in row),
         "Connection polynomial degree")
    traces = trace_basis(a, h)
    raw = reduce_z(observable, a, h)
    tau = h.scale(trace(raw, traces, h), F(1, d))
    need(tau == list(map(F, case["trace_mean"])), "Actual polynomial trace")
    m = hdegree // d
    need(degree(tau) <= m and all(degree(entry) <= m for entry in raw),
         "Initial reduction/trace degree")
    centered = [x[:] for x in raw]
    centered[0] = h.sub(centered[0], tau)
    need(trace(centered, traces, h) == [0], "Initial centered trace is not zero")
    q = d - 1
    numerators = [centered]
    for k in range(q):
        old = numerators[-1]
        cpart = matrix_vector(connection, old, h)
        new = [h.add(h.sub(h.mul(sigma, h.derivative(x)),
                          h.scale(h.mul(h.derivative(sigma), x), k)), y)
               for x, y in zip(old, cpart)]
        need(all(degree(entry) <= m + (k + 1) * (c - 1) for entry in new),
             "Derivative numerator degree")
        need(trace(new, traces, h) == [0], "Derivative numerator trace")
        numerators.append(new)

    selected_rows = list(range(1, d))
    cofactors = []
    for k in range(q + 1):
        minor = [[numerators[j][i] for j in range(q + 1) if j != k]
                 for i in selected_rows]
        cofactors.append(h.scale(polynomial_matrix_det(minor, h), (-1) ** k))
    need(cofactors[-1] != [0], "Claimed minimal order is not certified")
    for i in range(d):
        total = [F(0)]
        for k in range(q + 1):
            total = h.add(total, h.mul(cofactors[k], numerators[k][i]))
        need(total == [0], "Cofactor relation fails outside selected rows")
    operator = [h.mul(uk, h.power(sigma, k)) for k, uk in enumerate(cofactors)]
    M = q * m + (c - 1) * q * (q + 1) // 2
    need(all(degree(pk) - k <= M for k, pk in enumerate(operator)),
         "Cofactor operator weighted degree bound")
    delta = int(sigma[0] != 0)
    if not delta:
        need(all(pk == [0] or order(pk) >= k for k, pk in enumerate(operator)),
             "Critical-origin operator order bound")

    precision, checked = 40, 20
    inverse = h.inverse(a, precision)
    branch = h.compose(observable, inverse, precision)
    y = h.sub(branch, tau)
    derivative = y[:]
    for k in range(q + 1):
        lhs = on_branch(numerators[k], inverse, checked, h)
        rhs = h.mul(h.power(sigma, k), derivative, checked)
        need(lhs == rhs, "Actual derivative numerators do not match the inverse branch")
        derivative = h.derivative(derivative)
    need(apply_operator(operator, y, h, checked) == [0], "Cofactor ODE on actual branch")

    contacts = []
    for N in range(max(1, m), 13):
        polynomial = h.trim(branch[:N + 1])
        difference = h.sub(branch, polynomial)
        nu = order(difference)
        need(nu is not None and nu > N, "Insufficient inverse-series precision")
        rational_difference = h.sub(polynomial, tau)
        if rational_difference == [0]:
            need(nu <= hdegree, "Exceptional Q=tau polynomial-degree bound")
            contacts.append({"N": N, "contact": nu, "case": "Q=tau"})
            continue
        forcing = apply_operator(operator, rational_difference, h)
        need(forcing != [0], "Nonzero polynomial observation solved minimal ODE")
        need(degree(forcing) <= N + M, "Forcing degree bound")
        need(order(forcing) >= nu - delta * q, "Contact/operator order loss")
        need(nu <= N + M + delta * q, "Uniform contact bound")
        contacts.append({"N": N, "contact": nu, "bound": N + M + delta * q,
                         "forcing_degree": degree(forcing),
                         "forcing_order": order(forcing)})
    special = h.sub(observable, h.compose(tau, a))
    need(special != [0] and order(special) <= hdegree, "Q=tau exceptional branch")
    need(order(y) == order(special), "Exceptional contact under substitution")

    factorial_controls = []
    if case.get("factorial_support"):
        S = case["factorial_support"]
        need(m == 0 and observable[:2] == [0, 1], "Factorial positive-frame hypotheses")
        need(M + delta * q <= S * d * (d - 1) // 2, "Uniform factorial window")
        for ell in range(13):
            actual = (ell + 1) * h.coeff(branch, ell + 1)
            if case.get("cyclic_order"):
                kappa = case["cyclic_order"]
                expected = 0 if ell % kappa else comb((kappa + 1) * (ell // kappa), ell // kappa)
            else:
                exponent = case["ordinary_common_exponent"]
                expected = comb((exponent + 1) * ell, ell)
            need(actual == expected, "Actual normalized factorial moment control")
            factorial_controls.append(str(actual))

    negative = {}
    if case["name"] == "cubic_zero_critical_value":
        bad_sigma = [F(-4, 27), F(1)]
        _, bad_remainder = division(h.compose(bad_sigma, a), aprime, h)
        need(bad_remainder != [0], "Omitting the zero critical value was not rejected")
        negative["omitted_zero_critical_value_remainder"] = bad_remainder
    if case["name"] == "quadratic_cubic_observable":
        need(degree(tau) == 1, "Nonconstant trace control absent")
        uncentered = apply_operator(operator, branch, h, checked)
        need(uncentered != [0] and uncentered == apply_operator(operator, tau, h, checked),
             "Omitting trace centering was not rejected")
        negative["uncentered_ODE_residual"] = uncentered

    return {
        "name": case["name"], "a": a, "H": observable, "sigma": sigma,
        "critical_value_resultant": cv, "critical_value_interpolation_points": list(range(d)),
        "a_degree": d, "observable_degree": hdegree, "critical_value_count": c,
        "polynomial_connection": connection, "connection_multiplier": connection_multiplier,
        "basis_traces": traces, "trace_mean": tau, "initial_degree_bound_m": m,
        "derivative_numerators": numerators, "selected_rows": selected_rows,
        "minimal_order_q": q, "nonzero_first_q_minor": h.scale(cofactors[-1], (-1) ** q),
        "cofactors": cofactors, "operator_coefficients": operator,
        "weighted_degree_bound_M": M, "origin_is_critical_value": not bool(delta),
        "contact_additive_bound": M + delta * q,
        "branch_identity_coefficient_degree_checked": checked,
        "truncation_contacts": contacts, "exceptional_Q_equals_trace_contact": order(special),
        "factorial_moment_controls": factorial_controls, "negative_controls": negative,
    }


def run_checks():
    h = load_helper()
    cases = [
        {"name": "quadratic_linear_observable", "a": [0, 1, -1], "H": [0, 1],
         "sigma": [F(-1, 4), 1], "critical_resultant_scale": -4,
         "trace_mean": [F(1, 2)], "factorial_support": 1, "ordinary_common_exponent": 1},
        {"name": "cubic_zero_critical_value", "a": [0, 1, -2, 1], "H": [0, 1, F(-1, 2)],
         "sigma": [0, F(-4, 27), 1], "critical_resultant_scale": 27,
         "trace_mean": [F(1, 3)], "factorial_support": 1, "ordinary_common_exponent": 2},
        {"name": "quadratic_cubic_observable", "a": [0, 1, -1], "H": [0, 0, 0, 1],
         "sigma": [F(-1, 4), 1], "critical_resultant_scale": -4,
         "trace_mean": [F(1, 2), F(-3, 2)]},
        {"name": "quartic_cyclic_factorial_frame", "a": [0, 1, 0, 0, -1], "H": [0, 1],
         "sigma": [F(-27, 256), 0, 0, 1], "critical_resultant_scale": -256,
         "trace_mean": [0], "factorial_support": 3, "cyclic_order": 3},
    ]
    records = [check_case(case, h) for case in cases]
    return {"status": "PASS", "source_sha256": PINS,
            "scope": "Four exact finite actual-connection/cofactor controls; "
                     "the written proof, not these controls, gives unrestricted quantifiers.",
            "cases": records, "parameter_cases": len(records),
            "complete_critical_value_resultant_determinants": sum(len(c["a"]) - 1 for c in cases),
            "minimal_order_cofactor_equations": len(records),
            "inverse_derivative_vector_identities": sum(c["minimal_order_q"] + 1 for c in records),
            "coefficient_truncation_contacts": sum(len(c["truncation_contacts"]) for c in records),
            "normalized_factorial_moment_controls": sum(len(c["factorial_moment_controls"]) for c in records),
            "negative_controls": ["zero critical value cannot be omitted",
                                  "nonconstant trace must be subtracted"]}


def qjson(value):
    if isinstance(value, F):
        return {"numerator": value.numerator, "denominator": value.denominator}
    raise TypeError(type(value).__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=NOTE)
    parser.add_argument("--output", type=Path,
                        default=HERE / "sfc_trace_centered_window_verification.json")
    args = parser.parse_args()
    need(NOTE_SHA != "PENDING_FINAL_PROOF_PIN", "Final reviewed proof pin has not been installed")
    need(sha256(args.note.read_bytes()).hexdigest() == NOTE_SHA, "Proof note pin mismatch")
    receipt = run_checks()
    receipt["note_sha256"] = NOTE_SHA
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True, default=qjson) + "\n")
    print("PASS: 4 actual connections; 4 certified minimal cofactor ODEs; "
          "11 complete critical-value determinant points; 48 contacts; 39 factorial moments.")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as exc:
        raise SystemExit("REJECT: " + str(exc)) from exc
