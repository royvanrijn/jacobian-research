#!/usr/bin/env python3
"""Exact controls for uniform inverse-power windows and rational block closure.

All arithmetic is standard-library rational arithmetic.  The universal
statements are proved in the pinned note, not inferred from these inputs.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import combinations, permutations
import json
from math import comb, factorial, gcd as integer_gcd, prod
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
NOTE = HERE / "SFC_ALGEBRAIC_UNIFORM_WINDOWS.md"
NOTE_SHA = "5066e9fef20d387d80bf865cb07f19563bb5907ace2c8e5c0032b29028f3b2f4"
HELPER = HERE / "verify_sfc_inverse_power_contact.py"
HELPER_SHA = "cf9b7d8f7df5109c20ad623cfbddccd2d2902a57a3b4253ac988be8e31bca02d"


def need(value, message):
    if not value:
        raise ValueError(message)


need(sha256(HELPER.read_bytes()).hexdigest() == HELPER_SHA, "helper hash mismatch")
spec = importlib.util.spec_from_file_location("inverse_power_arithmetic", HELPER)
u = importlib.util.module_from_spec(spec)
spec.loader.exec_module(u)
trim, add, scale, mul, power = u.trim, u.add, u.scale, u.mul, u.power
der, compose, evaluate = u.derivative, u.compose, u.evaluate


def sub(a, b):
    return add(a, scale(b, -1))


def deg(a):
    return -1 if trim(a) == [0] else len(trim(a)) - 1


def division(a, b):
    a, b = list(map(F, trim(a))), list(map(F, trim(b)))
    need(b != [0], "zero polynomial denominator")
    q = [F(0)] * max(1, len(a) - len(b) + 1)
    while a != [0] and len(a) >= len(b):
        k, value = len(a) - len(b), a[-1] / b[-1]
        q[k] += value
        for j, c in enumerate(b):
            a[j + k] -= value * c
        a = trim(a)
    return trim(q), a


def exact(a, b):
    q, r = division(a, b)
    need(r == [0], "nonexact polynomial division")
    return q


def gcd(a, b):
    while trim(b) != [0]:
        a, b = b, division(a, b)[1]
    return scale(a, F(1) / a[-1])


def interpolate(values):
    out = [F(0)]
    for i, value in enumerate(values):
        basis, denominator = [F(1)], 1
        for j in range(len(values)):
            if j != i:
                basis = mul(basis, [-j, 1])
                denominator *= i - j
        out = add(out, scale(basis, F(value, denominator)))
    return out


def critical_polynomial(a):
    # Res_z(A(z)-t,A'(z)) has t-degree at most D-1.
    values = []
    for t in range(deg(a)):
        aa = a[:]
        aa[0] -= t
        values.append(u.resultant(aa, der(a)))
    result = interpolate(values)
    need(deg(result) == deg(a) - 1, "critical resultant degree")
    radical = exact(result, gcd(result, der(result)))
    return scale(radical, F(1) / radical[-1])


def basis_table(a, maximum):
    """Reduce successive powers by multiplication by z, independent of long division."""
    d = deg(a)
    table = [[[F(1)]] + [[F(0)]] * (d - 1)]
    for _ in range(maximum):
        previous = table[-1]
        last = previous[-1]
        row = [[F(0)]] + [x[:] for x in previous[:-1]]
        for j in range(d):
            row[j] = sub(row[j], scale(last, F(a[j], a[-1])))
        row[0] = add(row[0], scale(mul(last, [0, 1]), F(1, a[-1])))
        table.append(row)
    return table


def reduce_polynomial(p, table):
    row = [[F(0)] for _ in table[0]]
    for k, coefficient in enumerate(p):
        row = [add(x, scale(y, coefficient)) for x, y in zip(row, table[k])]
    return row


def reconstruct(row, a):
    out = [F(0)]
    for j, value in enumerate(row):
        out = add(out, [F(0)] * j + compose(value, a))
    return out


def newton_traces(a, maximum):
    d = deg(a)
    coefficients = [[F(a[d - j], a[-1])] for j in range(1, d + 1)]
    coefficients[-1] = [F(a[0], a[-1]), F(-1, a[-1])]
    traces = [[F(d)]]
    for k in range(1, maximum + 1):
        value = [F(0)]
        for j in range(1, min(k - 1, d) + 1):
            value = sub(value, mul(coefficients[j - 1], traces[k - j]))
        if k <= d:
            value = sub(value, scale(coefficients[k - 1], k))
        traces.append(value)
    return traces


def poly_det(matrix):
    """Leibniz determinant; small minors only, no polynomial Bareiss reuse."""
    n = len(matrix)
    out = [F(0)]
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        value = [F((-1) ** inversions)]
        for i in range(n):
            value = mul(value, matrix[i][p[i]])
        out = add(out, value)
    return out


def full_relation(vectors, minors):
    return [sum_polynomials(mul(m, v[row]) for m, v in zip(minors, vectors))
            for row in range(len(vectors[0]))]


def sum_polynomials(values):
    out = [F(0)]
    for value in values:
        out = add(out, value)
    return out


def operator_case(a, b):
    a = list(map(F, a))
    d, m = deg(a), b // deg(a)
    need(d >= 2 and a[0] == 0 and a[1] == 1, "normalized nonlinear polynomial required")
    sigma = critical_polynomial(a)
    c = deg(sigma)
    ja = exact(compose(sigma, a), der(a))
    need(deg(ja) == (c - 1) * d + 1, "connection quotient degree")
    table = basis_table(a, max(b, deg(ja) + d - 2))
    columns = []
    for j in range(d):
        numerator = [F(0)] if j == 0 else [F(0)] * (j - 1) + scale(ja, j)
        row = reduce_polynomial(numerator, table)
        need(reconstruct(row, a) == numerator, "connection basis reconstruction")
        need(all(deg(x) <= c - 1 for x in row), "connection degree")
        columns.append(row)
    matrix = [[columns[j][i] for j in range(d)] for i in range(d)]

    traces = newton_traces(a, max(b, d - 1))
    raw = table[b]
    trace_from_basis = sum_polynomials(mul(raw[j], traces[j]) for j in range(d))
    need(trace_from_basis == traces[b], "independent Newton trace")
    tau = scale(traces[b], F(1, d))
    need(deg(tau) <= m, "trace degree")
    centered = [x[:] for x in raw]
    centered[0] = sub(centered[0], tau)
    need(any(x != [0] for x in centered), "centered observation is zero")
    need(reconstruct(centered, a) == sub([0] * b + [1], compose(tau, a)),
         "centered polynomial reconstruction")
    need(all(deg(x) <= m for x in centered), "centered coordinate degree")
    need(sum_polynomials(mul(centered[j], traces[j]) for j in range(d)) == [0],
         "zero trace")

    vectors = [centered]
    certificate = None
    for k in range(d - 1):
        vector = vectors[-1]
        cv = [sum_polynomials(mul(matrix[i][j], vector[j]) for j in range(d))
              for i in range(d)]
        following = [add(sub(mul(sigma, der(x)), scale(mul(der(sigma), x), k)), y)
                     for x, y in zip(vector, cv)]
        need(all(deg(x) <= m + (k + 1) * (c - 1) for x in following),
             "derivative numerator degree")
        vectors.append(following)
        q = k + 1
        selected = None
        for rows in combinations(range(d), q):
            if poly_det([[vectors[j][i] for j in range(q)] for i in rows]) != [0]:
                selected = rows
                break
        need(selected is not None, "first derivative block not independent")
        minors = [scale(poly_det([[vectors[j][i] for j in range(q + 1) if j != omitted]
                                  for i in selected]), (-1) ** omitted)
                  for omitted in range(q + 1)]
        need(minors[-1] != [0], "zero leading cofactor")
        if all(x == [0] for x in full_relation(vectors, minors)):
            certificate = selected, minors
            break
    need(certificate is not None, "trace-zero rank bound failed")
    rows, minors = certificate
    q = len(vectors) - 1
    epsilon = int(sigma[0] != 0)
    M = q * m + (c - 1) * q * (q + 1) // 2
    bound = M + epsilon * q
    operator = [mul(x, power(sigma, k)) for k, x in enumerate(minors)]
    for k, p in enumerate(operator):
        need(deg(minors[k]) <= q * m + (c - 1) * (q * (q + 1) // 2 - k),
             "minor degree")
        need(deg(p) - k <= M, "operator degree")
        if not epsilon and p != [0]:
            need(u.valuation(p) >= k, "zero-fiber operator valuation")
    altered = [x[:] for x in minors]
    altered[0] = add(altered[0], [1])
    need(any(x != [0] for x in full_relation(vectors, altered)),
         "altered cofactor was not rejected")

    maximum = b + 2 + bound
    inverse_power = power(u.inverse(a, maximum), b, maximum)
    for N in (b, b + 1):
        approximation = inverse_power[:N + 1]
        shifted = sub(approximation, tau)
        need(shifted != [0], "actual truncation equals its lower-degree trace")
        value = [F(0)]
        deriv = shifted
        for p in operator:
            value = add(value, mul(p, deriv))
            deriv = der(deriv)
        need(value != [0], "unexpected polynomial solution")
        need(deg(value) <= N + M, "polynomial image degree")
        contact = u.valuation(sub(inverse_power, approximation))
        need(contact is not None and contact <= N + bound, "finite contact control")
    exceptional = sub(inverse_power, tau)
    need(u.valuation(exceptional) is not None and u.valuation(exceptional) <= b,
         "Q=tau exceptional contact")
    return {"a": a, "b": b, "D": d, "m": m, "sigma": sigma, "c": c,
            "epsilon": epsilon, "tau": tau, "connection": matrix,
            "derivative_numerators": vectors, "independent_rows": rows,
            "cofactors": minors, "operator": operator,
            "minimal_order_q": q, "M": M, "window_constant": bound}


def ordinary_normalized_moment(a, b, multiplicities, coefficients, n):
    exponents, values = [], []
    for v, c in zip(multiplicities, coefficients):
        exponents += [a * v] + [0] * (b * v - 1)
        values += [c] * (b * v)
    product = [F(1)]
    for r, c in zip(exponents, values):
        local = [F(factorial(n * r + k), factorial(k)) * F(c) ** k
                 for k in range(n + 1)]
        product = mul(product, local, n)
    raw = factorial(n) * u.coefficient(product, n)
    denominator = factorial(n) * prod(factorial(n * r) for r in exponents)
    return F(raw, denominator)


def block_case(a, b, multiplicities, coefficients, refined=False):
    h, nu = len(multiplicities), sum(multiplicities)
    d, support = 1 + a * nu, b * nu
    left = 2 * (b - b // d)
    right = (a * h - 1) * d + 2 if refined else a * h * d
    need(left >= right, "sufficient block condition failed")
    r = [F(1)]
    for v, c in zip(multiplicities, coefficients):
        r = mul(r, power([1, -c], v))
    aa = [F(0)] * (1 + a * deg(r) + 1)
    for k, value in enumerate(r):
        aa[1 + a * k] = value
    maximum_moment = 4
    maximum_degree = a * maximum_moment + b
    series = power(u.inverse(aa, maximum_degree), b, maximum_degree)
    moments = []
    for n in range(maximum_moment + 1):
        grouped = [F(1)]
        for v, c in zip(multiplicities, coefficients):
            shape = (a * n + b) * v
            grouped = mul(grouped, [F(comb(shape + k - 1, k)) * F(c) ** k
                                    for k in range(n + 1)], n)
        g = u.coefficient(grouped, n)
        need(g == ordinary_normalized_moment(a, b, multiplicities, coefficients, n),
             "ordinary factorial/grouped moment identity")
        need(u.coefficient(series, a * n + b) == F(b, a * n + b) * g,
             "inverse-power Lagrange identity")
        moments.append(g)
    need(all(value == 0 for i, value in enumerate(series) if (i - b) % a),
         "inverse-power sparse congruence")
    coefficient_constant = (d - 1) * (b // d)
    if refined:
        coefficient_constant += (a * h - 1) * d * (d - 1) // 2 + d - 1
    else:
        coefficient_constant += a * h * d * (d - 1) // 2
    need(coefficient_constant <= a * support, "support-size absorption")
    for start in (1, 2, 11):
        N = a * (start - 1) + b
        need(N >= b and N + a * (support + 1) > N + coefficient_constant,
             "sparse contact conversion")
    return {"a": a, "b": b, "multiplicities": multiplicities,
            "coefficients": coefficients, "D": d, "support": support,
            "refined_zero_free_test": refined, "condition_left": left,
            "condition_right": right, "coefficient_window_bound": coefficient_constant,
            "moment_window_bound": coefficient_constant // a,
            "direct_normalized_moments": moments}


def encode(value):
    if isinstance(value, F):
        return {"numerator": value.numerator, "denominator": value.denominator}
    raise TypeError(type(value).__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=NOTE)
    parser.add_argument("--output", type=Path,
                        default=HERE / "sfc_rational_block_sfc_verification.json")
    args = parser.parse_args()
    need(sys.flags.optimize == 0, "optimized execution is not an approved replay")
    need(sha256(args.note.read_bytes()).hexdigest() == NOTE_SHA, "note hash mismatch")
    cases = [([0, 1, -1], 2), ([0, 1, -1], 4),
             ([0, 1, 0, -1], 1), ([0, 1, 0, -1], 2),
             ([0, 1, -2, 1], 2), ([0, 1, -3, 2], 4),
             ([0, 1, -1, -1, 1], 5)]
    records = [operator_case(a, b) for a, b in cases]
    # Outside-domain q=0 control: centering the constant observation is zero.
    try:
        operator_case([0, 1, -1], 0)
    except ValueError as exc:
        need("centered observation is zero" in str(exc), "wrong zero-center rejection")
    else:
        raise ValueError("zero-centered observation accepted")
    # An actual failure if the N>=b condition is omitted without a replacement.
    small = records[1]
    need(small["m"] == 2 and small["window_constant"] == 3,
         "small-degree countercontrol parameters")
    need(4 > 0 + small["window_constant"], "small-degree countercontrol")

    blocks = [block_case(1, 4, (1, 1), (1, 2)),
              block_case(1, 4, (1, 1), (1, -1), refined=True),
              block_case(1, 5, (2, 1), (1, -1)),
              block_case(1, 6, (1, 1, 1), (1, 2, 3), refined=True),
              block_case(2, 12, (1, 1), (1, -2))]
    for h in range(2, 9):
        d, b = h + 1, (h + 1) ** 2 + 1
        need(2 * (b - b // d) >= h * d, "primitive family sufficient inequality")
        need(integer_gcd(h * b, d) == 1, "primitive family degree gcd")
    receipt = {"status": "PASS", "note_sha256": NOTE_SHA,
               "helper_sha256": HELPER_SHA,
               "scope": "Finite exact certificates; the note proves the uniform theorems",
               "operator_case_count": len(records), "operator_cases": records,
               "block_case_count": len(blocks), "block_cases": blocks,
               "actual_factorial_moment_checks": 5 * len(blocks),
               "controls": ["nonzero first-block minors and full coordinate relation",
                            "altered cofactors rejected", "zero-centered observation rejected",
                            "N<m countercontrol requires the separate elementary bound",
                            "Q=tau exception checked", "both zero-critical-value cases checked"],
               "general_sharp_inverse_power_claimed": False,
               "unrestricted_SFC_claimed": False}
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True, default=encode) + "\n")
    print(f"PASS: {len(records)} inverse-power operators; {len(blocks)} block families; "
          f"{5 * len(blocks)} actual factorial identities; exact controls")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as exc:
        raise SystemExit("REJECT: " + str(exc)) from exc
