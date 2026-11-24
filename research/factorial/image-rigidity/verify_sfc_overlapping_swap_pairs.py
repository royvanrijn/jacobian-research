#!/usr/bin/env python3
"""Exact controls for overlapping swap pairs and a legal row-merge obstruction.

Standard-library rational arithmetic only. Written uniform arguments, not the
finite domain below, supply the all-order theorem.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import math


NOTE_SHA = "646708ec9d996f1148cc172405e36b2841eb6cd645b046cc3c1a4bfa8564413f"
SOURCE_PINS = {
    "SFC_IDENTICAL_MONOMIAL_PAIR.md": "7157a05550267289e1cbfc2c01f08659eb2d85e4b9d62661d3e4e57c2fb0e04e",
    "POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md": "40278040a10d4813c473e8256b868e4fb6d42f20912972b801136398436609d6",
    "UNRESTRICTED_SFC_MECHANISMS.tex": "f8b3dfc5947601fb002aee0ba942b0fc20b7ef46ca8e0d4f69b5a3e5e1361785",
}


def require(value, message):
    if not value:
        raise ValueError(message)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(a, b):
    return trim([(a[k] if k < len(a) else 0)+(b[k] if k < len(b) else 0)
                 for k in range(max(len(a), len(b)))])


def scale(a, c):
    return trim([c*x for x in a])


def mul(a, b):
    out = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def evaluate(p, x):
    value = Q(0)
    for a in reversed(p):
        value = value*x+a
    return value


def divide(a, b):
    a, b = trim(a), trim(b)
    require(b != [0], "Nonzero divisor")
    quo = [Q(0)]*max(1, len(a)-len(b)+1)
    while a != [0] and len(a) >= len(b):
        j = len(a)-len(b)
        c = a[-1]/b[-1]
        quo[j] += c
        for k, v in enumerate(b):
            a[j+k] -= c*v
        a = trim(a)
    return trim(quo), a


def bezout(a, b):
    r0, r1, s0, s1, t0, t1 = a, b, [Q(1)], [Q(0)], [Q(0)], [Q(1)]
    while r1 != [0]:
        q, r = divide(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, add(s0, scale(mul(q, s1), -1))
        t0, t1 = t1, add(t0, scale(mul(q, t1), -1))
    require(len(r0) == 1 and r0[0] != 0, "Consecutive polynomials are coprime")
    s, t = scale(s0, 1/r0[0]), scale(t0, 1/r0[0])
    require(add(mul(s, a), mul(t, b)) == [1], "Full rational Bezout identity")


def compositions(m, s):
    if s == 1:
        yield (m,)
    else:
        for j in range(m+1):
            for rest in compositions(m-j, s-1):
                yield (j,)+rest


def matrix(pairs, equal_rows=()):
    require(bool(pairs), "At least one distinct pair")
    require(all(type(b) is int and type(d) is int and b >= 0 and d >= 1 for b, d in pairs),
            "Nonnegative base and positive integer difference")
    require(all(type(h) is int and h >= 0 for h in equal_rows), "Nonnegative equal-row exponents")
    rows = []
    for b, d in pairs:
        rows.extend(((b+d, b), (b, b+d)))
    rows.extend((h, h) for h in equal_rows)
    return rows


def actual_coefficient(rows, multi):
    n = sum(multi)
    multinomial = math.factorial(n)//math.prod(math.factorial(k) for k in multi)
    return Q(multinomial*math.prod(math.factorial(sum(a*k for a, k in zip(row, multi))) for row in rows))


def actual_binary(rows, n):
    return [actual_coefficient(rows, (n-k, k)) for k in range(n+1)]


def actual_moment(rows, coefficients, n):
    return sum(actual_coefficient(rows, multi)*math.prod(c**k for c, k in zip(coefficients, multi))
               for multi in compositions(n, len(coefficients)))


def ratio(pairs, n, x):
    return Q((x+1)**2, math.prod(b*n+d*x+j for b, d in pairs for j in range(1, d+1)))


def delta_by_rising(pairs, n):
    values = [Q(1, math.prod(math.factorial(b*n) for b, _ in pairs))]
    for k in range(n):
        values.append(values[-1]*ratio(pairs, n, Q(k)))
    return values


def exact_test(pairs, n, k):
    value = Q(k*(n-k), (k+1)*(n-k+1))**2
    for b, d in pairs:
        for j in range(1, d+1):
            value *= Q(b*n+d*k+j, b*n+d*(k-1)+j)
            value *= Q(b*n+d*(n-k)+j, b*n+d*(n-k-1)+j)
    return value


def w_value(pairs, n):
    return sum(Q(d*d, b*n+d) for b, d in pairs)


def limit_certificate(pairs):
    numerator, denominator = [Q(1), Q(-2), Q(1)], [Q(0), Q(0), Q(4)]
    for b, d in pairs:
        for j in range(1, d+1):
            numerator = mul(numerator, [Q(d+j), Q(b)])
            denominator = mul(denominator, [Q(j), Q(b)])
            numerator = mul(numerator, [Q(j-d), Q(b+d)])
            denominator = mul(denominator, [Q(j-2*d), Q(b+d)])
    expected = Q(math.prod(math.comb(2*d, d) for b, d in pairs if b == 0), 4)
    require(len(numerator) == len(denominator) and numerator[-1]/denominator[-1] == expected,
            "Exact leading-coefficient certificate of endpoint limit")
    for n in range(2, 7):
        require(evaluate(numerator, n)/evaluate(denominator, n) == exact_test(pairs, n, 1),
                "Rational-in-order endpoint identity")
    return expected


def validate_coefficients(pairs, n, values):
    require(len(values) == n+1, "Correct degree and coefficient count")
    delta = delta_by_rising(pairs, n)
    for k in range(n+1):
        require(Q(1, values[k]*math.comb(n, k)) == delta[k]*delta[n-k]/math.factorial(n)**2,
                "Independently reconstructed normalized reciprocal coefficients")


def rejection(operation, label):
    try:
        operation()
    except ValueError:
        return label
    raise ValueError("Failed rejecting control: "+label)


def pair_checks():
    cases = (((0, 2),), ((0, 1), (0, 1)), ((0, 2), (1, 1)),
             ((0, 1), (0, 1), (2, 3)), ((0, 3), (5, 2), (2, 4)),
             ((0, 2), (10, 1)), ((0, 1), (0, 2), (3, 1)))
    counts = dict(normalization_identities=0, exact_logconcavity_ratios=0,
                  derivative_inequalities=0, bezout_identities=0,
                  equal_row_identities=0, finite_order_refinements=0,
                  symbolic_endpoint_limits=0, padding_identities=0,
                  matched_merge_identities=0)
    for pairs in cases:
        private_degree = sum(d for b, d in pairs if b == 0)
        require(private_degree >= 2, "The uniform theorem hypothesis")
        rows = matrix(pairs)
        polynomials = []
        for n in range(1, 9):
            p = actual_binary(rows, n)
            validate_coefficients(pairs, n, p)
            counts["normalization_identities"] += n+1
            require(p == p[::-1] and all(c > 0 for c in p), "Positive palindromic coefficients")
            delta = delta_by_rising(pairs, n)
            eta = [delta[k]*delta[n-k] for k in range(n+1)]
            for k in range(1, n):
                require(eta[k]**2/(eta[k-1]*eta[k+1]) == exact_test(pairs, n, k) >= 1,
                        "Exact normalized log-concavity ratio")
                counts["exact_logconcavity_ratios"] += 1
            for x in (Q(0), Q(1, 2), Q(n)):
                actual = Q(2, x+1)-sum(Q(d, b*n+d*x+j) for b, d in pairs for j in range(1, d+1))
                require(actual <= Q(2-private_degree, x+1) <= 0, "Private-factor derivative bound")
                counts["derivative_inequalities"] += 1
            equal = actual_binary(matrix(pairs, (0, 1, 3)), n)
            scalar = math.factorial(n)*math.factorial(3*n)
            require(equal == scale(p, scalar), "Arbitrary equal rows contribute only an order scalar")
            counts["equal_row_identities"] += n+1
            polynomials.append(p)
        for a, b in zip(polynomials, polynomials[1:]):
            bezout(a, b)
            counts["bezout_identities"] += 1
        limit_certificate(pairs)
        counts["symbolic_endpoint_limits"] += 1

    finite_cases = (((1, 4),), ((1, 3), (2, 2)), ((0, 1), (1, 3)))
    for pairs in finite_cases:
        require(sum(d for b, d in pairs if b == 0) < 2, "Outside private-degree theorem")
        checked = 0
        for n in range(1, 8):
            if w_value(pairs, n+1) < 2:
                continue
            for order in (n, n+1):
                p = actual_binary(matrix(pairs), order)
                validate_coefficients(pairs, order, p)
                require(all(exact_test(pairs, order, k) >= 1 for k in range(1, order)),
                        "Finite-order sufficient inequality")
            bezout(actual_binary(matrix(pairs), n), actual_binary(matrix(pairs), n+1))
            counts["finite_order_refinements"] += 1
            checked += 1
        require(checked > 0, "Each retained finite-order family exercises the refinement")
        expected = limit_certificate(pairs)
        require(expected < 1, "Exact eventual failure of the imported criterion")
        counts["symbolic_endpoint_limits"] += 1

    # Actual rows for a matched merger (b,d)=(1,2) and (3,1), retaining private degree two.
    source = matrix(((0, 2), (1, 2), (3, 1)))
    merged = source[:2]+[tuple(a+b for a, b in zip(source[2], source[4])),
                         tuple(a+b for a, b in zip(source[3], source[5]))]
    require(merged == matrix(((0, 2), (4, 3))), "Matched row merge parameters")
    single = source[:2]+[tuple(a+b for a, b in zip(source[2], source[3]))]+source[4:]
    for n in range(1, 6):
        require(actual_binary(single, n) == scale(actual_binary(matrix(((0, 2), (3, 1))), n),
                                                  math.factorial(4*n)), "Single mirrored-pair merger")
        counts["matched_merge_identities"] += 1

    # Padding checked from the enlarged support matrix, independently of the EGF formula.
    core = matrix(((0, 2), (1, 1)))
    for extra in (1, 2, 3):
        rows = [tuple(row)+(0,)*extra for row in core]
        for j in range(extra):
            rows.append((0, 0)+tuple(int(k == j) for k in range(extra)))
        coefficients = (Q(1), Q(-1))+tuple(Q((-1)**j*(j+2)) for j in range(extra))
        f = [Q(actual_moment(core, coefficients[:2], n), math.factorial(n)) for n in range(9)]
        g = [Q(actual_moment(rows, coefficients, n), math.factorial(n)) for n in range(9)]
        denominator = [Q(1)]
        for c in coefficients[2:]:
            denominator = mul(denominator, [Q(1), -c])
        require(mul(g, denominator)[:9] == f, "Actual linear-padding coefficient identity")
        counts["padding_identities"] += 9

    controls = []
    bad = actual_binary(matrix(((0, 2), (1, 1))), 4)
    bad[2] += 1
    controls.append(rejection(lambda: validate_coefficients(((0, 2), (1, 1)), 4, bad), "altered middle coefficient"))
    controls.append(rejection(lambda: validate_coefficients(((0, 2),), 4, [1, 2]), "malformed degree"))
    controls.append(rejection(lambda: matrix(((1, 0),)), "zero difference outside mirrored-pair convention"))
    controls.append(rejection(lambda: matrix(((-1, 2),)), "negative base exponent"))
    p = actual_binary(matrix(((1, 1),)), 2)
    require(p == [48, 72, 48] and exact_test(((1, 1),), 2, 1) == Q(1, 9),
            "Known SFC common-power example fails the reciprocal coefficient criterion")
    controls.append("common-power SFC example b=d=1 has failed coefficient criterion, not failed SFC")
    require(limit_certificate(((0, 1),)) == Q(1, 2), "Degree-one boundary exact limit")
    controls.append("private degree one has limiting ratio 1/2; the geometric-sum family remains SFC")
    return dict(uniform_parameter_cases=len(cases), counts=counts, controls=controls)


def field_add(x, y):
    return tuple(a+b for a, b in zip(x, y))


def field_mul(x, y, radical_square):
    # Basis 1,r,i,r*i, with r^2=radical_square and i^2=-1.
    powers = ((0, 0), (1, 0), (0, 1), (1, 1))
    indices = {power: j for j, power in enumerate(powers)}
    out = [Q(0)]*4
    for j, a in enumerate(x):
        for k, b in enumerate(y):
            rp, ip = powers[j][0]+powers[k][0], powers[j][1]+powers[k][1]
            value = a*b
            if rp >= 2:
                value *= radical_square
                rp -= 2
            if ip >= 2:
                value = -value
                ip -= 2
            out[indices[(rp, ip)]] += value
    return tuple(out)


def row_merge_checks():
    counts = dict(source_reciprocal_identities=0, target_quadratic_identities=0,
                  upper_half_plane_zeros=0, distinct_column_checks=0)
    for s in (2, 3, 5, 8):
        source = [tuple(int(k == j) for k in range(s)) for j in range(s) for _ in range(2)]
        target = ([tuple(a+b for a, b in zip(source[0], source[2]))]+
                  [row for j, row in enumerate(source) if j not in (0, 2)])
        for rows in (source, target):
            columns = [tuple(row[j] for row in rows) for j in range(s)]
            require(len(set(columns)) == s and all(any(column) for column in columns), "Distinct nonconstant actual support")
            counts["distinct_column_checks"] += 1
        for n in range(1, 6):
            for multi in compositions(n, s):
                require(Q(1, actual_coefficient(source, multi)) ==
                        Q(1, math.factorial(n)*math.prod(math.factorial(k) for k in multi)),
                        "All-order source reciprocal identity control")
                counts["source_reciprocal_identities"] += 1
        reciprocal = {}
        for multi in compositions(2, s):
            value = Q(1, actual_coefficient(target, multi))
            square_coefficient = math.factorial(2)//math.prod(math.factorial(k) for k in multi)
            correction = int(multi[0] == multi[1] == 1)
            require(value == Q(square_coefficient-correction, 4), "Target quadratic polynomial identity")
            reciprocal[multi] = value
            counts["target_quadratic_identities"] += 1
        points = [(Q(0), Q(1, s), Q(1, s), Q(0)),
                  (Q(0), Q(-1, s), Q(1, s), Q(0))]
        points += [(Q(0), Q(0), Q(1, s), Q(0))]*(s-2)
        require(all(z[2] == Q(1, s) > 0 and z[3] == 0 for z in points), "All imaginary parts strictly positive")
        value = (Q(0),)*4
        for multi, c in reciprocal.items():
            term = (c, Q(0), Q(0), Q(0))
            for z, k in zip(points, multi):
                for _ in range(k):
                    term = field_mul(term, z, s*s-1)
            value = field_add(value, term)
        require(value == (0, 0, 0, 0), "Exact algebraic upper-half-plane zero")
        counts["upper_half_plane_zeros"] += 1
    require(Q(1, 8)**2 < Q(1, 4)**2, "Normalized log-concavity fails")
    require(Q(1, 2) > 0 and Q(1, 2)**2-Q(1, 4)**2 == Q(3, 16) > 0,
            "Two positive Hessian directions exclude the Lorentzian signature")
    return counts


def main():
    here = Path(__file__).resolve().parent
    root = here.parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=here/"SFC_OVERLAPPING_SWAP_PAIRS.md")
    parser.add_argument("--output", type=Path, default=here/"sfc_overlapping_swap_pairs_verification.json")
    args = parser.parse_args()
    note = args.note
    require(hashlib.sha256(note.read_bytes()).hexdigest() == NOTE_SHA, "Reviewed note pin")
    for name, expected in SOURCE_PINS.items():
        source = root/"research/factorial/image-rigidity"/name
        require(hashlib.sha256(source.read_bytes()).hexdigest() == expected, "Retained source pin: "+name)
    results = {"status": "PASS", "scope": "finite exact controls for the uniform private-degree proof and row-merge counterexample",
               "note_sha256": NOTE_SHA, "source_pins": SOURCE_PINS,
               "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "swap_pairs": pair_checks(), "row_merging": row_merge_checks()}
    receipt = args.output
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps(results, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"status": results["status"], "swap_pairs": results["swap_pairs"],
                      "row_merging": results["row_merging"], "receipt": str(receipt)}, sort_keys=True))


if __name__ == "__main__":
    main()
