#!/usr/bin/env python3
"""Exact finite controls for WEIGHTED_RADIAL_LINEAR_SFC_TRANSFER.md.

Standard library only. Formal inversion uses coefficientwise composition,
independently of the Lagrange sum and actual sparse moment expansions. The
written proof, not finite extrapolation, supplies unrestricted quantifiers.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys


DIRECTORY = Path("research/factorial/image-rigidity")
INSPECTED_BASE = "19eb9912fb60142f8b620dac24896931d4ff3e46"
SOURCE_PINS = {
    "WEIGHTED_RADIAL_LINEAR_SFC_TRANSFER.md": "f9ee0dd06d5e6bc85a5b8d35aacb6d55f5c8ec18e483a59c00d3275177495645",
    "POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md": "40278040a10d4813c473e8256b868e4fb6d42f20912972b801136398436609d6",
    "RIGIDITY_CRITICAL_VALUES.tex": "8a2b6faa7996398d8c0474b34550a56cab6a771436a49fb2d3dff7a393e9c4dd",
    "sources/rigidity.typ": "8a6e152c797c227529e88572db140fcfb5b6c3ad61650c2a94fc6304cb334af6",
    "UNRESTRICTED_SFC_MECHANISMS.tex": "f8b3dfc5947601fb002aee0ba942b0fc20b7ef46ca8e0d4f69b5a3e5e1361785",
}
PARAMETER_CASES = (
    ((1,), (2,)), ((2,), (-3,)), ((4,), (1,)),
    ((1, 2), (5, -3)), ((2, 1), (2, -3)),
    ((2, 2), (1, -1)), ((2, 3), (2, 2)),
    ((1, 3), (0, -2)), ((3, 2), (2, 0)),
    ((1, 2, 3), (1, -2, 3)), ((2, 1, 2), (1, 1, -2)),
    ((2, 3, 1), (0, 2, -2)), ((2, 3), (0, 0)),
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(repo):
    result = {}
    for name, expected in SOURCE_PINS.items():
        path = repo / DIRECTORY / name
        require(path.is_file(), "Missing proof/source: " + name)
        actual = file_hash(path)
        require(actual == expected, "Changed proof/source: " + name)
        result[name] = actual
    return result


def trim(poly):
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def coefficient(poly, degree):
    return poly[degree] if degree < len(poly) else 0


def multiply(a, b, degree):
    result = [0] * min(len(a) + len(b) - 1, degree + 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b[:max(0, degree + 1 - i)]):
                if bj:
                    result[i + j] += ai * bj
    return trim(result)


def compose(a, b, degree):
    result = [0]
    for term in reversed(a):
        result = multiply(result, b, degree)
        result[0] += term
    return trim(result)


def inverse_by_composition(a, degree):
    require(a[:2] == [0, 1], "Normalized polynomial required")
    b = [0, 1] + [0] * (degree - 1)
    for k in range(2, degree + 1):
        b[k] = -coefficient(compose(a, b[:k], k), k)
    require(compose(a, b, degree) == [0, 1], "a(b)=z")
    require(compose(b, a, degree) == [0, 1], "b(a)=z")
    return b


def validate(q, shapes, coefficients):
    require(type(q) is int and q >= 1, "q is a positive integer")
    require(len(shapes) == len(coefficients) >= 1, "Nonempty equal dimensions")
    require(all(type(p) is int and p >= 1 for p in shapes),
            "Every Gamma shape is a positive integer")


def inverse_polynomial(q, shapes, coefficients):
    validate(q, shapes, coefficients)
    result = [1]
    degree = q * sum(shapes)
    for p, c in zip(shapes, coefficients):
        factor = [1] + [0] * (q - 1) + [-c]
        for _ in range(p):
            result = multiply(result, factor, degree)
    return [0] + result


def rising(p, k):
    result = 1
    for offset in range(k):
        result *= p + offset
    return result


def sparse_multiply(a, b):
    result = {}
    for alpha, ca in a.items():
        for beta, cb in b.items():
            gamma = tuple(x + y for x, y in zip(alpha, beta))
            result[gamma] = result.get(gamma, 0) + ca * cb
    return {exponent: c for exponent, c in result.items() if c}


def weighted_phase(q, shapes, coefficients):
    validate(q, shapes, coefficients)
    result = {}
    for i, c in enumerate(coefficients):
        if c:
            alpha = [q * p for p in shapes]
            alpha[i] += 1
            result[tuple(alpha)] = c
    return result


def actual_moments(polynomial, shapes, maximum):
    power = {(0,) * len(shapes): 1}
    result = []
    for m in range(maximum + 1):
        result.append(sum(c * math.prod(rising(p, e) for p, e in zip(shapes, alpha))
                          for alpha, c in power.items()))
        if m != maximum:
            power = sparse_multiply(power, polynomial)
    return result


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for k in range(total + 1):
            for rest in compositions(total - k, length - 1):
                yield (k,) + rest


def lagrange_sum(q, shapes, coefficients, m):
    total = Fraction(0)
    for multiindex in compositions(m, len(shapes)):
        term = Fraction(1)
        for p, c, k in zip(shapes, coefficients, multiindex):
            term *= Fraction(rising(p * (q * m + 1), k) * c ** k,
                             math.factorial(k))
        total += term
    return total / (q * m + 1)


def derivative_factorization(q, shapes, coefficients):
    a = inverse_polynomial(q, shapes, coefficients)
    actual = trim([k * value for k, value in enumerate(a) if k])
    groups = {}
    for p, c in zip(shapes, coefficients):
        if c:
            groups[c] = groups.get(c, 0) + p
    groups = sorted(groups.items())
    h, upper = len(groups), len(a) - 2
    common = [1]
    for c, p in groups:
        for _ in range(p - 1):
            common = multiply(common, [1] + [0] * (q - 1) + [-c], upper)
    residual = [1]
    for c, _ in groups:
        residual = multiply(residual, [1] + [0] * (q - 1) + [-c], q * h)
    residual += [0] * (q * h + 1 - len(residual))
    for index, (c, p) in enumerate(groups):
        term = [1]
        for j, (d, _) in enumerate(groups):
            if j != index:
                term = multiply(term, [1] + [0] * (q - 1) + [-d], q * h)
        for k, value in enumerate(term):
            residual[k + q] -= q * p * c * value
    trim(residual)
    require(multiply(common, residual, max(0, upper)) == actual,
            "Complete derivative factorization, including repeated factors")
    require(len(residual) - 1 == q * h, "Residual degree equals qh")
    leading = (-1) ** h * math.prod(c for c, _ in groups) * (1 + q * sum(p for _, p in groups))
    require(residual[-1] == leading, "Nonzero highest derivative-residual coefficient")
    return {"groups": groups, "degree_a": len(a) - 1,
            "degree_common_zero_fiber_factor": len(common) - 1,
            "degree_residual": len(residual) - 1,
            "common_factor": common, "residual": residual}


def family_checks():
    rows = []
    identities = lagrange = structural = inactive = windows = marked = 0
    maximum = 4
    for q in (1, 2, 3, 4):
        for shapes, coefficients in PARAMETER_CASES:
            a = inverse_polynomial(q, shapes, coefficients)
            b = inverse_by_composition(a, q * maximum + 1)
            f = weighted_phase(q, shapes, coefficients)
            moments = actual_moments(f, shapes, maximum)
            for m, moment in enumerate(moments):
                factor = math.factorial(m) * (q * m + 1) * math.prod(
                    rising(p, q * p * m) for p in shapes)
                require(moment == factor * b[q * m + 1], "Gamma moment / inverse identity")
                identities += 1
                require(b[q * m + 1] == lagrange_sum(q, shapes, coefficients, m),
                        "Independent Lagrange sum / composition inverse identity")
                lagrange += 1
            for ell, value in enumerate(b):
                if ell % q != 1 % q:
                    require(value == 0, "Structural inverse zero class")
                    structural += 1
            require(moments[0] == b[1] == 1, "m=0 identity")
            active = [(p, c) for p, c in zip(shapes, coefficients) if c]
            h = len(set(c for c in coefficients if c))
            if active:
                for start in range(1, maximum - h + 2):
                    require(any(moments[start:start + h]), "Finite distinct-value window")
                    windows += 1
                active_shapes, active_coefficients = map(tuple, zip(*active))
                reduced = actual_moments(weighted_phase(q, active_shapes, active_coefficients),
                                         active_shapes, maximum)
                for m in range(maximum + 1):
                    missing = math.prod(rising(p, q * p * m)
                                        for p, c in zip(shapes, coefficients) if not c)
                    require(moments[m] == missing * reduced[m], "Inactive shape factors")
                    inactive += 1
            else:
                require(moments[1:] == [0] * maximum, "Zero phase boundary")
            power = {(0,) * len(shapes): 1}
            for m in range(maximum + 1):
                ordinary_marked = sum(c * math.prod(math.factorial(e + p - 1)
                                                    for e, p in zip(alpha, shapes))
                                      for alpha, c in power.items())
                require(ordinary_marked == math.prod(math.factorial(p - 1) for p in shapes) * moments[m],
                        "Exact fixed-marker identity")
                marked += 1
                if m < maximum:
                    power = sparse_multiply(power, f)
            rows.append({"q": q, "shapes": shapes, "coefficients": coefficients,
                         "h": h, "weighted_moments_through_four": moments,
                         "inverse_coefficients": b,
                         "derivative": derivative_factorization(q, shapes, coefficients)})
    return {"parameter_cases": len(rows), "moment_identities": identities,
            "independent_lagrange_identities": lagrange, "marked_identities": marked,
            "structural_zero_coefficients": structural,
            "inactive_variable_identities": inactive, "finite_window_controls": windows,
            "cases": rows}


def block_sum_phase(q, shapes, coefficients):
    validate(q, shapes, coefficients)
    dimension = sum(shapes)
    zero = (0,) * dimension
    common, linear = {zero: 1}, {}
    offset = 0
    for p, c in zip(shapes, coefficients):
        block = {}
        for j in range(offset, offset + p):
            alpha = list(zero)
            alpha[j] = 1
            block[tuple(alpha)] = 1
            if c:
                linear[tuple(alpha)] = c
        for _ in range(q * p):
            common = sparse_multiply(common, block)
        offset += p
    return sparse_multiply(common, linear)


def block_checks():
    inputs = ((1, (1, 2), (5, -3)), (2, (1, 2), (5, -3)),
              (1, (2, 2), (1, -1)), (1, (2, 2), (3, 3)),
              (1, (1, 2, 1), (1, -2, 0)), (2, (3,), (-2,)),
              (1, (2, 1), (0, -1)), (1, (1, 2), (0, 0)))
    rows, identities = [], 0
    for q, shapes, coefficients in inputs:
        g = block_sum_phase(q, shapes, coefficients)
        ordinary = actual_moments(g, (1,) * sum(shapes), 3)
        weighted = actual_moments(weighted_phase(q, shapes, coefficients), shapes, 3)
        require(ordinary == weighted, "Ordinary block expansion / weighted identity")
        identities += len(ordinary)
        support_formula = 0
        for i, c in enumerate(coefficients):
            if c:
                support_formula += math.comb((q + 1) * shapes[i], shapes[i] - 1) * math.prod(
                    math.comb((q + 1) * p - 1, p - 1) for j, p in enumerate(shapes) if j != i)
        active = sum(c != 0 for c in coefficients)
        h = len(set(c for c in coefficients if c))
        require(len(g) == support_formula >= active >= h, "Actual support count and window comparison")
        rows.append({"q": q, "shapes": shapes, "coefficients": coefficients,
                     "ordinary_variables": sum(shapes), "actual_support": len(g),
                     "support_formula": support_formula, "h": h,
                     "moments_through_three": ordinary})
    return {"cases": rows, "parameter_cases": len(rows), "moment_identities": identities}


def window_checks():
    rows = []
    for q, h, n in ((1, 1, 1), (1, 3, 4), (2, 1, 1), (2, 3, 2), (3, 2, 5), (5, 4, 7)):
        first, last = q * (n - 1) + 2, q * (n + h)
        selected = {q * m + 1 for m in range(n, n + h)}
        covered = {ell for ell in range(first, last + 1)
                   if ell % q != 1 % q or ell in selected}
        require(covered == set(range(first, last + 1)), "Complete inverse zero interval")
        require(first >= 2 and len(covered) == q * (h + 1) - 1 >= q * h,
                "Enough inverse zeros for the nonzero-critical-value bound")
        rows.append({"q": q, "h": h, "start": n, "first": first, "last": last,
                     "forced": len(covered), "required_upper_bound": q * h})
    return rows


def rejecting_checks():
    negatives = []
    for name, q, shapes, coefficients in (
        ("q_zero", 0, (1,), (1,)), ("q_noninteger", Fraction(3, 2), (1,), (1,)),
        ("q_boolean", True, (1,), (1,)), ("shape_zero", 1, (0,), (1,)),
        ("shape_boolean", 1, (True,), (1,)), ("shape_noninteger", 1, (Fraction(3, 2),), (1,)),
        ("empty_shapes", 1, (), ()), ("dimension_mismatch", 1, (1, 2), (1,))):
        try:
            inverse_polynomial(q, shapes, coefficients)
        except ValueError:
            negatives.append(name)
        else:
            raise ValueError("Malformed input accepted: " + name)
    q, shapes, coefficients = 2, (1, 2), (5, -3)
    f = weighted_phase(q, shapes, coefficients)
    weighted = actual_moments(f, shapes, 2)
    ordinary = actual_moments(f, (1, 1), 2)
    require(ordinary[1] == 0 and weighted[1] == -720,
            "Wrong-functional scope countercontrol is exactly 0 versus -720")
    negatives.append("ordinary_zero_does_not_imply_weighted_zero")
    b = inverse_by_composition(inverse_polynomial(q, shapes, coefficients), 5)
    factor = (q + 1) * math.prod(rising(p, q * p) for p in shapes)
    require(weighted[1] != factor * (b[q + 1] + 1), "Changed inverse coefficient rejected")
    negatives.append("changed_inverse_coefficient")
    wrong = 2 * math.prod(rising(p, q * p) for p in shapes) * b[q + 1]
    require(weighted[1] != wrong, "m+1 cannot replace qm+1")
    negatives.append("wrong_normalizer")
    bad_exponents = {(2, 6): 5, (1, 7): -3}
    bad_moment = actual_moments(bad_exponents, shapes, 1)[1]
    require(bad_moment != weighted[1], "Exponent/shape matching cannot be omitted")
    negatives.append("mismatched_exponent_shape")
    return {"passed": negatives, "count": len(negatives),
            "scope_control": {"q": q, "shapes": shapes, "coefficients": coefficients,
                              "ordinary_first_moment": ordinary[1],
                              "weighted_first_moment": weighted[1]}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sources = check_sources(args.repo_root.resolve())
    family, blocks = family_checks(), block_checks()
    windows, negatives = window_checks(), rejecting_checks()
    result = {
        "status": "PASS", "inspected_base": INSPECTED_BASE,
        "assurance": "Exact finite controls support the written all-order proof; no formal verification or external review is claimed.",
        "scope": "Integer-Gamma matched exponent family, and its ordinary disjoint-block-sum lift. Excludes ordinary unmarked unequal-exponent transfer and full SFC.",
        "source_sha256": sources, "checker_sha256": file_hash(Path(__file__).resolve()),
        "weighted_family": family, "ordinary_blocks": blocks,
        "window_controls": windows, "negative_controls": negatives,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "weighted_cases": family["parameter_cases"],
                      "weighted_moment_identities": family["moment_identities"],
                      "ordinary_block_cases": blocks["parameter_cases"],
                      "ordinary_block_identities": blocks["moment_identities"],
                      "negative_controls": negatives["count"], "output": str(args.output)}))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as error:
        print("FAIL: " + str(error), file=sys.stderr)
        raise SystemExit(1)
