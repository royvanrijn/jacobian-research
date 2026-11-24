#!/usr/bin/env python3
"""Exact finite controls for separating-row SFC amplification.

The universal theorem is proved in the pinned note. This script verifies
finite counting, factorial-gap and homogeneous-ideal certificates, including
negative controls. It does not infer an all-order theorem from finite tests.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
import importlib.util
from itertools import product
import json
from math import factorial, prod
from pathlib import Path

if not __debug__:
    raise SystemExit("Refusing optimized Python: dependency assertions must remain enabled.")

HERE = Path(__file__).resolve().parent
NOTE_NAME = "SFC_SEPARATING_ROW_AMPLIFICATION.md"
NOTE_PIN = "f10dc8c30e6a7f209bd2bca447ab3848e0c33c821958e49fb17e72945bbe2d7b"
DEPENDENCIES = {
    "SFC_BALANCED_PRIVATE_AMPLIFICATION.md": "f4ebaa9028cb04258cf65be178685650d53e2852644ed0678f1a655d134c065e",
    "verify_sfc_balanced_private_amplification.py": "eca0b5a325b31b347a058b1a6bbaadfcb5bc34e3b5a900d45ff42cb4da0ecca4",
    "sfc_balanced_private_amplification_source_ledger.json": "5a54143df4353a4b031a0c434d12262586e206e677d0923e173bdd6f399b123d",
}


def require(test, message):
    if not test:
        raise AssertionError(message)


@lru_cache(maxsize=None)
def compositions(m, s):
    if s == 1:
        return ((m,),)
    return tuple((j,) + k for j in range(m + 1)
                 for k in compositions(m - j, s - 1))


@lru_cache(maxsize=None)
def fac(n):
    return factorial(n)


def mult(k):
    return fac(sum(k)) // prod(fac(x) for x in k)


def qrecord(x):
    x = Q(x)
    return {"numerator": x.numerator, "denominator": x.denominator}


def profiles(matrix):
    s = len(matrix[0])
    require(all(len(row) == s for row in matrix), "Ragged exponent matrix")
    require(all(isinstance(x, int) and x >= 0 for row in matrix for x in row),
            "Nonintegral or negative exponent")
    return [tuple(sorted(row[j] for row in matrix if row[j])) for j in range(s)]


def balanced(matrix):
    p = profiles(matrix)
    require(all(x == p[0] for x in p), "Unbalanced pure profiles")
    return p[0]


def separating(code, layers):
    return len({tuple(word[l] for l in layers) for word in code}) == len(code)


def certificate(code, capacities, packing, threshold=True):
    s = len(code)
    require(s > 0 and len(capacities) > 0, "Empty code or layer list")
    require(all(len(word) == len(capacities) for word in code), "Ragged code")
    require(all(isinstance(d, int) and d >= 1 for d in capacities), "Invalid capacity")
    loads = [Q(0)] * len(capacities)
    W = Q(0)
    for layers, weight in packing:
        weight = Q(weight)
        require(weight >= 0 and len(set(layers)) == len(layers), "Invalid packing term")
        require(all(0 <= l < len(capacities) for l in layers), "Invalid layer index")
        if weight:
            require(separating(code, layers), "Positive-weight layer set does not separate")
        W += weight
        for l in layers:
            loads[l] += weight
    require(all(x <= d for x, d in zip(loads, capacities)), "Packing exceeds capacity")
    if threshold:
        require(W >= 6 * s, "Insufficient total amplification")
    return {"codewords": [list(x) for x in code], "capacities": capacities,
            "packing": [{"layers": list(B), "weight": qrecord(w)} for B, w in packing],
            "loads": [qrecord(x) for x in loads], "total_weight": qrecord(W)}


def code_matrix(code, capacities):
    return [[d * int(word[l] == symbol) for word in code]
            for l, d in enumerate(capacities)
            for symbol in sorted({word[l] for word in code})]


def forms(matrix, maximum):
    balanced(matrix)
    s = len(matrix[0])
    out = {}
    for m in range(1, maximum + 1):
        H = prod(fac(m * row[0]) for row in matrix)
        out[m] = {}
        for k in compositions(m, s):
            mixed = prod(fac(sum(a * b for a, b in zip(row, k))) for row in matrix)
            out[m][k] = Q(mult(k) * mixed, H)
    return out


def literal_expansion(matrix, maximum):
    """Multiply actual labelled monomials, independently of the composition formula."""
    s = len(matrix[0])
    features = [(tuple(int(j == h) for h in range(s)), tuple(row[j] for row in matrix))
                for j in range(s)]
    current = {((0,) * s, (0,) * len(matrix)): 1}
    expected = forms(matrix, maximum)
    count = 0
    for m in range(1, maximum + 1):
        nxt = {}
        for (k, alpha), value in current.items():
            for label, beta in features:
                key = (tuple(x + y for x, y in zip(k, label)),
                       tuple(x + y for x, y in zip(alpha, beta)))
                nxt[key] = nxt.get(key, 0) + value
        current = nxt
        H = prod(fac(m * row[0]) for row in matrix)
        actual = {}
        for (k, alpha), value in current.items():
            actual[k] = actual.get(k, Q(0)) + Q(value * prod(fac(x) for x in alpha), H)
        require(actual == expected[m], "Literal factorial expansion differs")
        count += len(actual)
    return count


def counting_controls():
    counts = {"literal_word_injections": 0, "composition_inequalities": 0,
              "stretched_multinomial_inequalities": 0, "active_subcodes": 0}
    for s in range(2, 7):
        L = (s - 1).bit_length()
        code = tuple(tuple((j >> l) & 1 for l in range(L)) for j in range(s))
        for m in range(1, 4):
            seen = set()
            for word in product(range(s), repeat=m):
                image = tuple(tuple(code[j][l] for j in word) for l in range(L))
                require(image not in seen, "Word map is not injective")
                seen.add(image)
                counts["literal_word_injections"] += 1
        for mask in range(1, 1 << s):
            subcode = tuple(code[j] for j in range(s) if mask >> j & 1)
            require(separating(subcode, tuple(range(L))), "Active boundary lost separation")
            counts["active_subcodes"] += 1
        for m in range(1, 7):
            for k in compositions(m, s):
                coarsened = [tuple(sum(k[j] for j in range(s) if code[j][l] == a)
                                   for a in (0, 1)) for l in range(L)]
                require(mult(k) <= prod(mult(q) for q in coarsened), "Counting bound failed")
                counts["composition_inequalities"] += 1
                for q in coarsened:
                    for D in (1, 2, 4):
                        require(mult(tuple(D * x for x in q)) >= mult(q) ** D,
                                "Stretched multinomial bound failed")
                        counts["stretched_multinomial_inequalities"] += 1
    return counts


def row_gap_controls():
    count = 0
    for row in product(range(4), repeat=3):
        if not any(row):
            Ds = (1,)
        else:
            Ds = range(1, min(x for x in row if x) + 1)
        for D in Ds:
            for m in range(1, 6):
                for k in compositions(m, 3):
                    K = sum(v for a, v in zip(row, k) if a)
                    A = sum(a * v for a, v in zip(row, k))
                    lhs = prod(fac(m * a) ** v for a, v in zip(row, k)) * fac(D * K) ** m
                    rhs = fac(A) ** m * fac(D * m) ** K
                    require(lhs >= rhs, "Row factorial-gap inequality failed")
                    count += 1
    return count


def intrinsic_controls():
    count = 0
    records = []
    matrices = [
        ("triangle", [[1, 0, 1], [1, 1, 0], [0, 1, 1]], [0, 1]),
        ("unequal_positive_entries", [[1, 0, 2], [2, 1, 0], [0, 2, 1]], [0, 1]),
        ("binary_square", [[1, 1, 0, 0], [0, 0, 1, 1],
                            [1, 0, 1, 0], [0, 1, 0, 1]], [0, 2]),
    ]
    for name, unscaled, selected in matrices:
        for D in (1, 2, 3):
            matrix = [[D * x for x in row] for row in unscaled]
            profile = balanced(matrix)
            s = len(matrix[0])
            ell = len(selected)
            code = tuple(tuple(int(matrix[h][j] > 0) for h in selected) for j in range(s))
            require(separating(code, tuple(range(ell))), "Selected incidence rows do not separate")
            for m in range(1, 7):
                H = prod(fac(m * a) for a in profile)
                for k in compositions(m, s):
                    K = [sum(k[j] for j in range(s) if row[j]) for row in matrix]
                    mixed = prod(fac(sum(a * v for a, v in zip(row, k))) for row in matrix)
                    binary = Q(prod(fac(D * x) for x in K), fac(D * m) ** len(profile))
                    ratio = Q(mixed, H)
                    require(ratio <= binary, "Conditional Jensen comparison failed")
                    for h in selected:
                        require(binary <= Q(fac(D * K[h]) * fac(D * (m - K[h])), fac(D * m)),
                                "Bounded-load concentration failed")
                    require(ratio ** ell * mult(k) ** D <= 1, "Intrinsic suppression failed")
                    count += 1
            records.append({"name": name, "matrix": matrix, "separating_rows": selected})
    return {"cases": count, "inputs": records}


def merger_controls():
    code = [(0, 0), (0, 1), (1, 0), (1, 1)]
    groups = [list(range(4)), list(range(4))]
    # Merge equal final symbols one pair at a time, crossing the two layers.
    steps = [(0, 0, 1), (1, 0, 2), (0, 2, 3), (1, 1, 3)]
    outputs = []
    for stage in range(len(steps) + 1):
        current = [tuple(groups[l][j] for l in range(2)) for j in range(4)]
        cert = certificate(current, [24, 24], [((0, 1), Q(24))])
        matrix = code_matrix(current, [24, 24])
        balanced(matrix)
        outputs.append({"stage": stage, "certificate": cert, "matrix": matrix})
        if stage < len(steps):
            l, a, b = steps[stage]
            old, new = groups[l][b], groups[l][a]
            groups[l] = [new if x == old else x for x in groups[l]]
    require(all((groups[l][j] == groups[l][h]) == (code[j][l] == code[h][l])
                for l in range(2) for j in range(4) for h in range(4)), "Wrong final merger")
    return outputs


def asymptotic_controls():
    """Exact finite-product certificates for the leading ratio, without numerical limits."""
    records = []
    for a, b in [((2 * d, d), (d, 2 * d)) for d in (1, 2, 6)] + [
            ((2, 0, 1), (0, 1, 2)), ((1, 1, 0), (0, 1, 1)),
            ((2, 1, 0), (1, 2, 0))]:
        require(sorted(a) == sorted(b), "Asymptotic control is not balanced")
        numerator_factors, denominator_factors = [], []
        constant = 1
        beta = sum(y for x, y in zip(a, b) if x == 0)
        expected = Q(1)
        for x, y in zip(a, b):
            if x == 0:
                constant *= fac(y)
                expected *= fac(y)
            else:
                delta = y - x
                expected *= Q(x) ** delta
                if delta >= 0:
                    numerator_factors.extend((x, t) for t in range(1, delta + 1))
                else:
                    denominator_factors.extend((x, -t) for t in range(-delta))
        degree_difference = len(numerator_factors) - len(denominator_factors)
        leading = Q(constant * prod(x for x, _ in numerator_factors),
                    prod(x for x, _ in denominator_factors))
        require(degree_difference == -beta and leading == expected,
                "Incorrect leading power or coefficient")
        for m in range(1, 8):
            rational_product = Q(constant * prod(x * m + t for x, t in numerator_factors),
                                 prod(x * m + t for x, t in denominator_factors))
            actual = Q(prod(fac((m - 1) * x + y) for x, y in zip(a, b)),
                       prod(fac(m * x) for x in a))
            require(actual == rational_product, "Factorial ratio product mismatch")
            if len(a) == 2:
                require(beta == 0 and actual > leading, "Identical-pattern lower bound failed")
        records.append({"a": a, "b": b, "beta": beta, "leading_constant": qrecord(leading),
                        "numerator_factors_slope_intercept": numerator_factors,
                        "denominator_factors_slope_intercept": denominator_factors,
                        "constant_factor": constant, "degree_difference": degree_difference})
    return records


def negative_controls(helper):
    rejected = []
    tests = [
        ("nonseparating_code", lambda: certificate([(0,), (0,)], [12], [((0,), Q(12))])),
        ("over_capacity", lambda: certificate([(0,), (1,)], [11], [((0,), Q(12))])),
        ("below_threshold", lambda: certificate([(0,), (1,)], [11], [((0,), Q(11))])),
        ("unbalanced_endpoint_profile", lambda: balanced([[2, 0], [0, 1], [0, 1]])),
    ]
    for name, call in tests:
        try:
            call()
        except AssertionError:
            rejected.append(name)
        else:
            raise AssertionError("Invalid input accepted: " + name)
    # A common marker gives no C_k^(-D) suppression, even with distinct balanced base terms.
    matrix = [[1, 0], [1, 1], [0, 1], [12, 12]]
    mixed = forms(matrix, 2)[2][(1, 1)]
    require(mixed == 1 and mixed > Q(1, 2 ** 11), "Collapsed-code negative control changed")
    require(not separating([(1, 1), (1, 1)], (0, 1)), "Identical incidence patterns accepted")
    # Repeated exponent columns produce a genuinely deficient coefficient system.
    collision = forms([[1, 1]], 3)
    _, columns, rows, _ = helper.macaulay_rows(collision, 2, 2, 104729)
    selected, _ = helper.select_independent_rows(rows, 104729)
    require(len(selected) < len(columns), "Rank checker accepted a deficient system")
    return {"rejected_inputs": rejected, "collapsed_code_mixed_coefficient": qrecord(mixed),
            "deficient_rank": len(selected), "deficient_column_count": len(columns),
            "scope": "These controls reject invalid certificates; they are not SFC counterexamples."}


def rank_certificate(helper, name, matrix, n, prime=104729):
    s = len(matrix[0])
    actual = forms(matrix, n + s - 1)
    L, columns, rows, labels = helper.macaulay_rows(actual, s, n, prime)
    selected, pivots = helper.select_independent_rows(rows, prime)
    require(len(selected) == len(columns), "Macaulay rank failed: " + name)
    det = helper.determinant_mod([rows[i] for i in selected], prime)
    require(det != 0, "Selected full minor vanishes")
    return {"name": name, "exponent_matrix": matrix, "starting_order": n,
            "prime": prime, "homogeneous_degree": L, "row_count": len(rows),
            "column_count": len(columns), "selected_rows": selected,
            "selected_row_labels": [labels[i] for i in selected],
            "column_monomials": [list(k) for k in columns], "pivot_columns": pivots,
            "full_minor_determinant_mod_prime": det,
            "matrix_sha256": sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest(),
            "normalized_forms": {str(m): [{"powers": list(k), "coefficient": qrecord(v)}
                                           for k, v in actual[m].items()]
                                 for m in range(n, n + s)},
            "implication": "All monomials of the recorded degree lie in the rational moment ideal; the only common zero is the origin."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=HERE / NOTE_NAME)
    parser.add_argument("--output", type=Path, default=HERE / "sfc_separating_row_amplification_verification.json")
    args = parser.parse_args()
    require(sha256(args.note.read_bytes()).hexdigest() == NOTE_PIN, "Reviewed note pin mismatch")
    for name, pin in DEPENDENCIES.items():
        require(sha256((HERE / name).read_bytes()).hexdigest() == pin, "Dependency pin mismatch: " + name)
    path = HERE / "verify_sfc_balanced_private_amplification.py"
    spec = importlib.util.spec_from_file_location("balanced_rank_helper", path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)

    binary = [(0, 0), (0, 1), (1, 0), (1, 1)]
    parity = [(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)]
    certificate_records = [certificate(binary, [24, 24], [((0, 1), Q(24))]),
                           certificate(parity, [16, 16, 16],
                                       [((0, 1), Q(8)), ((0, 2), Q(8)), ((1, 2), Q(8))])]
    binary_matrix = code_matrix(binary, [24, 24])
    parity_matrix = code_matrix(parity, [16, 16, 16])
    triangle = [[36, 0, 36], [36, 36, 0], [0, 36, 36]]
    unequal = [[36, 0, 72], [72, 36, 0], [0, 72, 36]]
    all_matrices = [binary_matrix, parity_matrix, triangle, unequal]
    literal = sum(literal_expansion(A, 4) for A in all_matrices)
    rank_cases = [("binary_shared", binary_matrix, 1), ("parity_fractional", parity_matrix, 1),
                  ("intrinsic_triangle", triangle, 2), ("intrinsic_unequal_entries", unequal, 3)]
    rank_records = [rank_certificate(helper, *case) for case in rank_cases]
    receipt = {
        "status": "PASS", "schema_version": 1,
        "repository_base": "19eb9912fb60142f8b620dac24896931d4ff3e46",
        "note_sha256": NOTE_PIN, "dependency_sha256": DEPENDENCIES,
        "scope": "Finite exact controls and explicitly listed all-coefficient Macaulay certificates; universal all-order results require the pinned proof and its inspected Turan input.",
        "packing_certificates": certificate_records,
        "counting_controls": counting_controls(), "row_gap_controls": row_gap_controls(),
        "intrinsic_controls": intrinsic_controls(), "merger_stages": merger_controls(),
        "asymptotic_barrier_certificates": asymptotic_controls(),
        "literal_phase_coefficients": literal, "rank_certificates": rank_records,
        "negative_controls": negative_controls(helper),
        "full_SFC_solved": False, "arbitrary_row_merging_proved": False,
        "private_decoration_removal_proved": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "row_gap_controls": receipt["row_gap_controls"],
                      "literal_phase_coefficients": literal,
                      "macaulay_certificates": len(rank_records),
                      "largest_minor": max(r["column_count"] for r in rank_records),
                      "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
