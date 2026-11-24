#!/usr/bin/env python3
"""Exact finite certificates for balanced private SFC amplification.

Standard library only. The all-parameter proof uses the pinned note and
the inspected Turan theorem; finite Macaulay certificates are additional
evidence for their precisely listed supports and starting orders.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, factorial, isqrt, prod
from pathlib import Path


if not __debug__:
    raise SystemExit("Refusing optimized Python: exact assertions must stay enabled.")


DIRECTORY = Path("research/factorial/image-rigidity")
NOTE_NAME = "SFC_BALANCED_PRIVATE_AMPLIFICATION.md"
NOTE_PIN = "f4ebaa9028cb04258cf65be178685650d53e2852644ed0678f1a655d134c065e"
LEDGER_NAME = "sfc_balanced_private_amplification_source_ledger.json"
LEDGER_PIN = "5a54143df4353a4b031a0c434d12262586e206e677d0923e173bdd6f399b123d"
PDF_PIN = "ef098f0a0e5e27fd16721e44fb2197b163357de45b16111e7b82d179f36e7855"


def frac(x):
    x = F(x)
    return {"numerator": x.numerator, "denominator": x.denominator}


@lru_cache(maxsize=None)
def compositions(total, size):
    if size == 1:
        return ((total,),)
    return tuple((a,)+tail for a in range(total+1)
                 for tail in compositions(total-a, size-1))


def multinomial(k):
    return factorial(sum(k))//prod(factorial(j) for j in k)


def profiles(matrix):
    size = len(matrix[0])
    assert all(len(row) == size for row in matrix)
    assert all(isinstance(a, int) and a >= 0 for row in matrix for a in row)
    return [tuple(sorted(row[j] for row in matrix if row[j] > 0))
            for j in range(size)]


def require_balanced(matrix):
    p = profiles(matrix)
    assert all(x == p[0] for x in p), "Endpoint exponent profiles are not balanced"
    return p[0]


def balance(matrix):
    """Literal private copies of every other column's positive entries."""
    size = len(matrix[0])
    out = [list(row) for row in matrix]
    additions = []
    for j in range(size):
        for h in range(size):
            if h == j:
                continue
            for i, row in enumerate(matrix):
                exponent = row[h]
                if exponent:
                    added = [0]*size
                    added[j] = exponent
                    out.append(added)
                    additions.append({"target_column": j, "source_column": h,
                                      "source_row": i, "exponent": exponent,
                                      "new_row": len(out)-1})
    expected = tuple(sorted(a for row in matrix for a in row if a > 0))
    assert require_balanced(out) == expected
    return out, additions


def monomial_bases(matrix):
    return [{tuple(row[j] for row in matrix): F(1)}
            for j in range(len(matrix[0]))]


def poly_mul(a, b):
    out = {}
    for alpha, x in a.items():
        for beta, y in b.items():
            gamma = tuple(i+j for i, j in zip(alpha, beta))
            out[gamma] = out.get(gamma, F(0))+x*y
    return {k: v for k, v in out.items() if v}


def poly_pow(a, n):
    variables = len(next(iter(a)))
    out = {(0,)*variables: F(1)}
    for _ in range(n):
        out = poly_mul(out, a)
    return out


def moment(poly):
    return sum((c*prod(factorial(j) for j in alpha)
                for alpha, c in poly.items()), F(0))


def literal_phase_coefficients(bases, D, maximum):
    """Expand the actual augmented phase, retaining c-monomial labels."""
    size = len(bases)
    variables = len(next(iter(bases[0])))
    phase = {}
    for j, base in enumerate(bases):
        label = tuple(int(i == j) for i in range(size))
        marks = tuple(D*int(i == j) for i in range(size))
        for alpha, coeff in base.items():
            phase[(label, alpha+marks)] = coeff
    current = {((0,)*size, (0,)*(variables+size)): F(1)}
    outputs = [{(0,)*size: F(1)}]
    for m in range(1, maximum+1):
        nxt = {}
        for (k, alpha), coeff in current.items():
            for (l, beta), other in phase.items():
                key = (tuple(i+j for i, j in zip(k, l)),
                       tuple(i+j for i, j in zip(alpha, beta)))
                nxt[key] = nxt.get(key, F(0))+coeff*other
        current = {k: v for k, v in nxt.items() if v}
        out = {}
        for (k, alpha), coeff in current.items():
            out[k] = out.get(k, F(0))+coeff*prod(factorial(j) for j in alpha)
        outputs.append(out)
    return outputs


def normalized_forms(bases, D, maximum, check_literal=True):
    size = len(bases)
    variables = len(next(iter(bases[0])))
    powers = [[poly_pow(base, m) for m in range(maximum+1)] for base in bases]
    endpoints = [[moment(powers[j][m]) for j in range(size)]
                 for m in range(maximum+1)]
    assert all(all(x == row[0] for x in row) and row[0] > 0 for row in endpoints)
    literal = literal_phase_coefficients(bases, D, maximum) if check_literal else None
    forms = {0: {(0,)*size: F(1)}}
    controls = []
    for m in range(1, maximum+1):
        H = endpoints[m][0]
        forms[m] = {}
        for k in compositions(m, size):
            mixed = {(0,)*variables: F(1)}
            for j in range(size):
                mixed = poly_mul(mixed, powers[j][k[j]])
            shared = moment(mixed)
            assert 0 < shared <= H
            C = multinomial(k)
            marker = F(prod(factorial(D*j) for j in k), factorial(D*m))
            assert marker <= F(1, C**D)
            coefficient = C*shared/H*marker
            if literal is not None:
                actual = literal[m].get(k, F(0))/(H*factorial(D*m))
                assert actual == coefficient
            assert 0 < coefficient <= F(1, C**(D-1))
            if max(k) == m:
                assert coefficient == 1
            else:
                assert C >= m
            forms[m][k] = coefficient
        if m == 1:
            assert all(v == 1 for v in forms[m].values())
        controls.append({"m": m, "common_endpoint_moment": frac(H),
                         "coefficient_count": len(forms[m]),
                         "mixed_coefficient_sum": frac(sum(
                             v for k, v in forms[m].items() if max(k) != m))})
    if literal is not None:
        assert literal[0] == {(0,)*size: F(1)}
    return forms, controls


def form_receipt(poly):
    return [{"exponents": list(k), **frac(v)} for k, v in sorted(poly.items())]


def rational_constant_checks():
    rows = []
    for s in range(2, 81):
        R = F(s*(s+1), 2)*F(3*(s+1), 8*s)**s
        Rn = F((s+1)*(s+2), 2)*F(3*(s+2), 8*(s+1))**(s+1)
        ratio = F(3*(s+2)**2, 8*s*(s+1))*F(s*(s+2), (s+1)**2)**s
        assert Rn/R == ratio < 1
        assert 8*s*(s+1)-3*(s+2)**2 == (s-2)*(5*s+6) >= 0
        Q = F(s*(s+1), 2)*F(3, 8)**s
        Qn = F((s+1)*(s+2), 2)*F(3, 8)**(s+1)
        assert Qn/Q == F(3*(s+2), 8*s) <= F(3, 4)
        assert R <= F(243, 256) < 1
        assert Q <= F(27, 64) < 1
        for n in [2, 3, 7, 31]:
            D = 6*s
            error = F(comb(n+s-1, s-1), n**(D-1))
            # The rational quantity replaces e by the upper bound 3.
            lower = 2*F(s, 24*(n+s-1))**s
            assert error/lower <= R
        if s <= 8:
            rows.append({"s": s, "D": 6*s, "R": frac(R), "Q": frac(Q),
                         "R_successive_ratio": frac(ratio)})
    assert F(3)*F(9, 16)**2 == F(243, 256)
    assert F(3)*F(3, 8)**2 == F(27, 64)
    return {"rational_replay_s": [2, 80], "examples": rows,
            "no_floating_point_evaluation": True,
            "universal_inequality": "(s-2)(5s+6)>=0; s(s+2)<(s+1)^2 for s>=2",
            "universal_remaining_margin": frac(F(13, 256))}


def factorial_inequalities():
    jensen_count = stretch_count = 0
    for matrix in [
            [[2, 0, 1], [1, 2, 0], [0, 1, 2]],
            [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]],
            [[2, 2], [1, 1]],
    ]:
        size = len(matrix[0])
        require_balanced(matrix)
        for m in range(1, 6):
            for k in compositions(m, size):
                for row in matrix:
                    average = sum(a*b for a, b in zip(row, k))
                    lhs = factorial(average)**m
                    rhs = prod(factorial(m*row[j])**k[j] for j in range(size))
                    assert lhs <= rhs
                    values = [m*row[j] for j in range(size) for _ in range(k[j])]
                    old = prod(factorial(v) for v in values)
                    while max(values)-min(values) >= 2:
                        lo, hi = values.index(min(values)), values.index(max(values))
                        values[lo] += 1
                        values[hi] -= 1
                        new = prod(factorial(v) for v in values)
                        assert new < old
                        old = new
                    assert all(v == average for v in values)
                    assert old == lhs
                    jensen_count += 1
                for D in [1, 2, 6*size]:
                    lhs = factorial(D*m)//prod(factorial(D*j) for j in k)
                    rhs = multinomial(k)**D
                    assert lhs >= rhs
                    stretch_count += 1
    return {"jensen_and_literal_smoothing_cases": jensen_count,
            "stretched_multinomial_cases": stretch_count}


def mod_fraction(x, p):
    assert x.denominator % p != 0, "Normalization denominator vanishes modulo prime"
    return x.numerator*pow(x.denominator, -1, p) % p


def macaulay_rows(forms, size, n, p):
    orders = list(range(n, n+size))
    L = sum(m-1 for m in orders)+1
    columns = list(compositions(L, size))
    lookup = {k: j for j, k in enumerate(columns)}
    rows, labels = [], []
    for m in orders:
        for multiplier in compositions(L-m, size):
            row = [0]*len(columns)
            for k, value in forms[m].items():
                exponent = tuple(a+b for a, b in zip(k, multiplier))
                row[lookup[exponent]] = mod_fraction(value, p)
            rows.append(row)
            labels.append({"moment_order": m, "multiplier_exponents": list(multiplier)})
    return L, columns, rows, labels


def select_independent_rows(rows, p):
    basis = {}
    selected, pivots = [], []
    width = len(rows[0])
    for index, row in enumerate(rows):
        v = list(row)
        for pivot in sorted(basis):
            factor = v[pivot]
            if factor:
                b = basis[pivot]
                for j in range(pivot, width):
                    v[j] = (v[j]-factor*b[j]) % p
        pivot = next((j for j, value in enumerate(v) if value), None)
        if pivot is None:
            continue
        inv = pow(v[pivot], -1, p)
        basis[pivot] = [(x*inv) % p for x in v]
        selected.append(index)
        pivots.append(pivot)
        if len(selected) == width:
            break
    return selected, pivots


def determinant_mod(matrix, p):
    """Independent square-minor elimination, with explicit row swaps."""
    n = len(matrix)
    assert all(len(row) == n for row in matrix)
    a = [list(row) for row in matrix]
    det = 1
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            a[pivot], a[j] = a[j], a[pivot]
            det = -det % p
        value = a[j][j]
        det = det*value % p
        inv = pow(value, -1, p)
        for i in range(j+1, n):
            if a[i][j]:
                factor = a[i][j]*inv % p
                for k in range(j+1, n):
                    a[i][k] = (a[i][k]-factor*a[j][k]) % p
                a[i][j] = 0
    return det


def macaulay_certificate(name, matrix, D, forms, n, p):
    size = len(matrix[0])
    L, columns, rows, labels = macaulay_rows(forms, size, n, p)
    selected, pivots = select_independent_rows(rows, p)
    assert len(selected) == len(columns), f"Macaulay rank failed for {name}, n={n}"
    minor = [rows[i] for i in selected]
    determinant = determinant_mod(minor, p)
    assert determinant != 0
    # A different row order must change the sign of the same full minor.
    swapped = list(minor)
    swapped[0], swapped[1] = swapped[1], swapped[0]
    assert determinant_mod(swapped, p) == (-determinant) % p
    return {
        "name": name, "active_support_size": size, "base_matrix": matrix,
        "private_exponent": D, "starting_order": n,
        "orders": list(range(n, n+size)), "homogeneous_degree_L": L,
        "prime": p, "row_count": len(rows), "column_count": len(columns),
        "full_column_rank": len(selected), "column_monomials": [list(k) for k in columns],
        "row_order": "moment_order ascending; multiplier compositions lexicographic",
        "selected_row_indices_zero_based": selected,
        "selected_row_labels": [labels[i] for i in selected],
        "elimination_pivot_columns": pivots,
        "full_minor_determinant_mod_prime": determinant,
        "matrix_sha256": sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest(),
        "normalized_moment_coefficients": {str(m): form_receipt(forms[m])
                                            for m in range(n, n+size)},
        "certificate_implication":
            "All degree-L coefficient monomials lie in the rational homogeneous ideal; no nonzero projective common zero.",
    }


def negative_controls():
    unbalanced = [[0, 1]]
    try:
        require_balanced(unbalanced)
    except AssertionError:
        pass
    else:
        raise AssertionError("Unbalanced profile was accepted")
    ratios = [F(factorial(m)) for m in [1, 2, 3]]
    assert ratios[1] != ratios[0]**2
    equal_degree_ratios = [F(factorial(2*m), factorial(m)**2) for m in [1, 2, 3]]
    assert equal_degree_ratios[1] != equal_degree_ratios[0]**2
    D = 12
    separate = F(factorial(D)**2, factorial(2*D))
    bound_value = F(1, 2**(D-1))
    merged = F(1)
    assert separate <= bound_value < merged
    # These homogeneous forms genuinely vanish on the nonzero point (0,1).
    deficient = {1: {(1, 0): F(1)}, 2: {(2, 0): F(1)}}
    _, columns, rows, _ = macaulay_rows(deficient, 2, 1, 104729)
    selected, _ = select_independent_rows(rows, 104729)
    assert len(selected) == 2 < len(columns) == 3
    return {
        "unbalanced_1_vs_Y_endpoint_ratios_m1_m2_m3": [frac(x) for x in ratios],
        "equal_private_degree_X2_vs_YZ_endpoint_ratios": [frac(x) for x in equal_degree_ratios],
        "merged_marker_control": {
            "phase_before": "Y*X1^12 + Z*X2^12",
            "phase_after": "X^12*(Y+Z)",
            "monomials_remain_distinct": True,
            "degree_two_mixed_coefficient_before": frac(separate),
            "claimed_private_bound": frac(bound_value),
            "degree_two_mixed_coefficient_after": frac(merged),
            "scope": "failure of coefficient suppression, not an SFC counterexample",
        },
        "deficient_macaulay_system_rank": len(selected),
        "deficient_macaulay_system_columns": len(columns),
        "moment_zero_scope": "P0=1; no false identification with the s-term power sum p0=s is used",
    }


def polynomial_template_controls():
    templates = [
        ("positive_linear_permutations",
         [{(0, 0): F(1), (1, 0): F(1), (0, 1): F(2)},
          {(0, 0): F(1), (1, 0): F(2), (0, 1): F(1)}]),
        ("nonnegative_with_negative_coefficients",
         [{(2, 0): F(1), (1, 1): F(-2), (0, 2): F(1), (1, 0): F(1)},
          {(2, 0): F(1), (1, 1): F(-2), (0, 2): F(1), (0, 1): F(1)}]),
    ]
    out = []
    for name, bases in templates:
        assert {(b, a): value for (a, b), value in bases[0].items()} == bases[1]
        forms, controls = normalized_forms(bases, 12, 3)
        out.append({"name": name, "private_exponent": 12,
                    "positivity_reason":
                        "1+Y+2Z and its permutation" if name.startswith("positive")
                        else "(Y-Z)^2+Y and (Y-Z)^2+Z on the positive orthant",
                    "base_polynomials": [[{"exponents": list(k), **frac(v)}
                                          for k, v in sorted(base.items())] for base in bases],
                    "controls": controls,
                    "normalized_forms": {str(m): form_receipt(forms[m]) for m in range(1, 4)}})
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path)
    parser.add_argument("--note", type=Path, help="Override bytes to check against the fixed reviewed note pin")
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repo_root is None:
        roots = [p for p in Path(__file__).resolve().parents
                 if (p/DIRECTORY).is_dir() and (p/"AGENTS.md").is_file()]
        assert roots, "Repository root not found"
        repo = roots[0]
    else:
        repo = args.repo_root.resolve()
    note = args.note or repo/DIRECTORY/NOTE_NAME
    assert sha256(note.read_bytes()).hexdigest() == NOTE_PIN, "Reviewed note pin mismatch"
    ledger_path = repo/DIRECTORY/LEDGER_NAME
    assert sha256(ledger_path.read_bytes()).hexdigest() == LEDGER_PIN, "Reviewed source ledger pin mismatch"
    ledger = json.loads(ledger_path.read_text())
    assert ledger["source"]["pdf_sha256"] == PDF_PIN
    source = args.source_pdf or repo/ledger["source"]["audit_copy"]
    if args.source_pdf is not None or source.exists():
        assert sha256(source.read_bytes()).hexdigest() == PDF_PIN, "Primary source PDF pin mismatch"
        source_status = {"bytes_rechecked": True, "sha256": PDF_PIN}
    else:
        source_status = {"bytes_rechecked": False, "sha256": PDF_PIN,
                         "scope": "Inspected source pin retained in checked ledger; local audit PDF absent"}
    prime = 104729
    assert all(prime % d for d in range(2, isqrt(prime)+1))
    original = [[0, 1, 3], [2, 0, 1]]
    balanced, provenance = balance(original)
    specifications = [
        ("two_shared_permutations", [[1, 2], [2, 1]], [1, 2, 3, 4]),
        ("three_cyclic_columns", [[2, 0, 1], [1, 2, 0], [0, 1, 2]], [1, 2, 3]),
        ("four_permutation_columns",
         [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]], [1]),
        ("privately_balanced_arbitrary_support", balanced, [1, 2]),
    ]
    certificates, coefficient_controls = [], []
    for name, matrix, starts in specifications:
        profile = require_balanced(matrix)
        size = len(matrix[0])
        D = 6*size
        forms, controls = normalized_forms(monomial_bases(matrix), D, max(starts)+size-1)
        coefficient_controls.append({"name": name, "base_matrix": matrix,
                                     "common_profile": list(profile), "private_exponent": D,
                                     "controls": controls})
        for n in starts:
            certificates.append(macaulay_certificate(name, matrix, D, forms, n, prime))
    # Identical base columns still produce distinct actual monomials via their markers.
    identical = [[2, 2], [1, 1]]
    _, identical_controls = normalized_forms(monomial_bases(identical), 12, 3)
    receipt = {
        "status": "PASS", "schema": 1,
        "scope": "Exact finite transcription controls and listed full-minor certificates; the unrestricted theorem is the pinned proof plus the inspected Turan input.",
        "reviewed_note_sha256": NOTE_PIN, "source_ledger_sha256": LEDGER_PIN,
        "primary_source": source_status,
        "factorial_inequalities": factorial_inequalities(),
        "rational_constants": rational_constant_checks(),
        "balancing": {"original_matrix": original, "balanced_matrix": balanced,
                      "private_copy_provenance": provenance,
                      "common_positive_exponent_profile": list(require_balanced(balanced))},
        "actual_coefficient_controls": coefficient_controls,
        "identical_base_column_controls": identical_controls,
        "polynomial_template_controls": polynomial_template_controls(),
        "macaulay_certificates": certificates,
        "negative_controls": negative_controls(),
        "field_scope": "Finite nonzero modular minors imply rational full rank; rational equations transfer to characteristic-zero coefficient fields.",
        "limits": ["Private markers remain private.",
                   "Endpoint profiles, or polynomial-template endpoint moments, are balanced.",
                   "Amplification is proved; de-amplification and unrestricted SFC are not."],
    }
    receipt["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    output = args.output or repo/DIRECTORY/"sfc_balanced_private_amplification_verification.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True)+"\n")
    print(f"PASS: {len(certificates)} full-rank Macaulay minors; "
          f"largest {max(c['column_count'] for c in certificates)} columns; "
          "exact coefficients, factorial bounds, 6s constants, balancing and negative controls.")
    print(output)


if __name__ == "__main__":
    main()
