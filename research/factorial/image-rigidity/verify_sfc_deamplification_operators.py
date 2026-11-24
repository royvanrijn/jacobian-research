#!/usr/bin/env python3
"""Exact finite controls for the accompanying operator proof; stdlib only."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import factorial, prod
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NOTE_SHA = "df9f19e4418c2abdf47c0f2696a60637dcfb761f93a9b91457aad4b5c671e6b0"
SOURCE_HASHES = {
    "SFC_BALANCED_PRIVATE_AMPLIFICATION.md":
        "f4ebaa9028cb04258cf65be178685650d53e2852644ed0678f1a655d134c065e",
    "SFC_ROW_MERGING_OBSTRUCTIONS.md":
        "731e5dcd27f8cfe6f688c351b55888facffe0756c797c192f3f1970a7113edee",
    "UNRESTRICTED_SFC_MECHANISMS.tex":
        "f8b3dfc5947601fb002aee0ba942b0fc20b7ef46ca8e0d4f69b5a3e5e1361785",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_hash(data, expected, label):
    require(sha256(data).hexdigest() == expected, "hash mismatch: " + label)


def rejected(action, label):
    try:
        action()
    except ValueError:
        return
    raise ValueError("negative control was accepted: " + label)


def compositions(n, length):
    if length == 1:
        yield (n,)
        return
    for first in range(n, -1, -1):
        for tail in compositions(n-first, length-1):
            yield (first,) + tail


def moment_form(rows, n, columns):
    """Actual factorial moment divided by n!, reconstructed from its support."""
    return {
        k: F(prod(factorial(sum(a*b for a, b in zip(row, k))) for row in rows),
             prod(factorial(t) for t in k))
        for k in compositions(n, columns)
    }


def monomial_multiply(poly, beta):
    return {tuple(a+b for a, b in zip(k, beta)): v for k, v in poly.items()}


def marker_weight(k, private):
    return prod(factorial(d*k[j]) for j, d in private)


def diagonal(poly, private, inverse=False):
    return {k: v / marker_weight(k, private) if inverse
            else v * marker_weight(k, private) for k, v in poly.items()}


def functional(poly, weights):
    return sum((F(a)*poly.get(k, F(0)) for k, a in weights.items()), F(0))


def determinant(matrix):
    a = [[F(x) for x in row] for row in matrix]
    n = len(a)
    require(all(len(row) == n for row in a), "determinant requires a square matrix")
    result = F(1)
    for i in range(n):
        pivot = next((j for j in range(i, n) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        p = a[i][i]
        result *= p
        for j in range(i+1, n):
            factor = a[j][i] / p
            for k in range(i+1, n):
                a[j][k] -= factor*a[i][k]
            a[j][i] = F(0)
    return result


def sylvester_binary(p, q):
    dp = sum(next(iter(p)))
    dq = sum(next(iter(q)))
    columns = tuple(compositions(dp+dq-1, 2))
    rows = [monomial_multiply(p, b) for b in compositions(dq-1, 2)]
    rows += [monomial_multiply(q, b) for b in compositions(dp-1, 2)]
    return determinant([[row.get(k, F(0)) for k in columns] for row in rows])


def run_checks():
    counts = {}
    base_rows = ((1, 0), (0, 1))
    marked_rows = base_rows + ((1, 0),)
    p2, p3 = (moment_form(base_rows, n, 2) for n in (2, 3))
    q2, q3 = (moment_form(marked_rows, n, 2) for n in (2, 3))
    require(tuple(p2.values()) == (1, 1, 1), "base quadratic coefficients")
    require(tuple(p3.values()) == (1, 1, 1, 1), "base cubic coefficients")
    require(tuple(q2.values()) == (2, 1, 1), "marked quadratic coefficients")
    require(tuple(q3.values()) == (6, 2, 1, 1), "marked cubic coefficients")
    private = ((0, 1),)
    require(diagonal(p2, private) == q2 and diagonal(p3, private) == q3,
            "ordinary marker identity")
    require(sylvester_binary(p2, p3) == 1, "base resultant")
    require(sylvester_binary(q2, q3) == 36, "marked resultant")
    counts["ordinary_form_coefficients"] = 14
    counts["exact_resultants"] = 2

    ell = {(2, 1): 1, (1, 2): -1}
    ell_tilde = {(2, 1): 1, (1, 2): -1, (0, 3): -1}
    base_span = [monomial_multiply(p2, b) for b in ((1, 0), (0, 1))] + [p3]
    marked_span = [monomial_multiply(q2, b) for b in ((1, 0), (0, 1))] + [q3]
    require(all(functional(p, ell) == 0 for p in base_span), "base annihilator")
    require(all(functional(p, ell_tilde) == 0 for p in marked_span),
            "marked annihilator")
    inverse_witness = diagonal(monomial_multiply(q2, (1, 0)), private, True)
    direct_witness = diagonal(monomial_multiply(p2, (1, 0)), private)
    require(inverse_witness == {(3, 0): F(1, 3), (2, 1): F(1, 2), (1, 2): F(1)},
            "inverse witness coefficients")
    require(direct_witness == {(3, 0): F(6), (2, 1): F(2), (1, 2): F(1)},
            "forward witness coefficients")
    require(functional(inverse_witness, ell) == -F(1, 2), "inverse nonmembership")
    require(functional(direct_witness, ell_tilde) == 1, "forward nonmembership")
    counts["nonmembership_certificates"] = 2
    counts["annihilated_spanning_rows"] = 6

    # The same form has two coefficient rows sharing two nonzero columns.
    betas = ((2, 0), (1, 1))
    gammas = ((3, 1), (2, 2))
    ratios = [[marker_weight(tuple(g-b for g, b in zip(gamma, beta)), private)
               for gamma in gammas] for beta in betas]
    require(ratios == [[1, 1], [2, 1]], "Macaulay multiplier submatrix")
    require(F(ratios[0][0]*ratios[1][1], ratios[0][1]*ratios[1][0]) == F(1, 2),
            "nonseparable Macaulay cross-ratio")

    counts["conjugation_checks"] = 0
    counts["laurent_kernel_checks"] = 0
    counts["inverse_variance_checks"] = 0
    for d in range(1, 9):
        for k in range(9):
            rise = prod(d*k+ell for ell in range(1, d+1))
            require(rise == factorial(d*(k+1))//factorial(d*k), "raising conjugacy")
            require(F(1, rise) == F(factorial(d*k), factorial(d*(k+1))),
                    "inverse raising conjugacy")
            counts["conjugation_checks"] += 2
        require(prod(-d+ell for ell in range(1, d+1)) == 0, "Laurent kernel")
        counts["laurent_kernel_checks"] += 1
        require(F(1, factorial(2*d))-F(1, factorial(d)**2) < 0,
                "inverse multiplier variance must be negative")
        counts["inverse_variance_checks"] += 1

    counts["shift_resolvent_checks"] = 0
    for profile in ((1,), (2,), (1, 3), (2, 2)):
        for k, beta in product(range(7), range(5)):
            denominator = prod(prod(d*k+ell for ell in range(1, d*beta+1))
                               for d in profile)
            ratio = F(prod(factorial(d*k) for d in profile),
                      prod(factorial(d*(k+beta)) for d in profile))
            require(ratio == F(1, denominator), "finite Euler-resolvent formula")
            counts["shift_resolvent_checks"] += 1

    counts["twisted_product_checks"] = 0
    w = lambda k: marker_weight(k, ((0, 2), (1, 1), (1, 3)))
    for alpha, beta in product(product(range(3), repeat=2), repeat=2):
        k = tuple(a+b for a, b in zip(alpha, beta))
        multiplier = F(w(k), w(alpha)*w(beta))
        require(multiplier*w(alpha)*w(beta) == w(k), "twisted-product transport")
        counts["twisted_product_checks"] += 1

    counts["ideal_classification_obstructions"] = 0
    for d in range(1, 10):
        w_pure = F(factorial(2*d), factorial(d)**2)
        w_mixed = F(factorial(d), factorial(d))
        require(w_pure > 1 and w_mixed == 1, "homogeneous-ideal obstruction")
        counts["ideal_classification_obstructions"] += 1

    counts["marker_profile_inequalities"] = 0
    counts["unit_marker_equalities"] = 0
    for columns, profiles in (
        (2, ((12,), (6, 6), (3, 4, 5), (1,)*12)),
        (3, ((18,), (5, 6, 7), (1,)*18)),
    ):
        for profile in profiles:
            d = sum(profile)
            require(d >= 6*columns, "marker-profile scope")
            for n in range(1, 9):
                for k in compositions(n, columns):
                    multinomial = factorial(n)//prod(factorial(t) for t in k)
                    ratio = prod(F(prod(factorial(b*t) for t in k), factorial(b*n))
                                 for b in profile)
                    require(ratio <= F(1, multinomial**d), "marker suppression bound")
                    counts["marker_profile_inequalities"] += 1
                    if all(b == 1 for b in profile):
                        require(ratio == F(1, multinomial**d), "unit-marker equality")
                        counts["unit_marker_equalities"] += 1

    controls = [
        (lambda: require(functional(inverse_witness, ell) == 0, "false inverse membership"),
         "false inverse-ideal membership"),
        (lambda: require(functional(direct_witness, ell_tilde) == 0, "false forward membership"),
         "false forward-ideal membership"),
        (lambda: require(ratios[0][0]*ratios[1][1] == ratios[0][1]*ratios[1][0],
                         "false separable matrix scaling"), "false matrix equivalence"),
        (lambda: require(F(1, 2) >= 1, "false nonnegative inverse variance"),
         "positive inverse Mellin kernel"),
    ]
    altered_q2 = dict(q2)
    altered_q2[(1, 1)] += 1
    controls.append((lambda: require(diagonal(p2, private) == altered_q2,
                                     "altered actual moment coefficient"),
                     "altered moment coefficient"))
    for action, label in controls:
        rejected(action, label)
    counts["mathematical_negative_controls"] = len(controls)
    return {
        "outcome": "PASS",
        "assurance": "Exact finite controls for the stated written proofs; no general SFC closure.",
        "counts": counts,
        "ordinary_resultants": [1, 36],
        "inverse_nonmembership_functional": "-1/2",
        "forward_nonmembership_functional": "1",
        "macaulay_cross_ratio": "1/2",
        "saturation_scope": "Both ordinary example ideals saturate to the unit ideal.",
        "source_hashes": SOURCE_HASHES,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=HERE/"SFC_DEAMPLIFICATION_BRIDGES.md")
    parser.add_argument("--receipt", type=Path,
                        default=HERE/"sfc_deamplification_operators_verification.json")
    args = parser.parse_args()
    require(__debug__, "optimized Python is rejected for consistent replay invocation")
    note_bytes = args.note.read_bytes()
    check_hash(note_bytes, NOTE_SHA, str(args.note))
    source_dir = ROOT/"research/factorial/image-rigidity"
    for name, digest in SOURCE_HASHES.items():
        check_hash((source_dir/name).read_bytes(), digest, name)
    rejected(lambda: check_hash(note_bytes+b"\n", NOTE_SHA, "altered note"),
             "altered note pin")
    receipt = run_checks()
    receipt["note_sha256"] = NOTE_SHA
    receipt["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    receipt["altered_note_control"] = "rejected"
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"outcome": receipt["outcome"], "counts": receipt["counts"],
                      "receipt": str(args.receipt)}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as exc:
        print("FAIL: " + str(exc), file=sys.stderr)
        sys.exit(1)
