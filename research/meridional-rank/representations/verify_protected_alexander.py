#!/usr/bin/env python3
"""Check the integral protected Alexander relations and scan two-block closures."""

from __future__ import annotations

import argparse
from functools import reduce
from hashlib import sha256
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import platform

import sympy as sp

from verify_multiple_protected_obstructions import (
    U,
    cycles,
    hurwitz,
    prefix,
    strand_permutation,
    tail,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
STRANDS = range(8, 13)
EXPONENT_PARAMETERS = range(16)
PRIME_LIMIT = 101


def fox_matrix(word, size, t, prime=None):
    """Propagate the abelianized Fox rows before subtracting the identity."""
    inverse = 1 / t if prime is None else pow(t, -1, prime)
    rows = [[int(i == j) for j in range(size)] for i in range(size)]
    for crossing in word:
        i = abs(crossing) - 1
        a, c = rows[i:i + 2]
        if crossing > 0:
            first = c
            second = [inverse * x + (1 - inverse) * y for x, y in zip(a, c)]
        else:
            first = [(1 - t) * x + t * y for x, y in zip(a, c)]
            second = a
        if prime is None:
            first = [sp.expand(x) for x in first]
            second = [sp.expand(x) for x in second]
        else:
            first = [x % prime for x in first]
            second = [x % prime for x in second]
        rows[i:i + 2] = first, second
    return rows


def direct_fox(word, size, t):
    """Differentiate the literal free words using their prefix exponents."""
    rows = []
    for image in hurwitz([(i,) for i in range(1, size + 1)], word):
        row = [0] * size
        exponent = 0
        for letter in image:
            if letter > 0:
                row[letter - 1] += t**exponent
                exponent += 1
            else:
                exponent -= 1
                row[-letter - 1] -= t**exponent
        assert exponent == 1
        rows.append([sp.expand(value) for value in row])
    return rows


def closure_matrix(word, size, alpha, prime):
    rows = fox_matrix(word, size, alpha, prime)
    return [[(value - int(i == j)) % prime for j, value in enumerate(row)]
            for i, row in enumerate(rows)]


def integral_checks():
    t = sp.Symbol("t")
    protected = sp.Matrix(fox_matrix(U, 4, t)) - sp.eye(4)
    assert protected + sp.eye(4) == sp.Matrix(direct_fox(U, 4, t))
    assert protected * sp.ones(4, 1) == sp.zeros(4, 1)
    f = t**6 - 3*t**5 + 4*t**4 - 7*t**3 + 7*t**2 - 4*t + 1
    h = t**3 - t**2 + 2*t - 1
    assert sp.expand(t**4 * protected[0, 0]) == f
    assert sp.expand(t * protected[0, 3] + (t - 1)**2*h) == 0
    assert f.subs(t, 1) == -1
    assert sp.expand(f - (t**3 - 2*t**2 - 2)*h - (3*t**2 - 1)) == 0
    assert sp.expand(3*h - (t - 1)*(3*t**2 - 1) - (7*t - 4)) == 0
    assert sp.expand(3*(7*t + 4)*(7*t - 4) - 49*(3*t**2 - 1)) == 1

    minors = {}
    for i, j in combinations(range(4), 2):
        value = protected.extract([0, 3], [i, j]).det()
        numerator, _ = sp.fraction(sp.cancel(value))
        minors[i + 1, j + 1] = sp.Poly(numerator, t, domain=sp.ZZ)
    resultant = int(minors[1, 3].resultant(minors[2, 3]))
    assert resultant == -83
    common_factor = reduce(sp.gcd, (
        sp.Poly(p.as_expr(), t, modulus=83) for p in minors.values()
    ))
    assert common_factor.as_expr() == t + 24

    base_word = prefix(2) + tail(2, 8)
    assert fox_matrix(base_word, 8, t) == direct_fox(base_word, 8, t)
    assert module_dimension(closure_matrix([1, 1, 1], 2, 2, 3), 3) == 1
    assert module_dimension(closure_matrix([1, 1, 1], 2, 1, 3), 3) == 0
    return {
        "first_protected_row_entries_1_and_4": [str(protected[0, j]) for j in (0, 3)],
        "integral_unit_ideal_checked": True,
        "minor_numerator_resultant_columns_13_and_23": resultant,
        "all_minor_numerator_gcd_mod_83": str(common_factor.as_expr()),
        "exceptional_local_rank_one": {"prime": 83, "alpha": 59},
        "independent_free_word_controls": ["U", "K_(2,8,0)"],
        "trefoil_t_equals_1_and_t_equals_minus_1_controls": "passed",
    }


def multiply(left, right, prime):
    columns = list(zip(*right))
    return [[sum(a * b for a, b in zip(row, col)) % prime for col in columns]
            for row in left]


def power(matrix, exponent, prime):
    result = [[int(i == j) for j in range(len(matrix))] for i in range(len(matrix))]
    while exponent:
        if exponent % 2:
            result = multiply(result, matrix, prime)
        exponent //= 2
        if exponent:
            matrix = multiply(matrix, matrix, prime)
    return result


def module_dimension(relations, prime):
    # Rows of the full closure matrix sum to zero. Dropping the last
    # coordinate identifies that hyperplane with the module's generators.
    rows = [row[:-1] for row in relations]
    rank = 0
    for j in range(len(rows) - 1):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][j], -1, prime)
        for i in range(rank + 1, len(rows)):
            if rows[i][j]:
                factor = rows[i][j] * inverse % prime
                for c in range(j, len(rows) - 1):
                    rows[i][c] = (rows[i][c] - factor * rows[rank][c]) % prime
        rank += 1
    return len(rows) - 1 - rank


def scan():
    summaries = []
    primes = [p for p in range(2, PRIME_LIMIT + 1)
              if all(p % d for d in range(2, isqrt(p) + 1))]
    for b in STRANDS:
        block, end = prefix(2), tail(2, b)
        mobile_count = b - 4
        assert all(len(cycles(strand_permutation(
            b, block + end * (1 + mobile_count * ell)
        ))) == 1 for ell in EXPONENT_PARAMETERS)
        best = {"dimension": -1}
        by_ell = [-1] * len(EXPONENT_PARAMETERS)
        exceptional = []
        for p in primes:
            for alpha in range(1, p):
                prefix_rows = fox_matrix(block, b, alpha, p)
                tail_rows = fox_matrix(end, b, alpha, p)
                step = power(tail_rows, mobile_count, p)
                matrix = multiply(tail_rows, prefix_rows, p)
                is_exceptional = (p, alpha) == (83, 59)
                local_rank = 1 if is_exceptional else 2
                for ell in EXPONENT_PARAMETERS:
                    relations = [[(value - int(i == j)) % p
                                  for j, value in enumerate(row)]
                                 for i, row in enumerate(matrix)]
                    assert all(sum(row) % p == 0 for row in relations)
                    dimension = module_dimension(relations, p)
                    assert dimension <= b - 1 - 2 * local_rank
                    if alpha == 1:
                        assert dimension == 0
                    if is_exceptional:
                        exceptional.append(dimension)
                        literal = closure_matrix(
                            block + end * (1 + mobile_count * ell), b, alpha, p
                        )
                        assert relations == literal
                    by_ell[ell] = max(by_ell[ell], dimension)
                    if dimension > best["dimension"]:
                        best = {"ell": ell, "prime": p, "alpha": alpha,
                                "dimension": dimension}
                    if ell < EXPONENT_PARAMETERS[-1]:
                        matrix = multiply(step, matrix, p)
        p, alpha, ell = best["prime"], best["alpha"], best["ell"]
        assert module_dimension(closure_matrix(
            block + end * (1 + mobile_count * ell), b, alpha, p
        ), p) == best["dimension"]
        summaries.append({
            "strands": b,
            "maximum": best,
            "maximum_dimension_for_each_ell": by_ell,
            "dimensions_at_exceptional_pair_for_each_ell": exceptional,
            "best_meridional_lower_bound": best["dimension"] + 1,
        })
    return {
        "strands": list(STRANDS),
        "ell": list(EXPONENT_PARAMETERS),
        "primes": primes,
        "nonzero_parameters_tested_per_closure_word": sum(p - 1 for p in primes),
        "closure_words_checked": len(STRANDS) * len(EXPONENT_PARAMETERS),
        "specializations_checked": len(STRANDS) * len(EXPONENT_PARAMETERS)
        * sum(p - 1 for p in primes),
        "summaries": summaries,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    result = {
        "status": "passed",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "sources": {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                    for path in (Path(__file__).resolve(),
                                 HERE / "verify_multiple_protected_obstructions.py")},
        "integral_controls": integral_checks(),
        "prime_field_search": scan(),
    }
    if args.write_receipt:
        args.write_receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "specializations_checked": result["prime_field_search"]["specializations_checked"],
        "maxima": [entry["maximum"] | {"strands": entry["strands"]}
                   for entry in result["prime_field_search"]["summaries"]],
    }, indent=2))


if __name__ == "__main__":
    main()
