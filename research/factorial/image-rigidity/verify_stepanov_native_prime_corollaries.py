#!/usr/bin/env python3
"""Exact guards for STEPANOV_NATIVE_PRIME_OFFSET_COROLLARIES.md.

Universal conclusions follow from the written algebra and declared inputs,
not from extrapolation of these finite controls.
"""

import argparse
import hashlib
import json
from fractions import Fraction
from math import comb, factorial, isqrt, prod
from pathlib import Path

import sympy as sp

if not __debug__:
    raise RuntimeError("Assertions must execute; do not run with -O or -OO")

SOURCE_PINS = {
    "DIAGONAL_PRIME_GAP_SFC.tex": "db034b5fb29ece9c3a0cb8cbce4c82f692b933e3679b27ee823db5ca5e7f98eb",
    "LOWER_TWO_PRIME_GAP_SFC.tex": "8453b5b0ccd220573ef106fa26790e1f691885970b8c2116f036ed893ab2c200",
    "HIGHER_GAP_ADJOINT_RANK.tex": "ec19415938c9b0b92e7d365fb4a4c66b701fdf36d73883b68f7cfb2971d7ec28",
    "PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex": "00b5561db7b6f42dcf0c67a3f848e1243ae2b70e8e2ba976a6c0eed4cdc0a169",
    "RELATIVE_GAP_THREE_STRONG_FACTORIAL.md": "0e9ab6b15400c7e1b6a7d5630fc781d6bfab9f20483cfb76d7476fc8667bda2e",
}


def check_sources(repo):
    folder = repo / 'research/factorial/image-rigidity'
    actual = {name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
              for name in SOURCE_PINS}
    if actual != SOURCE_PINS:
        raise ValueError('Changed pinned mathematical input')
    return actual


def valuation(number, prime):
    number = abs(number)
    assert number != 0
    exponent = 0
    while number % prime == 0:
        exponent += 1
        number //= prime
    return exponent


def factor(number):
    number = abs(number)
    assert number != 0
    factors = {}
    divisor = 2
    while divisor * divisor <= number:
        while number % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            number //= divisor
        divisor += 1
    if number > 1:
        factors[number] = factors.get(number, 0) + 1
    return factors


def prime(number):
    return number >= 2 and all(number % d for d in range(2, isqrt(number) + 1))


def harmonic(number):
    return sum((Fraction(1, j) for j in range(1, number + 1)), Fraction())


def hasse_values(a):
    x, y = sp.symbols("x y")
    sigma = sp.Poly(sp.expand(prod(x + y + j for j in range(a))), y).coeff_monomial(y ** 2)
    values = []
    for j in range(a):
        direct = int(sigma.subs(x, -j))
        harmonic_value = ((-1) ** j * factorial(j) * factorial(a - 1 - j)
                          * (harmonic(a - 1 - j) - harmonic(j)))
        assert harmonic_value.denominator == 1
        assert direct == harmonic_value
        assert (direct == 0) == (a % 2 == 1 and j == (a - 1) // 2)
        values.append(direct)
    return sigma, values


def apply_u(vector, p, shift):
    result = [(shift - j) * vector[j] % p for j in range(p)]
    for j in range(1, p):
        result[j] = (result[j] + vector[j - 1]) % p
    return result


def modular_rank(columns, p):
    rows = [list(row) for row in zip(*columns)]
    rank = 0
    for column in range(len(columns)):
        pivot = next((row for row in range(rank, len(rows))
                      if rows[row][column] % p), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column] % p, -1, p)
        rows[rank] = [entry * inverse % p for entry in rows[rank]]
        for row in range(len(rows)):
            if row != rank:
                multiple = rows[row][column] % p
                rows[row] = [(left - multiple * right) % p
                             for left, right in zip(rows[row], rows[rank])]
        rank += 1
    return rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name(
        'stepanov_native_prime_corollaries_verification.json'))
    args = parser.parse_args()
    source_hashes = check_sources(args.repo_root)
    a, p, n = sp.symbols("a p n")
    c = (a + p) * n + p + 1
    assert sp.expand(c + a - 1 - (a + p) * (n + 1)) == 0
    checks = {"symbolic_last_factor_cancellation": True}

    rows = []
    all_values = {}
    all_sigmas = {}
    for offset in range(1, 15):
        sigma, values = hasse_values(offset)
        all_values[offset] = values
        all_sigmas[offset] = sigma
        factorizations = []
        exceptional = set()
        for value in values:
            if value == 0:
                factorizations.append(None)
                continue
            factors = factor(value)
            assert prod(base ** exponent for base, exponent in factors.items()) == abs(value)
            assert all(prime(base) for base in factors)
            exceptional.update(base for base in factors if base > offset)
            factorizations.append(factors)
        rows.append({"a": offset, "sigma": str(sigma), "values": values,
                     "factorizations": factorizations,
                     "exceptional_primes_above_a": sorted(exceptional)})
    checks["hasse_harmonic_certificates"] = rows
    expected = {3: [], 4: [11], 5: [], 6: [13, 137], 7: [11],
                8: [11, 29, 47], 10: [13, 17, 37, 7129]}
    for offset, exceptions in expected.items():
        assert rows[offset - 1]["exceptional_primes_above_a"] == exceptions

    residue_records = []
    valuation_guards = 0
    for gap in [2, 3, 5, 7, 11, 13]:
        for offset in range(1, gap):
            inverse = pow(offset + gap, -1, gap * gap)
            predicted = sorted({r * inverse % (gap * gap)
                                for r in range(1, offset)})
            actual = []
            refined = []
            for terminal_residue in range(gap * gap):
                # This positive representative has N>=2, including residue 0,1.
                terminal = terminal_residue + gap * gap
                lower = terminal - 1
                c_value = (offset + gap) * lower + gap + 1
                original = Fraction(gap * (lower + 1),
                                    prod(c_value + j for j in range(offset)))
                canceled = Fraction(gap,
                                    (offset + gap) * prod(c_value + j for j in range(offset - 1)))
                assert original == canceled
                val = valuation(original.numerator, gap) - valuation(original.denominator, gap)
                deep_indices = [j for j in range(offset - 1)
                                if (c_value + j) % (gap * gap) == 0]
                assert (val < 0) == bool(deep_indices)
                assert len(deep_indices) <= 1
                if val < 0:
                    actual.append(terminal_residue)
                    j = deep_indices[0]
                    assert (offset * lower + 1) % gap == (-j) % gap
                    if all_values[offset][j] % gap == 0:
                        refined.append(terminal_residue)
                valuation_guards += 1
            assert actual == predicted
            assert len(actual) == offset - 1
            if offset in [3, 6, 8, 10]:
                residue_records.append({"a": offset, "p": gap,
                                        "modulus": gap * gap,
                                        "integral_branch_residual_N": actual,
                                        "refined_residual_N": refined})
    checks["valuation_congruence_guards"] = valuation_guards
    checks["selected_residue_certificates"] = residue_records
    requested = {(row["a"], row["p"]): row for row in residue_records}
    assert requested[3, 5]["integral_branch_residual_N"] == [19, 22]
    assert requested[3, 5]["refined_residual_N"] == [22]
    assert requested[3, 7]["integral_branch_residual_N"] == [5, 10]
    assert requested[3, 7]["refined_residual_N"] == [5]

    # Every shift and every positive power class modulo p-1 for these primes.
    krylov_records = []
    for gap in [2, 3, 5, 7, 11, 13]:
        count = 0
        for shift in range(gap):
            powers = [[1] + [0] * (gap - 1)]
            for exponent in range(1, 3 * gap):
                powers.append(apply_u(powers[-1], gap, shift))
            assert modular_rank(powers[:gap], gap) == gap
            assert powers[gap] == powers[1]
            for exponent in range(1, gap):
                assert modular_rank(powers[exponent:exponent + gap - 1], gap) == gap - 1
                assert powers[exponent + gap - 1] == powers[exponent]
                count += 1
        krylov_records.append({"p": gap, "all_shifts": gap,
                               "all_positive_power_classes": gap - 1,
                               "rank_checks": count})
    checks["krylov_certificates"] = krylov_records

    # Compare actual integer recurrence with the reduced operator, rather than
    # constructing both sides with the same modular update.
    z = sp.Symbol("z")
    recurrence_cases = 0
    for gap in [2, 3, 5, 7]:
        for offset in range(1, gap):
            for lower in [1, 2, 3, 4]:
                gamma = offset * lower + 1
                polynomial = sp.Integer(1)
                vector = [1] + [0] * (gap - 1)
                for index in range(lower + gap):
                    coefficients = [int(polynomial.expand().coeff(z, j)) % gap
                                    for j in range(gap)]
                    assert coefficients == vector
                    polynomial = sp.expand((z + gap * index + gamma + gap - 1) * polynomial
                                           - z * sp.diff(polynomial, z))
                    vector = apply_u(vector, gap, offset * lower % gap)
                recurrence_cases += 1
    checks["actual_auxiliary_recurrence_cases"] = recurrence_cases

    # Direct coefficient-valuation guards for the lower reciprocal.
    lower_bounds = 0
    for gap in [2, 3, 5, 7]:
        for offset in range(1, gap):
            for lower in range(1, 7):
                gamma = offset * lower + 1
                for k in range(1, lower + 1):
                    coefficient = comb(lower, k) * prod(gamma + j for j in range(gap * k))
                    assert valuation(coefficient, gap) >= k
                    lower_bounds += 1
    checks["lower_coefficient_valuation_guards"] = lower_bounds

    # Kernel/Hasse evaluation used by the added digit. Constant rising-basis
    # coefficient is omitted exactly as in the actual R=0 jet equation.
    x = sp.Symbol("x")
    kernel_evaluations = 0
    for gap in [2, 3, 5, 7, 11, 13]:
        derivatives = [sp.Poly(sp.diff(prod(x + j for j in range(m)), x), x)
                       for m in range(1, gap + 1)]
        for residue in range(gap):
            kernel = [int(poly.eval(residue)) % gap for poly in derivatives]
            assert kernel[-1] == gap - 1
            for offset in range(1, gap + 2):
                row = [prod(residue + m + 1 + j for j in range(offset - 1 - m)) % gap
                       if m < offset else 0 for m in range(1, gap + 1)]
                contraction = sum(left * right for left, right in zip(row, kernel)) % gap
                assert contraction == int(all_sigmas[offset].subs(x, residue)) % gap
                kernel_evaluations += 1
    checks["primitive_kernel_hasse_evaluations"] = kernel_evaluations

    # The new p^3 classes are obtained by exact inversion. Checking every lift
    # of each former p^2 class also certifies the genuine extra exclusion:
    # p-1 lifts have kappa valuation -1 and just one survives this lemma.
    cubic_records = []
    expected_cubic = {(3, 5): [47], (3, 7): [103], (4, 11): [1065],
                      (6, 13): [694, 1272], (8, 11): [841],
                      (10, 13): [669, 2006], (10, 17): [364, 1274]}
    eliminated_lift_count = 0
    for (offset, gap), expected_classes in expected_cubic.items():
        indices = [j for j in range(offset - 1) if all_values[offset][j] % gap == 0]
        inverse2 = pow(offset + gap, -1, gap ** 2)
        inverse3 = pow(offset + gap, -1, gap ** 3)
        former = sorted((offset - 1 - j) * inverse2 % gap ** 2 for j in indices)
        refined = sorted((offset - 1 - j) * inverse3 % gap ** 3 for j in indices)
        assert refined == expected_classes
        eliminated = []
        for former_class in former:
            lifts = [former_class + j * gap ** 2 for j in range(gap)]
            assert sum(terminal in refined for terminal in lifts) == 1
            for terminal in lifts:
                assert terminal >= 2
                c_value = (offset + gap) * terminal - offset + 1
                kappa = Fraction(gap, (offset + gap)
                                 * prod(c_value + j for j in range(offset - 1)))
                val = valuation(kappa.numerator, gap) - valuation(kappa.denominator, gap)
                sigma_residue = int(all_sigmas[offset].subs(x, offset * (terminal - 1) + 1)) % gap
                assert sigma_residue == 0
                if terminal in refined:
                    assert val <= -2
                else:
                    assert val == -1
                    eliminated.append(terminal)
        eliminated_lift_count += len(eliminated)
        cubic_records.append({"a": offset, "p": gap,
                              "hasse_zero_indices": indices,
                              "former_modulus": gap ** 2, "former_classes": former,
                              "refined_modulus": gap ** 3, "refined_classes": refined,
                              "eliminated_lifts_mod_p3": sorted(eliminated)})
    checks["cubic_modulus_certificates"] = cubic_records
    checks["boundary_minus_one_eliminated_lifts"] = eliminated_lift_count
    checks["first_residual_after_digit"] = {"a3_p5": 47, "a3_p7": 103,
                                            "requires_finite_prefix": False}

    # Positive rational p-integral shapes have the same lower-root valuation
    # premise. These finite coefficient guards support its written block proof.
    rational_guards = 0
    for gamma, gaps in [(Fraction(1, 2), [3, 5, 7]),
                        (Fraction(7, 3), [2, 5, 7]),
                        (Fraction(11, 5), [2, 3, 7])]:
        for gap in gaps:
            assert gamma.denominator % gap
            for lower in range(1, 5):
                for k in range(1, lower + 1):
                    coefficient = comb(lower, k) * prod(gamma + j for j in range(gap * k))
                    assert valuation(coefficient.numerator, gap) - valuation(coefficient.denominator, gap) >= k
                    rational_guards += 1
    checks["rational_shape_coefficient_valuation_guards"] = rational_guards

    # Negative scope control: mere ramified integrality does not imply the
    # lower reciprocal-root valuation >=1 used in the new argument.
    square_root = sp.sqrt(3)
    assert sp.expand((1 + square_root) * (2 + square_root) * (3 + square_root)) == 24 + 14 * square_root
    term_valuations = [Fraction(valuation(24, 3)),
                       Fraction(valuation(14, 3)) + Fraction(1, 2)]
    assert term_valuations == [Fraction(1), Fraction(1, 2)]
    assert min(term_valuations) < 1
    checks["ramified_scope_negative_control"] = {
        "p": 3, "n": 1, "gamma": "1+sqrt(3)",
        "minimal_polynomial_of_sqrt3": "s^2-3",
        "rising_factorial": "24+14*sqrt(3)",
        "summand_valuations": [str(v) for v in term_valuations],
        "lower_reciprocal_root": "-24-14*sqrt(3)",
        "root_valuation": "1/2",
        "refutes": "v_p(u)>=1 from mere ramified integrality; not SFC",
    }

    assert source_hashes == check_sources(args.repo_root)
    receipt = {"status": "PASS", "sympy_version": sp.__version__,
               "scope": "finite exact guards and constants for retained-theorem corollaries plus the written ordinary extra-digit refinement; no unrestricted SFC conclusion",
               "source_sha256": source_hashes,
               "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "checks": checks}
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({"status": "PASS", "output": str(args.output),
                      "check_groups": list(checks)}))


if __name__ == "__main__":
    main()
