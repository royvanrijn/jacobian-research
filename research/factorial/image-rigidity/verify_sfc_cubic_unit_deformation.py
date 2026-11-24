#!/usr/bin/env python3
"""Exact controls for the written uniform determinant and kernel-transfer proof."""

import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, factorial, prod
from pathlib import Path
import platform


if not __debug__:
    raise RuntimeError("This assertion-based checker requires assertions enabled; do not use python -O.")


# Reviewed proof bytes. These must be changed only after reviewing the changed
# statement/dependency; recording whatever bytes happen to be present is not a pin.
EXPECTED_SOURCE_SHA256 = {
    "research/factorial/image-rigidity/SFC_CUBIC_UNIT_DEFORMATION_TRANSFER.md":
        "0dff58f64cd64b3e324e71c17d7314de849eaeccf797e4177414a89c7684b57a",
    "research/factorial/image-rigidity/STEPANOV_ADJACENT_POWERS_OF_THREE_SFC.md":
        "8dbda8fec87afcf02fe7e0ade2a09284c6b080f2212d61f33ee990f428ac804d",
    "research/factorial/image-rigidity/DYADIC_BINOMIAL_SFC.tex":
        "9d4a5a3084c3a2d9956c83c444e216f528efa523d21006e20a0da00bf4011f61",
}


def verify_reviewed_sources(repo_root):
    records = []
    for relative, expected in EXPECTED_SOURCE_SHA256.items():
        actual = sha256((repo_root / relative).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError(f"Reviewed source hash mismatch: {relative}; expected {expected}, got {actual}")
        records.append({"path": relative, "sha256": actual,
                        "expected_sha256": expected, "reviewed_pin_verified": True})
    own_path = "research/factorial/image-rigidity/verify_sfc_cubic_unit_deformation.py"
    records.append({"path": own_path,
                    "sha256": sha256((repo_root / own_path).read_bytes()).hexdigest(),
                    "role": "current checker bytes; external manifests may pin this hash"})
    return records


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a or [0]


def add(a, b):
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def scale(a, c):
    return trim([c * x for x in a])


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, av in enumerate(a):
        for j, bv in enumerate(b):
            result[i + j] += av * bv
    return trim(result)


def power(a, exponent):
    result = [1]
    for _ in range(exponent):
        result = mul(result, a)
    return result


def derivative(a):
    return trim([i * a[i] for i in range(1, len(a))])


def evaluate(a, x):
    result = 0
    for c in reversed(a):
        result = result * x + c
    return result


def residue(value, modulus):
    value = Fraction(value)
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def reduce_poly(a, prime):
    return trim([residue(c, prime) for c in a])


def remainder(a, b, prime):
    a, b = reduce_poly(a, prime), reduce_poly(b, prime)
    inverse = pow(b[-1], -1, prime)
    while a != [0] and len(a) >= len(b):
        shift = len(a) - len(b)
        a = reduce_poly(sub(a, [0] * shift + scale(b, a[-1] * inverse % prime)), prime)
    return a


def gcd_mod(a, b, prime):
    a, b = reduce_poly(a, prime), reduce_poly(b, prime)
    while b != [0]:
        a, b = b, remainder(a, b, prime)
    return reduce_poly(scale(a, pow(a[-1], -1, prime)), prime)


def valuation(value, prime):
    value = Fraction(value)
    assert value
    numerator, denominator = abs(value.numerator), value.denominator
    result = 0
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def resultant(a, b):
    n, m = len(a) - 1, len(b) - 1
    rows = []
    for shift in range(m):
        rows.append([0] * shift + list(reversed(a)) + [0] * (m - shift - 1))
    for shift in range(n):
        rows.append([0] * shift + list(reversed(b)) + [0] * (n - shift - 1))
    matrix = [[Fraction(c) for c in row] for row in rows]
    sign, previous = 1, Fraction(1)
    for k in range(len(matrix) - 1):
        pivot_row = next((i for i in range(k, len(matrix)) if matrix[i][k]), None)
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != k:
            matrix[k], matrix[pivot_row] = matrix[pivot_row], matrix[k]
            sign = -sign
        pivot = matrix[k][k]
        for i in range(k + 1, len(matrix)):
            for j in range(k + 1, len(matrix)):
                matrix[i][j] = (matrix[i][j] * pivot - matrix[i][k] * matrix[k][j]) / previous
        for i in range(k + 1, len(matrix)):
            matrix[i][k] = Fraction(0)
        previous = pivot
    return sign * matrix[-1][-1]


def model(n, prime=3):
    return [(-1) ** k * comb(n + (prime - 1) * k, prime * k)
            for k in range(n, -1, -1)]


@lru_cache(maxsize=None)
def unit_block(prime, m):
    return prod(j for j in range(1, prime * m + 1) if j % prime)


def normalized_kernel(n, units, zeta):
    eta = -2 * Fraction(zeta) ** 2
    result = [Fraction(0)] * (n + 1)
    for k in range(n + 1):
        # Direct factorial expression, independent of the displayed binomial-unit formula.
        result[n - k] = (Fraction(comb(n, k) * factorial(n + 2 * k), factorial(n))
                         * Fraction(units(n + 2 * k), units(n)) / (3 * eta) ** k)
        claimed = (comb(n + 2 * k, 3 * k) * unit_block(3, k)
                   * Fraction(units(n + 2 * k), units(n)) / eta ** k)
        assert result[n - k] == claimed
        assert result[n - k].denominator % 3 != 0
        assert residue(result[n - k], 9) == residue(model(n)[n - k], 9)
    assert result[-1] == 1 and valuation(result[0], 3) == 0
    return result


def model_and_determinant_guards():
    records = []
    for n in [1, 2, 3, 4, 5, 6, 12, 16]:
        f, g = model(n), model(n + 1)
        common = gcd_mod(f, g, 3)
        expected = len(common) - 1
        assert gcd_mod(common, derivative(f), 3) == [1]
        value = resultant(f, g)
        assert valuation(value, 3) == expected
        perturbed_f = add(f, [9 * (-1) ** j * (n + 2) * (j + 1) for j in range(n)])
        perturbed_g = add(g, [9 * (n + 3) * (j + 2) for j in range(n + 1)])
        changed = resultant(perturbed_f, perturbed_g)
        assert valuation(changed, 3) == expected
        records.append({"n": n, "residual_gcd_degree": expected,
                        "model_resultant_v3": valuation(value, 3),
                        "perturbed_resultant_v3": valuation(changed, 3)})
    large = []
    for n, expected in [(39, 10), (52, 13)]:
        common = gcd_mod(model(n), model(n + 1), 3)
        assert len(common) - 1 == expected
        assert gcd_mod(common, derivative(model(n)), 3) == [1]
        large.append({"n": n, "degree": expected, "gcd_ascending_mod3": common})

    ss = [[3], [-1], [1]]
    vv = [[3], [0], [0, 2]]
    for m in range(3, 76):
        ss.append(add(scale(ss[m - 1], -1), mul([0, 1], ss[m - 3])))
    for m in range(3, 26):
        vv.append(add(mul([0, 1], vv[m - 2]), mul([0, 0, 1], vv[m - 3])))
    for n in range(1, 26):
        assert ss[3 * n] == add(model(n), mul([0, 2], model(n - 1)))
        assert ss[3 * n] == add(sub(power(ss[n], 3), scale(mul(ss[n], vv[n]), 3)), [0] * n + [3])
    return {"exact_resultant_cases": records, "growing_residual_degree_controls": large,
            "both_trace_identities_checked_for_N": [1, 25]}


def kernel_guards():
    records = []
    for label, zeta, units in [
        ("constant units", 1, lambda L: 1),
        ("nonconstant unit deformation", 2, lambda L: 2 ** L * (1 + 9 * L * L))]:
        for n in range(1, 9):
            f = normalized_kernel(n, units, zeta)
            g = normalized_kernel(n + 1, units, zeta)
            degree = len(gcd_mod(model(n), model(n + 1), 3)) - 1
            assert valuation(resultant(f, g), 3) == degree
            raw_f = [comb(n, k) * factorial(n + 2 * k) * units(n + 2 * k) for k in range(n + 1)]
            raw_g = [comb(n + 1, k) * factorial(n + 1 + 2 * k) * units(n + 1 + 2 * k)
                     for k in range(n + 2)]
            raw_value = resultant(raw_f, raw_g)
            expected = n * (n + 1) + (n + 1) * valuation(factorial(n), 3) + n * valuation(factorial(n + 1), 3) + degree
            assert valuation(raw_value, 3) == expected
            records.append({"kernel": label, "n": n, "normalized_resultant_v3": degree,
                            "raw_resultant_v3": expected})
    assert normalized_kernel(1, lambda L: 1, 1) == [-1, 1]
    assert normalized_kernel(2, lambda L: 1, 1) == [10, -4, 1]
    return {"cases": records, "r_zero_extra_unit_scaling": {"eta": -2,
            "G1": [-1, 1], "G2": [10, -4, 1], "old_eta_one_coefficient_mod9": 2}}


def complete_homogeneous(weights, maximum):
    result = [Fraction(1)] + [Fraction(0)] * maximum
    for weight in weights:
        for degree in range(1, maximum + 1):
            result[degree] += weight * result[degree - 1]
    return result


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, length - 1):
                yield (first,) + rest


def direct_linear_moment(weights, degree):
    result = Fraction(0)
    for exponent in compositions(degree, len(weights)):
        multinomial = factorial(degree) // prod(factorial(e) for e in exponent)
        coefficient = multinomial * prod(c ** e for c, e in zip(weights, exponent))
        result += coefficient * prod(factorial(e) for e in exponent)
    return result


def composition_guards():
    records = []
    for weights in [[1, 3, -3], [2, 9, 0], [1, 3, 6, 9]]:
        denominator = [1]
        for weight in weights:
            denominator = mul(denominator, [1, -weight])
        c = sum(weights)
        assert reduce_poly(denominator, 9) == reduce_poly([1, -c], 9)
        hs = complete_homogeneous(weights, 81)
        for degree in range(82):
            assert residue(hs[degree], 9) == pow(c, degree, 9)
        for degree in range(9):
            assert direct_linear_moment(weights, degree) == factorial(degree) * hs[degree]
        for r in range(3):
            q = 3 ** r
            c_q = 3 ** ((q - 1) // 2)
            zeta = 2 ** ((q - 1) // 2) * c ** q
            def units(L):
                return Fraction(factorial(q * L), c_q ** L * factorial(L)) * hs[q * L]
            for n in range(1, 3):
                normalized_kernel(n, units, zeta)
                normalized_kernel(n + 1, units, zeta)
            records.append({"weights": weights, "r": r, "q": q,
                            "c_mod9": c % 9, "effective_zeta_mod9": zeta % 9,
                            "lower_n_checked": [1, 2]})
    dyadic = []
    for weights in [[1, 2, 4], [3, -2, 6]]:
        hs = complete_homogeneous(weights, 48)
        for q in [1, 2, 4]:
            for L in range(13):
                effective = Fraction(factorial(q * L), 2 ** ((q - 1) * L) * factorial(L)) * hs[q * L]
                assert valuation(effective, 2) == 0
            dyadic.append({"weights": weights, "q": q, "kernel_L_checked": [0, 12]})
    return {"ordinary_cubic_composition_cases": records,
            "dyadic_transfer_unit_guards": dyadic,
            "dyadic_all_order_input": "DYADIC_BINOMIAL_SFC.tex, lem:units"}


def gamma_marker_guards():
    count = 0
    for r in [1, 2, 3]:
        q = 3 ** r
        c_q = 3 ** ((q - 1) // 2)
        for delta in range(q):
            epsilon = int(q // 3 <= delta < 2 * q // 3)
            zeta = 2 ** ((q - 1) // 2) * (1 + 3 * epsilon)
            for L in range(7):
                actual = Fraction(factorial(q * L + delta), factorial(delta) * c_q ** L * factorial(L))
                assert valuation(actual, 3) == 0
                assert residue(actual, 9) == pow(zeta, L, 9)
                count += 1
    markers = []
    for q, coefficients in [(3, {0: 2, 1: -1, 2: 1}), (9, {0: 1, 3: 2, 8: -1})]:
        mass = sum(coefficients.values())
        h = sum(Fraction(a, mass) for e, a in coefficients.items() if q // 3 <= e < 2 * q // 3)
        zeta = 1 + 3 * h
        for L in range(9):
            multiplier = sum(Fraction(a * factorial(q * L + e), mass * factorial(e) * factorial(q * L))
                             for e, a in coefficients.items())
            assert valuation(multiplier, 3) == 0
            assert residue(multiplier, 9) == residue(zeta ** L, 9)
        markers.append({"q": q, "divided_power_coefficients": coefficients,
                        "mass": mass, "zeta_mod9": residue(zeta, 9), "L_checked": [0, 8]})
    return {"integer_Gamma_cases": count, "r_checked": [1, 2, 3],
            "all_beta_in_1_through_q_checked": True, "L_checked": [0, 6],
            "marker_controls": markers}


def all_prime_guards():
    residuals = []
    for p in [3, 5, 7, 11, 13, 17, 19]:
        n = 2 * p - 2
        factor = power([1, 1], p - 2)
        a = add([1, 1], [0] * p + [1])
        b = add([-1, -1], [0] * (p + 1) + [1])
        assert reduce_poly(model(n, p), p) == reduce_poly(mul(factor, a), p)
        assert reduce_poly(model(n + 1, p), p) == reduce_poly(mul(factor, b), p)
        assert gcd_mod(model(n, p), model(n + 1, p), p) == reduce_poly(factor, p)
        assert sub(b, mul([0, 1], a)) == scale(power([1, 1], 2), -1)
        residuals.append({"p": p, "n": n, "common_residual_multiplicity_at_minus_one": p - 2})

    normalization_count = 0
    for p in [3, 5, 7, 11, 13]:
        for r in ([1, 2] if p <= 7 else [1]):
            a = p ** r
            for n in range(1, 4):
                for k in range(n + 1):
                    raw = comb(n, k) * factorial(a * (n + (p - 1) * k)) // factorial(a * n)
                    actual, rem = divmod(raw, p ** (a * k))
                    assert rem == 0
                    unit = Fraction(unit_block(p, k))
                    for j in range(r):
                        unit *= Fraction(unit_block(p, p ** j * (n + (p - 1) * k)),
                                         unit_block(p, p ** j * n))
                    assert actual == comb(n + (p - 1) * k, p * k) * unit
                    assert actual % (p * p) == ((-1) ** k * comb(n + (p - 1) * k, p * k)) % (p * p)
                    normalization_count += 1
    return {"uniform_residual_formula_prime_guards": residuals,
            "actual_adjacent_prime_power_normalization_coefficient_guards": normalization_count}


def exact_weighted_moment(n, gap, phase_parameter, units):
    return sum((Fraction(comb(n, k) * factorial(n + gap * k), factorial(n))
                * units(n + gap * k) * phase_parameter ** k for k in range(n + 1)), Fraction())


def rational_record(value):
    value = Fraction(value)
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def counterexamples():
    modified = {3: Fraction(220), 4: Fraction(-12860)}
    phase = Fraction(-1, 6)
    baseline = [exact_weighted_moment(n, 2, phase, lambda L: Fraction(1)) for n in [3, 4]]
    assert baseline == [-219, 12861]
    normalized_bad = []
    error_depths = []
    for n in [3, 4]:
        assert residue(modified[n], 3) == 1
        assert exact_weighted_moment(n, 2, phase, lambda L: modified.get(L, Fraction(1))) == 0
        coefficients = [Fraction(0)] * (n + 1)
        for k in range(n + 1):
            coefficients[n-k] = (Fraction(comb(n, k) * factorial(n + 2*k), factorial(n))
                                  * modified.get(n + 2*k, Fraction(1))
                                  / modified[n] / Fraction(-6) ** k)
        assert coefficients[-1] == 1 and valuation(coefficients[0], 3) == 0
        assert evaluate(coefficients, 1) == 0
        errors = sub(coefficients, model(n))
        error_depths.append(min(valuation(c, 3) for c in errors if c))
        normalized_bad.append(coefficients)
    assert min(error_depths) == 1
    assert resultant(*normalized_bad) == 0

    x = 8
    values = [evaluate(model(n, 5), x) % 125 for n in [156, 157]]
    assert values == [25, 100]
    phase5 = Fraction(-1, 960)
    baseline5 = {n: exact_weighted_moment(n, 4, phase5, lambda L: Fraction(1)) for n in [156, 157]}
    units5 = {n: 1 - value for n, value in baseline5.items()}
    for n in [156, 157]:
        assert residue(units5[n], 25) == 1
        assert valuation(units5[n], 5) == 0
        assert exact_weighted_moment(n, 4, phase5, lambda L: units5.get(L, Fraction(1))) == 0
        assert residue(baseline5[n], 25) == pow(x, -n, 25) * evaluate(model(n, 5), x) % 25
    return {
        "mod3_is_insufficient": {"u3": 220, "u4": -12860, "all_other_units": 1,
                                "phase": "X*(1-X^2/6)", "orders": [3, 4],
                                "baseline_normalized_moments": [int(value) for value in baseline],
                                "modified_normalized_moments": [0, 0],
                                "normalized_common_root": 1,
                                "minimum_coefficient_error_valuations": error_depths,
                                "normalized_resultant": 0},
        "all_prime_mod_p_squared_shortcut_is_false": {
            "p": 5, "orders": [156, 157], "model_point": 8,
            "model_values_mod125": values, "phase": "X*(1-X^4/960)",
            "all_other_units": 1, "modified_units": {str(n): rational_record(value) for n, value in units5.items()},
            "modified_units_mod25": [1, 1], "modified_normalized_moments": [0, 0],
            "ordinary_SFC_counterexample": False}}


def run(repo_root):
    pins = verify_reviewed_sources(repo_root)
    models = model_and_determinant_guards()
    kernels = kernel_guards()
    compositions_report = composition_guards()
    shapes = gamma_marker_guards()
    primes = all_prime_guards()
    controls = counterexamples()
    return {"status": "PASS", "base_commit": "19eb9912fb60142f8b620dac24896931d4ff3e46",
            "python_version": platform.python_version(), "dependencies": "Python standard library only",
            "scope": "exact finite determinant and kernel controls, universal-identity guards, and explicit counterexamples",
            "unrestricted_quantifiers_proved_by_finite_sampling": False,
            "universal_proof": "SFC_CUBIC_UNIT_DEFORMATION_TRANSFER.md, Sections 1-6",
            "source_pins": pins, "model_and_resultant_guards": models,
            "unit_kernel_guards": kernels, "ordinary_composition_guards": compositions_report,
            "Gamma_and_marker_guards": shapes, "odd_prime_guards": primes,
            "counterexamples": controls,
            "assurance": "internal exact replay; written proof supplies all-parameter statements"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = run(args.repo_root.resolve())
    output = args.output or Path(__file__).with_name("sfc_cubic_unit_deformation_verification.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "model_resultant_pairs": 8,
                      "raw_unit_kernel_resultant_pairs": 16,
                      "ordinary_cubic_composition_cases": 9,
                      "counterexamples_verified": 2, "receipt": str(output)}))
