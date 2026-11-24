#!/usr/bin/env python3
"""Exact finite replay of factorial-kernel density and its rank obstructions.

Only Python's standard library is used. Large constructed factorials are
never formed. The unrestricted statements have proofs in the pinned note;
the receipt records exactly the finite arithmetic replay performed here.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

if not __debug__:
    raise SystemExit("Refusing optimized Python: this checker uses assertions.")

NOTE_SHA256 = "df9f19e4418c2abdf47c0f2696a60637dcfb761f93a9b91457aad4b5c671e6b0"
HERE = Path(__file__).resolve().parent


def prime(p: int) -> bool:
    return p >= 2 and all(p % d for d in range(2, math.isqrt(p) + 1))


def vp(n: int, p: int) -> int | None:
    if n == 0:
        return None
    n = abs(n)
    out = 0
    while n % p == 0:
        out += 1
        n //= p
    return out


def digit_sum(n: int, p: int) -> int:
    out = 0
    while n:
        n, a = divmod(n, p)
        out += a
    return out


def factorial_valuation(n: int, p: int) -> int:
    out = 0
    while n:
        n //= p
        out += n
    return out


class KernelMod:
    def __init__(self, p: int, precision: int):
        assert prime(p) and precision >= 1
        self.p = p
        self.precision = precision
        self.modulus = p**precision
        self.prefix = [1] * self.modulus
        for j in range(1, self.modulus):
            self.prefix[j] = self.prefix[j - 1]
            if j % p:
                self.prefix[j] = self.prefix[j] * j % self.modulus
        self.block = self.prefix[-1]

    def W(self, n: int) -> int:
        assert n >= 0
        blocks, rest = divmod(n, self.modulus)
        return pow(self.block, blocks, self.modulus) * self.prefix[rest] % self.modulus

    def unit_factorial(self, n: int) -> int:
        out = 1
        while n:
            out = out * self.W(n) % self.modulus
            n //= self.p
        assert math.gcd(out, self.p) == 1
        return out

    def kernel(self, d: int, k: int) -> tuple[int, int]:
        assert d >= 0 and k >= 0
        value = factorial_valuation(d * k, self.p) - k * factorial_valuation(d, self.p)
        digit_value = (k * digit_sum(d, self.p) - digit_sum(d * k, self.p)) // (self.p - 1)
        assert value == digit_value and value >= 0
        if value >= self.precision:
            return value, 0
        unit = self.unit_factorial(d * k)
        denominator = pow(self.unit_factorial(d), k, self.modulus)
        unit = unit * pow(denominator, -1, self.modulus) % self.modulus
        return value, pow(self.p, value, self.modulus) * unit % self.modulus

    def R(self, n: int, k: int) -> int:
        denominator = pow(self.W(n), k, self.modulus)
        return self.W(k * n) * pow(denominator, -1, self.modulus) % self.modulus


def difference_valuation(residue: int, p: int, precision: int) -> dict:
    if residue == 1:
        return {"kind": "lower_bound", "value": precision}
    value = vp(residue - 1, p)
    assert value is not None and value < precision
    return {"kind": "exact", "value": value}


def constructed_d(p: int, r: int, count: int | None = None, spacing: int | None = None) -> int:
    if count is None:
        count = (p - 1) * p ** (r - 1)
    if spacing is None:
        spacing = r
    return sum(p ** (spacing * i) for i in range(1, count + 1))


def compositions(total: int, slots: int):
    if slots == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, slots - 1):
            yield (first,) + rest


def multinomial(k: tuple[int, ...]) -> int:
    return math.factorial(sum(k)) // math.prod(math.factorial(j) for j in k)


def base_coefficient(rows: list[list[int]], k: tuple[int, ...]) -> int:
    return multinomial(k) * math.prod(
        math.factorial(sum(a * b for a, b in zip(row, k))) for row in rows
    )


def determinant(matrix: list[list[int]]) -> int:
    n = len(matrix)
    assert all(len(row) == n for row in matrix)
    if n == 0:
        return 1
    return sum(
        (-1) ** j * matrix[0][j]
        * determinant([row[:j] + row[j + 1:] for row in matrix[1:]])
        for j in range(n)
    )


QZERO = (Fraction(0), Fraction(0))
QONE = (Fraction(1), Fraction(0))
QX = (Fraction(0), Fraction(1))


def qadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def qmul(a, b):
    # x^2 = -(4/3)x - 1, exactly.
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0] - Fraction(4, 3) * a[1] * b[1])


def qpow(a, n: int):
    out = QONE
    for _ in range(n):
        out = qmul(out, a)
    return out


def qscale(a, value):
    return (a[0] * value, a[1] * value)


def qjson(a):
    return [str(v) for v in a]


def qmod(a, root: int, p: int) -> int:
    return sum(v.numerator * pow(v.denominator, -1, p) * root**i
               for i, v in enumerate(a)) % p


def qmoment(rows, coefficients, order: int, kernels):
    out = QZERO
    for k in compositions(order, len(coefficients)):
        term = QONE
        for coefficient, power in zip(coefficients, k):
            term = qmul(term, qpow(coefficient, power))
        scalar = base_coefficient(rows, k) * math.prod(kernels[j] for j in k)
        out = qadd(out, qscale(term, scalar))
    return out


def minor_minimum(matrix: list[list[int]], p: int) -> int | None:
    width = len(matrix[0])
    values = []
    for indices in itertools.combinations(range(len(matrix)), width):
        value = vp(determinant([matrix[i] for i in indices]), p)
        if value is not None:
            values.append(value)
    return min(values) if values else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=HERE / "SFC_DEAMPLIFICATION_BRIDGES.md")
    parser.add_argument("--output", type=Path, default=HERE / "sfc_deamplification_arithmetic_verification.json")
    args = parser.parse_args()
    actual_hash = hashlib.sha256(args.note.read_bytes()).hexdigest()
    if actual_hash != NOTE_SHA256:
        raise SystemExit(f"Pinned note mismatch: expected {NOTE_SHA256}, got {actual_hash}")

    receipt = {
        "status": "PASS",
        "scope": "Finite exact replay supporting the separately written all-parameter proofs; no SFC de-amplification claim.",
        "source_note": str(args.note),
        "source_note_sha256": actual_hash,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "algorithm": "Legendre valuation and recursive factorial-unit products modulo prime powers; constructed large factorials are never formed.",
    }

    # Independent small-input check against actual integer factorial quotients.
    small_count = 0
    for p in (2, 3, 5, 7):
        for t in range(1, 5):
            arithmetic = KernelMod(p, t)
            for d in range(0, 25):
                for k in range(0, 11):
                    exact = math.factorial(d * k) // math.factorial(d) ** k
                    value, residue = arithmetic.kernel(d, k)
                    assert value == vp(exact, p)
                    assert residue == exact % arithmetic.modulus
                    small_count += 1
            for n in range(0, 32):
                for k in range(0, p + 1):
                    assert arithmetic.R(n, k) == arithmetic.R(n + arithmetic.modulus, k)
    receipt["direct_factorial_crosscheck"] = {
        "primes": [2, 3, 5, 7], "precisions": [1, 2, 3, 4],
        "D_range": [0, 24], "k_range": [0, 10], "cases": small_count,
        "largest_factorial_argument": 240,
    }

    # Exact digit-window construction at several finite precisions.
    digit_cases = []
    digit_specs = [(2, r) for r in range(1, 5)]
    digit_specs += [(3, r) for r in range(1, 5)]
    digit_specs += [(5, r) for r in range(1, 4)]
    digit_specs += [(7, r) for r in range(1, 3)] + [(11, 1)]
    for p, r in digit_specs:
        d = constructed_d(p, r)
        count = (p - 1) * p ** (r - 1)
        arithmetic = KernelMod(p, r)
        higher = KernelMod(p, r + 2)
        rows = []
        for k in range(p):
            value, residue = arithmetic.kernel(d, k)
            assert value == 0 and residue == 1
            quot = d
            modulus = 1
            windows = []
            while quot:
                assert (k * d) // modulus == k * (d // modulus)
                windows.append(arithmetic.R(quot, k))
                modulus *= p
                quot //= p
            product = math.prod(windows) % arithmetic.modulus
            components = [arithmetic.R(p**j, k) for j in range(r)]
            cycle = math.prod(components) % arithmetic.modulus
            assert product == residue == pow(cycle, count, arithmetic.modulus)
            _, lifted = higher.kernel(d, k)
            difference = difference_valuation(lifted, p, r + 2)
            assert difference["value"] >= r
            rows.append({
                "k": k, "kernel_valuation": value,
                "kernel_mod_p_to_r": residue,
                "one_digit_window_factors": components,
                "window_cycle_mod_p_to_r": cycle,
                "kernel_mod_p_to_r_plus_2": lifted,
                "valuation_of_kernel_minus_one": difference,
            })
        digit_cases.append({
            "p": p, "r": r, "M": p - 1, "L": count,
            "D": str(d), "nonzero_base_p_positions": list(range(r, r * count + 1, r)),
            "rows": rows,
        })
    receipt["digit_window_cases"] = digit_cases

    bad_spacing = constructed_d(3, 2, spacing=1)
    bad_count = constructed_d(3, 2, count=2)
    mod9 = KernelMod(3, 2)
    assert bad_spacing == 1092 and mod9.kernel(bad_spacing, 2) == (0, 7)
    assert bad_count == 90 and mod9.kernel(bad_count, 2) == (0, 4)
    receipt["negative_digit_controls"] = [
        {"p": 3, "r": 2, "spacing": 1, "L": 6, "D": 1092, "k": 2, "kernel_mod_9": 7},
        {"p": 3, "r": 2, "spacing": 2, "L": 2, "D": 90, "k": 2, "kernel_mod_9": 4},
    ]

    # The necessary M<p condition and two tempting but wrong sequences.
    prime_boundary = []
    natural_sequences = []
    for p in (2, 3, 5, 7):
        for d in range(1, 65):
            value = factorial_valuation(p * d, p) - p * factorial_valuation(d, p)
            assert value == digit_sum(d, p) and value > 0
        prime_boundary.append({"p": p, "D_range": [1, 64], "claim": "v_p K_D(p)=s_p(D)>0"})
        for h in range(1, 5):
            power = p**h
            arithmetic = KernelMod(p, h + 2)
            for k in range(1, min(power, 10) + 1):
                value, residue = arithmetic.kernel(power, k)
                assert value == factorial_valuation(k, p)
                if k < p:
                    assert residue % p == math.factorial(k) % p
                complement_value, _ = arithmetic.kernel(power - 1, k)
                assert complement_value == h * (k - 1)
            natural_sequences.append({
                "p": p, "h": h, "k_range": [1, min(power, 10)],
                "power_sequence_valuation": "v_p(k!)",
                "complement_sequence_valuation": "h(k-1)",
            })
    receipt["prime_boundary_controls"] = prime_boundary
    receipt["natural_sequence_controls"] = natural_sequences

    # Exact moment-coefficient congruences, including an unbalanced fixed base.
    distinct_rows = [[1, 1, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 1]]
    arbitrary_rows = [[1, 0], [0, 2]]
    coefficient_cases = []
    for name, support, p in (
            ("distinct_balanced_YZ_YW_UV", distinct_rows, 7),
            ("unbalanced_X_Y_squared", arbitrary_rows, 5)):
        for r in (1, 2):
            arithmetic = KernelMod(p, r)
            d = constructed_d(p, r)
            values = {k: arithmetic.kernel(d, k)[1] for k in range(4)}
            witnesses = []
            for m in range(1, 4):
                for k in compositions(m, len(support[0])):
                    original = base_coefficient(support, k)
                    decorated = original * math.prod(values[j] for j in k) % arithmetic.modulus
                    assert decorated == original % arithmetic.modulus
                    witnesses.append({"order": m, "multi_index": list(k), "original_coefficient": original,
                                      "decorated_coefficient_modulus": decorated})
            coefficient_cases.append({"base": name, "rows": support, "p": p, "r": r,
                                      "D": str(d), "coefficients": witnesses})
    receipt["actual_moment_coefficient_congruences"] = coefficient_cases

    # The full degree-two Macaulay matrix is retained, not just its determinant.
    determinant_cases = []
    for d in range(0, 13):
        k = math.comb(2 * d, d)
        matrix = [[1, 1, 0], [0, 1, 1], [k, 2, k]]
        delta = determinant(matrix)
        assert delta == 2 * (k - 1)
        assert (delta == 0) == (d == 0)
        determinant_cases.append({"D": d, "matrix": matrix, "determinant": delta})
    determinant_precision = []
    for r in range(1, 5):
        d = constructed_d(3, r)
        arithmetic = KernelMod(3, r + 2)
        _, residue = arithmetic.kernel(d, 2)
        delta = 2 * (residue - 1) % arithmetic.modulus
        value = vp(delta, 3)
        assert value is None or value >= r
        determinant_precision.append({
            "p": 3, "r": r, "D": str(d), "determinant_mod_p_to_r_plus_2": delta,
            "valuation": {"kind": "lower_bound", "value": r + 2} if value is None
                         else {"kind": "exact", "value": value},
        })
    assert determinant_cases[12]["determinant"] == 5408310
    receipt["two_column_full_determinant_obstruction"] = {
        "qualification": "Unamplified base columns coincide; this is not an SFC counterexample.",
        "small_exact_cases": determinant_cases,
        "constructed_precision_cases": determinant_precision,
    }

    # Distinct balanced base columns: two vanishing original moments but not three.
    coefficients = [QONE, QX, (Fraction(-1), Fraction(-1))]
    assert qadd(qadd(coefficients[0], coefficients[1]), coefficients[2]) == QZERO
    assert qadd(qadd(qpow(coefficients[0], 2), qpow(coefficients[1], 2)),
                qpow(coefficients[2], 2)) == (Fraction(0), Fraction(-2, 3))
    original_moments = [qmoment(distinct_rows, coefficients, m, [1] * (m + 1)) for m in range(1, 4)]
    assert original_moments == [QZERO, QZERO, (Fraction(60), Fraction(20))]
    decorated_moments = []
    for d in range(1, 13):
        kernels = [math.factorial(d * k) // math.factorial(d) ** k for k in range(4)]
        first = qmoment(distinct_rows, coefficients, 1, kernels)
        second = qmoment(distinct_rows, coefficients, 2, kernels)
        assert first == QZERO
        assert second == (Fraction(0), -Fraction(8, 3) * (kernels[2] - 1))
        assert second != QZERO
        decorated_moments.append({"D": d, "second_moment": qjson(second)})
    residual = []
    for root in (3, 5):
        assert (3 * root**2 + 4 * root + 3) % 7 == 0
        assert (6 * root + 4) % 7 != 0
        assert root != 0
        third = qmod(original_moments[2], root, 7)
        assert third in (1, 6)
        assert qmod((Fraction(0), Fraction(-8, 3)), root, 7) != 0
        residual.append({"p": 7, "root": root, "third_moment_mod_7": third})
    receipt["distinct_support_prefix_obstruction"] = {
        "base_rows": distinct_rows,
        "coefficient_field": "Q[x]/(3x^2+4x+3)",
        "coefficients_in_basis_1_x": [qjson(c) for c in coefficients],
        "original_moments_orders_1_2_3": [qjson(c) for c in original_moments],
        "decorated_second_moment_formula": "-(8x/3)(K_D(2)-1)",
        "small_actual_decorated_cases": decorated_moments,
        "residue_field_controls": residual,
        "qualification": "Only the prefix of length two is poorly conditioned; the support-size three window has a unit third moment at both places above seven.",
    }

    # Exact closed pure-coefficient invariant for balanced profiles.
    pure_controls = []
    for h1, hm, b1, bm, lamj, lamh, m in (
            (1, 2, 1, 1, 1, 2, 2), (2, 24, 3, 5, 7, 11, 2),
            (6, 362880, 2, 7, 3, 5, 3), (1, 1, 0, 3, 2, 4, 4)):
        a1j, a1h = b1 * h1 * lamj, b1 * h1 * lamh
        amj, amh = bm * hm * lamj**m, bm * hm * lamh**m
        value = amj * a1h**m - amh * a1j**m
        assert value == 0
        pure_controls.append({"m": m, "pure_order_1": [a1j, a1h],
                              "pure_order_m": [amj, amh], "invariant": value})
    target_invariant = 2 * 2**2 - 24 * 1**2
    assert target_invariant == -16
    receipt["balanced_limit_invariant"] = {
        "formula": "A_mj*A_1h^m-A_mh*A_1j^m",
        "balanced_controls": pure_controls,
        "unbalanced_target": {"columns": ["X", "Y^2"], "pure_order_1": [1, 2],
                              "pure_order_2": [2, 24], "invariant": target_invariant},
        "scope": "Rules out balanced coefficientwise limits under diagonal coefficient rescaling and degreewise scalar normalization; does not address arbitrary coefficient mixing.",
    }

    # Representative exact matrix controls for the general polynomial-in-entries lemma.
    matrix_controls = []
    full_limit = [[3, 0], [0, 1], [0, 0]]
    deficient_limit = [[1, 2], [2, 4], [0, 0]]
    assert minor_minimum(full_limit, 3) == 1
    assert minor_minimum(deficient_limit, 3) is None
    for r in range(2, 6):
        p_to_r = 3**r
        full_approx = [[3 + p_to_r, p_to_r], [p_to_r, 1 + p_to_r], [p_to_r, 0]]
        deficient_approx = [[1, 2], [2, 4 + p_to_r], [0, 0]]
        assert minor_minimum(full_approx, 3) == 1
        assert minor_minimum(deficient_approx, 3) == r
        matrix_controls.append({"p": 3, "r": r, "full_limit_approximation": full_approx,
                                "full_limit_minimum": 1, "deficient_limit_approximation": deficient_approx,
                                "deficient_limit_minimum": r})
    receipt["rank_conditioning_controls"] = {
        "full_limit": full_limit, "deficient_limit": deficient_limit,
        "cases": matrix_controls,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"PASS: {small_count} direct factorial cross-checks; {len(digit_cases)} digit-window constructions; exact determinant, distinct-prefix, balanced-limit and rank controls.")
    print(f"Receipt: {args.output}")


if __name__ == "__main__":
    main()
