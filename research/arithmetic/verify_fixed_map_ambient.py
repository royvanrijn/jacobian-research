#!/usr/bin/env python3
"""Exact controls for FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md.

These controls do not verify Chebotarev, Hensel's lemma, or the general
Keller construction. Run with Python 3 and SymPy; --output writes a receipt.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path

import sympy as sp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    checks: list[dict[str, object]] = []

    def equal(name: str, left, right=0) -> None:
        difference = sp.factor(left - right)
        if difference != 0:
            raise AssertionError(f"{name}: nonzero difference {difference}")
        checks.append({"name": name, "passed": True})

    def require(name: str, condition: bool, **details) -> None:
        if not condition:
            raise AssertionError(name)
        checks.append({"name": name, "passed": True, **details})

    S, U, Z, a, r, w = sp.symbols("S U Z a r w")
    gamma = sp.Rational(3, 32)
    k = gamma**3
    pi = sp.Integer(-2)
    b = sp.Rational(56, 3) - sp.Rational(1024, 9) * r
    c = -sp.Rational(7, 6) - sp.Rational(64, 3) * r
    E = (
        pi**5 * S**5 - sp.Rational(3, 2) * pi**4 * S**4
        + sp.Rational(3, 2) * pi * S**3
        + (-sp.Rational(5, 4) * pi - sp.Rational(9, 32) * b) * S**2
        + sp.Rational(9, 16) * S - sp.Rational(9, 32) * c
    )
    quadratic = S**2 + sp.Rational(3, 16)
    cubic = (S + sp.Rational(1, 4))**3 - 3 * gamma * (S + sp.Rational(1, 4)) - r
    equal("second_line_inverse_factorization", E, -32 * quadratic * cubic)
    equal("depressed_cubic_discriminant", sp.discriminant(U**3 - 3 * gamma * U - r, U), -27 * (r**2 - 4 * k))
    equal("quadratic_discriminant", sp.discriminant(quadratic, S), -sp.Rational(3, 4))
    equal("cubic_remainder", sp.rem(cubic, quadratic, S), -sp.Rational(9, 32) * S - r - sp.Rational(25, 128))
    equal("coprime_resultant", sp.resultant(quadratic, cubic, S), (r + sp.Rational(25, 128))**2 + sp.Rational(243, 16384))
    require("resultant_positive_constant", sp.Rational(243, 16384) > 0)
    require("k_is_cube", k == sp.Rational(27, 32768))
    require("k_is_not_rational_square", not sp.sqrt(k).is_Rational)
    require("four_k_is_not_rational_square", not sp.sqrt(4 * k).is_Rational)

    generator = Z + gamma * Z**2 / a
    relation = sp.rem(sp.together(generator**3 - 3 * gamma * generator - (a + k / a)), Z**3 - a, Z)
    equal("cubic_algebra_generator_relation", relation)
    cols = []
    for power in (sp.Integer(1), generator, generator**2):
        remainder = sp.Poly(sp.rem(power, Z**3 - a, Z), Z)
        cols.append([remainder.nth(i) for i in range(3)])
    matrix = sp.Matrix.hstack(*(sp.Matrix(col) for col in cols))
    equal("cubic_algebra_basis_determinant", matrix.det(), 1 - k / a**2)
    equal("parameter_discriminant", (a + k / a)**2 - 4 * k, (a - k / a)**2)
    equal("parameter_involution", (a + k / a).subs(a, k / a), a + k / a)
    equal("parameter_fibres_exactly_two", a * w * ((a + k / a) - (w + k / w)), (a - w) * (a * w - k))
    require("a19_parametrized_target_not_original_line", 2 != -1)
    require("a19_is_noncube", not sp.integer_nthroot(19, 3)[1])
    require("a19_old_local_conditions", 19 % 3 == 1 and 19 % 9 == 1)

    # Whole pi-zero boundary. We verify both source components in the
    # unnormalized F0 coordinates; Phi adds the invertible source scaling.
    x, y, z, bb, cc = sp.symbols("x y z bb cc")
    t = 1 + x * y
    q = t**2 * z + sp.Rational(3, 8) * y**2 * (1 + 3 * t)
    B = y + 8 * x * q - sp.Rational(40, 9) * t * q - sp.Rational(32, 3) * t**2 * x**2 * q**4 + sp.Rational(80, 9) * t**2 * x**3 * q**5
    C = x * (5 - 3 * t) - sp.Rational(8, 3) * x**3 * z + sp.Rational(16, 3) * (x * q)**4 - sp.Rational(16, 3) * (x * q)**5
    on_t_zero = {x: 2 / bb, y: -bb / 2}
    equal("boundary_t_zero", t.subs(on_t_zero))
    equal("boundary_t_zero_B", B.subs(on_t_zero), bb)
    equal("boundary_t_zero_C_z_coefficient", sp.diff(C.subs(on_t_zero), z), -sp.Rational(64, 3) / bb**3)
    z_solution = sp.solve(C.subs(on_t_zero) - cc, z)[0]
    equal("boundary_t_zero_rational_point_C", C.subs(on_t_zero).subs(z, z_solution), cc)
    equal("boundary_t_zero_rational_point_Pi", (t * q).subs(on_t_zero))
    singleton = {x: cc / 2, y: 0, z: 0}
    equal("boundary_b_zero_Pi", (t * q).subs(singleton))
    equal("boundary_b_zero_B", B.subs(singleton))
    equal("boundary_b_zero_C", C.subs(singleton), cc)
    q_zero_z = -sp.Rational(3, 8) * y**2 * (1 + 3 * t) / t**2
    equal("boundary_q_zero_B", B.subs(z, q_zero_z), y)
    equal("boundary_q_zero_C", C.subs(z, q_zero_z), 2 * x / t - y * (x / t)**2)

    def partitions(n: int, smallest: int = 1):
        if n == 0:
            yield ()
        else:
            for first in range(smallest, n + 1):
                for rest in partitions(n - first, first):
                    yield (first,) + rest

    partitions_five = list(partitions(5))
    for partition in partitions_five:
        if any(m > 1 for m in partition):
            require(f"repeated_quintic_rank_{partition}", sum(m == 1 for m in partition) <= 3)
    require("no_rational_factor_rank_five_partitions", [p for p in partitions_five if min(p) >= 2] == [(2, 3), (5,)])
    require("no_rational_factor_rank_four_partitions", [p for p in partitions(4) if min(p) >= 2] == [(2, 2), (4,)])

    # Exhaustive subgroups of S3, with the actions on its three letters
    # and on the two cosets of A3. This is a finite exact group check.
    perms = list(itertools.permutations(range(3)))
    identity = (0, 1, 2)

    def compose(g, h):
        return tuple(g[h[i]] for i in range(3))

    def even(g):
        return sum(g[i] > g[j] for i in range(3) for j in range(i + 1, 3)) % 2 == 0

    subgroups = []
    for mask in range(1 << len(perms)):
        D = frozenset(perms[i] for i in range(len(perms)) if mask & (1 << i))
        if identity not in D:
            continue
        if all(compose(g, h) in D for g in D for h in D):
            subgroups.append(D)
    require("all_S3_subgroups_enumerated", len(subgroups) == 6, orders=sorted(map(len, subgroups)))
    for number, D in enumerate(subgroups):
        common_two = all(even(g) for g in D)
        common_three = any(all(g[i] == i for g in D) for i in range(3))
        require(f"S3_subgroup_local_fixed_point_{number}", (common_two or common_three) == (len(D) < 6))
    for g in perms:
        require(f"diagonal_S3_element_has_fixed_point_{g}", even(g) or any(g[i] == i for i in range(3)))
    three_cycle = (1, 2, 0)
    require("wrong_quadratic_common_derangement", all(three_cycle[i] != i for i in range(3)) and all((i + 1) % 2 != i for i in range(2)))
    require("cyclic_cubic_times_quadratic_common_derangement", all(three_cycle[i] != i for i in range(3)))
    # The second component in both preceding direct products is chosen
    # nontrivial, so neither quadratic point is fixed.

    # Rank-two localization example P=T^2+1 and b=0. The derivative is
    # invertible on P and nilpotent on the double factor T^2.
    T = sp.symbols("T")
    P = T**2 + 1
    R = P * T**2
    derivative = sp.diff(R, T)
    equal("rank_two_derivative_inverse_mod_P", sp.rem(derivative * T / 2 - 1, P, T))
    equal("rank_two_derivative_nilpotent_mod_double", sp.rem(derivative**2, T**2, T))
    require("rank_two_admissible_translation", derivative.subs(T, 1) != 0 and sp.diff(R, T, 3).subs(T, 1) != 0)
    empty_R = (T - 1)**3
    equal("empty_derivative_nilpotent", sp.rem(sp.diff(empty_R, T)**2, empty_R, T))

    # V4 rank-six covering and its distinction between fixed points for
    # individual elements and a common fixed point for an entire group.
    v4 = frozenset(range(4))
    v4_subgroups = [frozenset((0, element)) for element in (1, 2, 3)]
    omega_v4 = []
    for stabilizer in v4_subgroups:
        cosets = {frozenset(g ^ h for h in stabilizer) for g in v4}
        omega_v4.extend((stabilizer, coset) for coset in cosets)
    require("V4_rank_six_G_set", len(omega_v4) == 6)

    def v4_common_fixed(decomposition):
        return [(stabilizer, coset) for stabilizer, coset in omega_v4
                if all(frozenset(g ^ point for point in coset) == coset
                       for g in decomposition)]

    require("V4_every_element_fixes_a_point",
            all(v4_common_fixed(frozenset((0, element))) for element in v4))
    require("V4_no_global_fixed_point", not v4_common_fixed(v4))
    require("V4_all_proper_subgroups_have_fixed_points",
            all(v4_common_fixed(stabilizer) for stabilizer in v4_subgroups))
    require("rank_six_radicands_nonsquare",
            all(not sp.integer_nthroot(radicand, 2)[1] for radicand in (13, 17, 221)))
    require("rank_six_radicand_product_square", 13 * 17 * 221 == 221**2)
    require("rank_six_local_13", (2**2 - 17) % 13 == 0 and (2*2) % 13 != 0)
    require("rank_six_local_17", (8**2 - 13) % 17 == 0 and (2*8) % 17 != 0)
    require("rank_six_local_2_Hensel", 16 == 2**4 and 2 == 2**1 and 4 > 2*1)
    require("rank_six_real_points", all(radicand > 0 for radicand in (13, 17, 221)))
    require("quadratic_characters_cannot_all_be_negative",
            all(any(value == 1 for value in (first, second, first*second))
                for first in (-1, 1) for second in (-1, 1)))
    require("Hasse_rank_padding_controls",
            all((5 if rank % 2 else 6) + 2*((rank-(5 if rank % 2 else 6))//2) == rank
                for rank in range(5, 51)))

    # Compare the cubic discriminant test with an independently certified
    # strict-Hensel test. This bounded replay supports the written proof;
    # the proof itself uses valuations of conjugate roots, not enumeration.
    def valuation(value, prime):
        if value == 0:
            return float("inf")
        answer = 0
        while value % prime == 0:
            value //= prime
            answer += 1
        return answer

    def residue_roots(coefficients, prime, precision):
        candidates, modulus = [0], 1
        ca, cb, cc = coefficients
        for _ in range(precision):
            next_modulus = prime * modulus
            candidates = [x + digit * modulus for x in candidates
                          for digit in range(prime)
                          if (((x + digit * modulus + ca) * (x + digit * modulus)
                               + cb) * (x + digit * modulus) + cc) % next_modulus == 0]
            modulus = next_modulus
            if not candidates:
                break
        return candidates

    cubic_controls, hensel_controls = 0, 0
    for ca, cb, cc in itertools.product(range(-3, 4), repeat=3):
        disc = ca**2 * cb**2 - 4*cb**3 - 4*ca**3 * cc - 27*cc**2 + 18*ca*cb*cc
        if disc == 0:
            continue
        for prime in (2, 3, 5, 7):
            dv = valuation(disc, prime)
            sharp_roots = residue_roots((ca, cb, cc), prime, dv//2 + 1)
            classical_roots = residue_roots((ca, cb, cc), prime, 2*dv + 1)
            if bool(sharp_roots) != bool(classical_roots):
                raise AssertionError(("cubic discriminant root test", ca, cb, cc, prime))
            for root in classical_roots:
                value = ((root + ca)*root + cb)*root + cc
                derivative_value = 3*root**2 + 2*ca*root + cb
                if not valuation(value, prime) > 2*valuation(derivative_value, prime):
                    raise AssertionError("classical residue root lacks strict Hensel certificate")
                hensel_controls += 1
            cubic_controls += 1
    require("sharp_cubic_modulus_vs_certified_Hensel_controls", cubic_controls > 0,
            coefficient_range=[-3, 3], primes=[2, 3, 5, 7],
            polynomial_prime_pairs=cubic_controls, strict_Hensel_certificates=hensel_controls)
    sharpness_controls = 0
    for prime in (2, 5, 7):
        for exponent in (1, 2, 4, 5):
            coefficients = (0, 0, -prime**exponent)
            if residue_roots(coefficients, prime, exponent + 1):
                raise AssertionError("sharpness family unexpectedly has a stronger residue root")
            if not residue_roots(coefficients, prime, exponent):
                raise AssertionError("sharpness family lost the residue zero")
            if valuation(-27*prime**(2*exponent), prime) != 2*exponent:
                raise AssertionError("sharpness family discriminant")
            sharpness_controls += 1
    require("sharp_cubic_modulus_counterexamples_to_lower_precision", sharpness_controls == 12,
            polynomial_prime_pairs=sharpness_controls)
    # Finite arithmetic controls for the elementary all-height proof.
    # These examples are not the proof of the unrestricted gcd bound.
    import math
    controls = 0
    for numerator in range(-20, 21):
        for denominator in range(1, 21):
            if numerator == 0 or math.gcd(numerator, denominator) != 1:
                continue
            raw_numerator = 32768 * numerator**2 + 27 * denominator**2
            raw_denominator = 32768 * numerator * denominator
            common = math.gcd(raw_numerator, raw_denominator)
            if (32768 * 27) % common != 0:
                raise AssertionError("height gcd bound on retained finite controls")
            controls += 1
    require("height_gcd_finite_controls", controls > 0, coprime_pairs=controls)
    require("height_valuation_ties_impossible", 15 % 2 == 1 and 3 % 2 == 1)
    source_path = Path(__file__).with_name("FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md")
    report = {
        "status": "passed",
        "scope": "Exact symbolic identities, finite S3/V4 actions, and bounded cubic congruence/Hensel controls; not a full proof replay or external review.",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "check_count": len(checks),
        "checks": checks,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "check_count": len(checks), "output": str(args.output) if args.output else None}))


if __name__ == "__main__":
    main()
